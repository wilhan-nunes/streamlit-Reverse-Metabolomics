import tempfile
import shutil
from pathlib import Path
import requests
import time
import logging
import os
import pandas as pd
from datetime import datetime

# Setup logger
logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def convert_tsv_to_parquet(tsv_path, output_path) -> bool:
    """
    Convert an already-downloaded raw ReDU TSV dump to parquet.

    Args:
        tsv_path: Path to the raw TSV file (e.g. a manually downloaded
            https://redu.gnps2.org/dump dropped in place by the user).
        output_path: Path where the converted parquet file should be written.

    Returns:
        bool: True on success, False on failure.
    """
    tsv_path = Path(tsv_path)
    output_path = Path(output_path)

    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"Converting {tsv_path} to parquet format...")
        df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
        df.to_parquet(output_path, index=False)
        print(f"Conversion completed successfully: {output_path}")
        return True
    except Exception as e:
        print(f"Error converting {tsv_path} to parquet: {e}")
        return False


RAW_REDU_FILENAME = "redu.tsv"


def find_raw_redu_file(output_path):
    """
    Look for a manually-downloaded raw ReDU dump sitting next to output_path.

    Per the manual-download instructions, users must rename the file
    downloaded from https://redu.gnps2.org/dump to RAW_REDU_FILENAME and
    place it in the same directory as the expected parquet file.

    Args:
        output_path: The expected final parquet path.

    Returns:
        Path to the raw file if found, else None.
    """
    output_path = Path(output_path)
    raw_path = output_path.parent / RAW_REDU_FILENAME
    return raw_path if raw_path.is_file() else None


def download_redu_metadata(output_path):
    """
    Download ReDU metadata file from GNPS2 using streaming and save as parquet.
    Downloads to a temporary file first, then converts to parquet and moves to final location.

    Args:
        output_path: Path where to save the downloaded parquet file
    """
    url = "https://redu.gnps2.org/dump"
    output_path = Path(output_path)

    # Ensure output has .parquet extension
    if output_path.suffix != ".parquet":
        output_path = output_path.with_suffix(".parquet")

    # Create directory if it doesn't exist
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Create temporary file in the same directory as the output
    temp_file = None
    try:
        # Create temporary file in the same directory as the final output
        with tempfile.NamedTemporaryFile(
            dir=output_path.parent, delete=False, suffix=".tmp"
        ) as temp_file:
            temp_path = Path(temp_file.name)

            # Make request with stream enabled
            with requests.get(url, stream=True) as response:
                response.raise_for_status()

                # Get total file size if available
                total_size = int(response.headers.get("content-length", 0))

                # Stream the download to temporary file
                if total_size == 0:
                    temp_file.write(response.content)
                else:
                    downloaded = 0
                    start_time = time.time()
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            temp_file.write(chunk)
                            downloaded += len(chunk)

                            # Calculate progress
                            progress = (downloaded / total_size) * 100
                            elapsed_time = time.time() - start_time
                            speed = downloaded / (1024 * 1024 * elapsed_time)  # MB/s

                            print(
                                f"\rDownloading: {progress:.1f}% ({speed:.1f} MB/s)",
                                end="",
                            )
                    print()  # New line after download completes

        # Read TSV and convert to parquet
        success = convert_tsv_to_parquet(temp_path, output_path)

        # Clean up temporary TSV file
        temp_path.unlink()

        if not success:
            return False

        print(f"Download completed successfully: {output_path}")
        return True

    except requests.RequestException as e:
        print(f"Error downloading file: {e}")
        # Clean up temporary file if it exists
        if temp_file and Path(temp_file.name).exists():
            Path(temp_file.name).unlink()
        return False
    except Exception as e:
        print(f"Unexpected error: {e}")
        # Clean up temporary file if it exists
        if temp_file and Path(temp_file.name).exists():
            Path(temp_file.name).unlink()
        return False


def get_redu_path() -> Path:
    """Resolve where the ReDU metadata file lives.

    Server deployments get the historical default (the working directory).
    Local runs set REVMET_DATA_DIR so the download persists outside the repo.
    REDU_METADATA_PATH overrides both, for a manually downloaded file.
    """
    explicit = os.environ.get("REDU_METADATA_PATH")
    if explicit:
        return Path(explicit).expanduser()

    data_dir = os.environ.get("REVMET_DATA_DIR")
    if data_dir:
        return Path(data_dir).expanduser() / "REDU_metadata.parquet"

    return Path("REDU_metadata.parquet")


def manual_download_hint(file_path) -> str:
    """Instructions for users whose automatic download failed.

    Paths are wrapped in backticks: this string is rendered as Markdown
    (e.g. via st.error), which otherwise treats a bare backslash before
    punctuation as an escape and silently drops it — corrupting Windows
    paths like "C:\\Users\\name\\.reverse-metabolomics".
    """
    file_path = Path(file_path)
    return (
        f"Download https://redu.gnps2.org/dump manually, rename the downloaded "
        f"file to `{RAW_REDU_FILENAME}`, and place it in `{file_path.parent}` — "
        f"the app will detect and convert it automatically next time it loads. "
        f"Alternatively, convert it to parquet yourself and save it to "
        f"`{file_path}`, or point REDU_METADATA_PATH at an existing parquet file."
    )


def load_redu(max_age_days=30, columns_to_load=None):
    """
    Load ReDU metadata file, checking if it exists and how old it is.
    Downloads the file if it doesn't exist or is older than max_age_days.

    Args:
        max_age_days (int): Maximum age in days before re-downloading the file.
        columns_to_load (list, optional): List of columns to load. If None, loads only essential columns.

    Returns:
        tuple: (pd.DataFrame, datetime) - Processed ReDU metadata DataFrame and file modification datetime,
               or (None, None) if loading fails.
    """

    # Define essential columns used by the application
    if columns_to_load is None:
        columns_to_load = [
            "filename",
            "ATTRIBUTE_DatasetAccession",
            "NCBITaxonomy",
            "UBERONBodyPartName",
            "DOIDCommonName",
            "BiologicalSex",
            "AgeInYears",
        ]

    file_path = get_redu_path()

    if os.path.exists(file_path):
        # Get modification time (ctime on Linux is mtime)
        mtime = os.path.getmtime(file_path)
        age_seconds = time.time() - mtime
        age_days = age_seconds / (24 * 3600)

        if age_days > max_age_days:
            print(f"File is {age_days:.1f} days old, re-downloading...")
            # Keep serving the stale copy if the refresh fails.
            if download_redu_metadata(file_path):
                mtime = os.path.getmtime(file_path)
            else:
                print("Re-download failed, continuing with the existing file.")
    else:
        raw_file = find_raw_redu_file(file_path)
        if raw_file:
            print(f"Found manually downloaded file {raw_file}, converting...")
            if not convert_tsv_to_parquet(raw_file, file_path):
                print(f"Conversion failed. {manual_download_hint(file_path)}")
                return None, None
        else:
            print(f"File does not exist, downloading to {file_path}...")
            if not download_redu_metadata(file_path):
                print(f"Download failed. {manual_download_hint(file_path)}")
                return None, None
        mtime = os.path.getmtime(file_path)

    file_date = datetime.fromtimestamp(mtime)

    try:
        # Load only specified columns for memory efficiency
        df_redu = pd.read_parquet(file_path, columns=columns_to_load)

        # Process ReDU table
        df_redu["filename"] = df_redu["filename"].apply(
            lambda x: x[2:] if x.startswith("f.") else x
        )

        df_redu["filepath"] = (
            df_redu["ATTRIBUTE_DatasetAccession"].astype(str)
            + "/"
            + df_redu["filename"].apply(lambda x: Path(x).stem)
        )

        print(
            f"Loaded ReDU metadata with {len(df_redu):,} rows and {len(df_redu.columns)} columns"
        )
        print(f"Memory usage: {df_redu.memory_usage(deep=True).sum() / 1024**2:.2f} MB")

        return df_redu, file_date
    except Exception as e:
        print(f"Error loading ReDU metadata: {str(e)}")
        return None, None
