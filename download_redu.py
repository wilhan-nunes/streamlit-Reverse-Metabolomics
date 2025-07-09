import tempfile
import shutil
from pathlib import Path
import requests
import time


def download_redu_metadata(output_path):
    """
    Download ReDU metadata file from GNPS2 using streaming.
    Downloads to a temporary file first, then moves to final location on completion.

    Args:
        output_path: Path where to save the downloaded file
    """
    url = "https://redu.gnps2.org/dump"
    output_path = Path(output_path)

    # Create directory if it doesn't exist
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Create temporary file in the same directory as the output
    temp_file = None
    try:
        # Create temporary file in the same directory as the final output
        with tempfile.NamedTemporaryFile(
                dir=output_path.parent,
                delete=False,
                suffix='.tmp'
        ) as temp_file:
            temp_path = Path(temp_file.name)

            # Make request with stream enabled
            with requests.get(url, stream=True) as response:
                response.raise_for_status()

                # Get total file size if available
                total_size = int(response.headers.get('content-length', 0))

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

                            print(f"\rDownloading: {progress:.1f}% ({speed:.1f} MB/s)", end="")
                    print()  # New line after download completes

        # Move temporary file to final location only after successful download
        shutil.move(temp_path, output_path)
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