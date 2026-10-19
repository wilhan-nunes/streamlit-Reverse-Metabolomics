#!/usr/bin/env python3
"""Build the single-file browser edition with a ReDU copy built in.

ReDU sends no CORS headers, so the page cannot download its metadata. This script
writes the 7 ReDU columns the analysis uses into the page's empty
<script id="redu-embedded"> block (gzipped, base64, about 7 MB), along with the date
of that ReDU data. The result is one HTML file that works with no ReDU file at all;
users can still load a newer dump by hand.

    python embed_redu.py                          # download redu.gnps2.org/dump (~0.9 GB)
    python embed_redu.py --input redu.tsv         # a dump, a .tsv.gz, a compact export,
                                                  # or the Streamlit app's REDU_metadata.parquet
    python embed_redu.py --output somewhere.html  # default: build/reverse_metabolomics.html

Running it on an already-built page replaces the data in it, so it also refreshes.
Only the standard library is needed, except pandas + pyarrow for a .parquet input.
"""
import argparse
import base64
import csv
import email.utils
import gzip
import html as htmllib
import io
import os
import re
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
DUMP_URL = "https://redu.gnps2.org/dump"

# Must match REDU_COLS in reverse_metabolomics.html (the columns download_redu.load_redu keeps).
COLUMNS = [
    "filename",
    "ATTRIBUTE_DatasetAccession",
    "NCBITaxonomy",
    "UBERONBodyPartName",
    "DOIDCommonName",
    "BiologicalSex",
    "AgeInYears",
]

BLOCK_RE = re.compile(r'<script type="application/gzip" id="redu-embedded"[^>]*>.*?</script>', re.S)

# Some dump fields are long free text; the default 128 KB limit is too small.
csv.field_size_limit(2**31 - 1)


def open_dump_url(url):
    """Return (text stream, data date) for the live dump. The date is its Last-Modified."""
    resp = urllib.request.urlopen(url)
    modified = resp.headers.get("Last-Modified")
    date = email.utils.parsedate_to_datetime(modified) if modified else datetime.now(timezone.utc)
    return io.TextIOWrapper(resp, encoding="utf-8", newline=""), date


def tsv_rows(stream):
    """Yield the wanted columns of each row. Values are passed through untouched:
    the page applies pandas' missing-value rules itself when it parses them."""
    reader = csv.reader(stream, delimiter="\t")
    header = [h.strip() for h in next(reader)]
    missing = [c for c in COLUMNS if c not in header]
    if missing:
        sys.exit(f"Input is missing the ReDU column(s) {', '.join(missing)}.")
    idx = [header.index(c) for c in COLUMNS]
    last = max(idx)
    for fields in reader:
        if len(fields) <= last:
            continue  # blank or malformed line; the page skips these too
        yield [fields[i] for i in idx]


def parquet_rows(path):
    import pandas as pd

    df = pd.read_parquet(path, columns=COLUMNS)
    for row in df.itertuples(index=False):
        yield ["" if pd.isna(v) else str(v) for v in row]


def read_input(path):
    """Return (row iterator, data date, description) for --input or the live dump."""
    if path is None:
        print(f"Downloading {DUMP_URL} (about 0.9 GB)…", file=sys.stderr)
        stream, date = open_dump_url(DUMP_URL)
        return tsv_rows(stream), date, "redu.gnps2.org/dump"

    path = Path(path)
    date = datetime.fromtimestamp(path.stat().st_mtime, timezone.utc)
    if path.suffix == ".parquet":
        return parquet_rows(path), date, path.name
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    stream = gzip.open(path, "rt", encoding="utf-8", newline="") if gz else \
        open(path, encoding="utf-8", newline="")
    return tsv_rows(stream), date, path.name


def compact_gzip(rows):
    """Write the rows as a TSV in the same format as the page's "Export compact copy"."""
    buf = io.BytesIO()
    # mtime=0 keeps the output identical for identical data.
    with gzip.GzipFile(fileobj=buf, mode="wb", compresslevel=9, mtime=0) as gz:
        text = io.TextIOWrapper(gz, encoding="utf-8", newline="")
        writer = csv.writer(text, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        n = 0
        for n, row in enumerate(rows, 1):
            writer.writerow(row)
            if n % 200000 == 0:
                print(f"  {n:,} rows", file=sys.stderr)
        text.flush()
        text.detach()
    return buf.getvalue(), n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", help="ReDU dump (.tsv or .tsv.gz), compact export, or .parquet; default: download the live dump")
    ap.add_argument("--date", help="ReDU data date (YYYY-MM-DD) to record; default: the dump's Last-Modified, or the input file's modified time")
    ap.add_argument("--template", default=HERE / "reverse_metabolomics.html", type=Path)
    ap.add_argument("--output", default=HERE / "build" / "reverse_metabolomics.html", type=Path)
    args = ap.parse_args()

    page = args.template.read_text(encoding="utf-8")
    if len(BLOCK_RE.findall(page)) != 1:
        sys.exit(f"{args.template} needs exactly one <script id=\"redu-embedded\"> block.")

    t0 = time.time()
    rows, date, source = read_input(args.input)
    if args.date:
        date = datetime.strptime(args.date, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    data, n = compact_gzip(rows)
    if not n:
        sys.exit("No ReDU rows found in the input.")

    b64 = base64.encodebytes(data).decode("ascii")  # 76-character lines keep editors and diffs usable
    block = (
        '<script type="application/gzip" id="redu-embedded"'
        f' data-redu-date="{date.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}"'
        f' data-redu-rows="{n}" data-redu-source="{htmllib.escape(source)}">\n{b64}</script>'
    )
    out = BLOCK_RE.sub(lambda _m: block, page)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.output.with_suffix(".tmp")
    tmp.write_text(out, encoding="utf-8")
    os.replace(tmp, args.output)
    print(
        f"Wrote {args.output} ({len(out.encode()) / 1048576:.1f} MB): {n:,} ReDU rows dated "
        f"{date.date()}, from {source}, in {time.time() - t0:.0f} s.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
