# Running Reverse Metabolomics on your own computer

The tool is available online at
[reverse-metabolomics.gnps2.org](https://reverse-metabolomics.gnps2.org) — that
site needs no setup and is the right choice for most people.

Run it locally instead if you want your analysis to stay on your machine, you
need to work without a stable connection to the hosted site, or you want to
modify the code.

## Setup

Works on Windows, macOS, and Linux. You do **not** need Python installed first —
`uv` provides it.

**1. Install `uv`** (one time). It handles Python and every dependency for you,
inside its own cache, without changing anything else on your system.

- macOS / Linux: `curl -LsSf https://astral.sh/uv/install.sh | sh`
- Windows (PowerShell): `powershell -c "irm https://astral.sh/uv/install.ps1 | iex"`

On Windows, **close PowerShell and open a new window** afterwards, so it picks up
the updated PATH. Check it worked with `uv --version`.

**2. Get the code.** With git:

```bash
git clone https://github.com/wilhan-nunes/streamlit-Reverse-Metabolomics.git
cd streamlit-Reverse-Metabolomics
```

Without git (common on Windows): download the repository ZIP from GitHub via
**Code → Download ZIP**, unzip it, and open a terminal in the unzipped folder.

**3. Start the app:**

```bash
uv run run_local.py
```

That is the whole setup. The app picks an unused port, starts, and opens your
browser at it. To stop, press `Ctrl+C` in the terminal window.

To start it again later, run `uv run run_local.py` again.

Use `uv run run_local.py`, not `python run_local.py`, unless you have already
installed the requirements yourself — on a machine with no Python at all, `uv
run` is what supplies the interpreter. If you would rather manage the
environment yourself, `pip install -r requirements.txt` then `python
run_local.py` works too.

## First run is slow

Two things happen only the first time:

1. `uv` downloads Python and the dependencies (a few minutes).
2. The app downloads the ReDU metadata file from
   [redu.gnps2.org/dump](https://redu.gnps2.org/dump). This is a large file and
   can take several minutes.

The ReDU file is saved to `~/.reverse-metabolomics/REDU_metadata.parquet` — in
your home folder, not in the repository, so it survives re-cloning or moving
the code. It is re-downloaded automatically when it is more than 30 days old.
If a refresh fails, the app keeps using the copy it already has.

### If the automatic download fails

Download <https://redu.gnps2.org/dump> yourself, then point the app at it:

```bash
REDU_METADATA_PATH=/path/to/your/file.parquet uv run run_local.py
```

Note that the dump is a TSV file; the app stores it as Parquet. To convert it:

```python
import pandas as pd
pd.read_csv("dump.tsv", sep="\t", low_memory=False).to_parquet("REDU_metadata.parquet", index=False)
```

## What differs from the hosted site

| | Hosted | Local |
|---|---|---|
| Address | `0.0.0.0:5000` behind nginx | `127.0.0.1`, a free port |
| Reachable from other machines | Yes | **No** — localhost only |
| Usage analytics | Enabled | Disabled |
| ReDU metadata | Shared, cached server-side | Your own copy in `~/.reverse-metabolomics` |

The analysis itself is identical — it is the same `app.py`.

## Settings

| Environment variable | Effect |
|---|---|
| `REVMET_DATA_DIR` | Where the ReDU file is stored (default `~/.reverse-metabolomics`) |
| `REDU_METADATA_PATH` | Use this exact ReDU file; overrides `REVMET_DATA_DIR` |
| `REVMET_LOCAL=1` | Disables usage analytics; set automatically by `run_local.py` |

## Note on what stays online

Running locally does not make the tool offline. MASST searches are executed by
`api.fasst.gnps2.org` and dataset details come from MassIVE, MetaboLights, and
Metabolomics Workbench — searching public repositories is what the tool does,
so those queries necessarily leave your machine. What running locally changes
is that your input USI lists, merged tables, and figures are processed and
stored on your own computer.
