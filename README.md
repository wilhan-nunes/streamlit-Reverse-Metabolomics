# Reverse Metabolomics

A Streamlit app for exploring public mass spectrometry data through
[fastMASST](https://fasst.gnps2.org) spectral matching, combined with
[ReDU](https://redu.ucsd.edu) sample metadata (body part, disease, sex, age,
organism) to reveal where a set of molecules occurs across public datasets.

Hosted, ready-to-use version: **[reverse-metabolomics.gnps2.org](https://reverse-metabolomics.gnps2.org)**

Method described in:
> Charron-Lamoureux, V., Mannochio-Russo, H., Lamichhane, S. et al. A guide to
> reverse metabolomics—a framework for big data discovery strategy. *Nat
> Protoc* (2025). https://doi.org/10.1038/s41596-024-01136-2

## What it does

1. Take a list of USIs (Universal Spectrum Identifiers) — one per molecule of
   interest — either typed in or loaded from example data.
2. Run a fastMASST query against public repositories to find matching spectra.
3. Cross-reference matches with ReDU metadata to see which sample types
   (body part, disease, sex, age, organism) they came from.
4. Generate interactive heatmaps (raw, log-transformed, or ReDU-normalized
   counts), with clustering options.
5. Export result tables and figures for downstream use or publication.

## Running it

Most people should just use the hosted site above — no setup required.

To run it on your own machine instead — for development, or so your input
USI lists, merged tables, and figures stay on your computer rather than a
shared server — see **[RUNNING_LOCALLY.md](RUNNING_LOCALLY.md)**. Note that
this does not make the tool usable offline: it works entirely from USIs
(references, not uploaded data), and resolving them via fastMASST plus
looking up ReDU metadata both require an internet connection regardless of
where the app runs.

## Repository layout

| Path | Purpose |
|---|---|
| `app.py` | Main Streamlit application |
| `welcome.py` | Landing page content |
| `masst_sidebar.py`, `masst_client.py` | fastMASST query UI and client |
| `download_redu.py` | Downloads/caches and loads the ReDU metadata file |
| `utils/` | Query-parameter API and dataset metadata helpers |
| `run_local.py` | Entry point for running the app locally via `uv` |
| `Dockerfile`, `docker-compose*.yml`, `Makefile` | Server deployment |
| `tests/` | Automated tests, including external-service health checks |
| `README_API_Params.md` | URL query-parameter configuration reference |
| `RUNNING_LOCALLY.md` | Local setup instructions |

## Configuring via URL

The app's state (USIs, search parameters, analysis settings) can be
configured entirely through URL query parameters, for sharing a
pre-configured link. See
**[README_API_Params.md](README_API_Params.md)** for the full reference.

## Development

```bash
pip install -r requirements.txt
streamlit run app.py
```

Tests:

```bash
pytest tests/
```

Server deployment (Docker) is managed via the `Makefile` targets, e.g.
`make server-compose-production`.

## Citation

If you use this tool, please cite:

> Charron-Lamoureux, V., Mannochio-Russo, H., Lamichhane, S. et al. A guide to
> reverse metabolomics—a framework for big data discovery strategy. *Nat
> Protoc* (2025). https://doi.org/10.1038/s41596-024-01136-2

This application is part of the GNPS downstream analysis ecosystem known as
**MetaboApps**. See the
[full tool index](https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/toolindex/#gnps2-web-tools)
for related tools.
