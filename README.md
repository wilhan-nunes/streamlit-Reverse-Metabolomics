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

## Browser edition

[`browser_app/reverse_metabolomics.html`](browser_app/reverse_metabolomics.html)
is a single-file version of the app that does all of its analysis in the web
browser. It has no build step and needs no Python, and it is meant for the
GNPS2 web apps site (apps.gnps2.org). It shows every figure and table the
Streamlit app does, with a few additions:

- **No limit on the number of USIs** (the Streamlit app allows 20).
- **ReDU metadata stays on your machine.** You drop the dump from
  [redu.gnps2.org/dump](https://redu.gnps2.org/dump) into the page once, and
  a compact copy is kept in your browser for later visits.
- **Biological sex and age** are available as analysis variables, alongside
  body part and disease.
- **Collapsible search configuration**, and USI lists can be imported or
  exported as CSV/TSV.

It uses the same URL parameters as the Streamlit app, so shared links work in
both. FASST only accepts requests from pages served on `*.gnps2.org`, so test
it from the apps.gnps2.org sandbox rather than opening the file locally.

### Versioning and changelog

The browser edition has its own version, a release date (`YYYY-MM-DD`) like
the Streamlit app's. The version is set in one place, the
`<meta name="version">` tag in the HTML file, and the page footer reads it
from there. Every change is recorded in
**[browser_app/CHANGELOG.md](browser_app/CHANGELOG.md)**. To release a change,
move the changelog's *Unreleased* entries under a new dated heading and set
the same date in the meta tag.

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
| `browser_app/reverse_metabolomics.html` | Browser edition: the whole app in one HTML file |
| `browser_app/CHANGELOG.md` | Version history of the browser edition |
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

> Mannochio-Russo, H., Gonçalves Nunes, W.D., Zhao, H.N. et al. Bridging
> complexity and accessibility in metabolomics with MetaboApps. *Nat
> Methods* 23, 1475–1477 (2026). https://doi.org/10.1038/s41592-026-03186-8

This application is part of the GNPS downstream analysis ecosystem known as
**MetaboApps**. See the
[full tool index](https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/toolindex/#gnps2-web-tools)
for related tools.
