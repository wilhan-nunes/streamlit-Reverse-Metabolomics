# Changelog — Reverse Metabolomics, browser edition

All notable changes to `reverse_metabolomics.html` are recorded here.

Versions are release dates (`YYYY-MM-DD`), the same scheme the Streamlit app
uses. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

**Releasing a change:**

1. Move the entries under *Unreleased* into a new dated section.
2. Set the same date in the app's `<meta name="version" content="…">` tag. The
   footer reads its version from that tag, so there is nothing else to edit.

## [Unreleased]

## [2026-10-14]

First release: a single-file, in-browser port of the Streamlit app (version
2026-08-14).

### Added

- Every figure and table from the Streamlit app:
  - counts table
  - heatmap: raw, log2 or ReDU-normalized; ten clustering metrics; row and
    column clustering with dendrograms; custom size
  - single-compound bar chart
  - dataset distribution chart
  - merged-data table with mirror-plot links
  - dataset details from MetaboLights and Metabolomics Workbench
- CSV downloads for every table and SVG downloads for every chart.
- No limit on the number of USIs (the Streamlit app allowed 20). Searches run
  three at a time, with per-USI progress and a Cancel button.
- Collapsible search configuration that folds into a one-line summary after a
  search.
- USI list import (paste, or a CSV/TSV file) and export.
- ReDU metadata is read in the browser from the dump at
  `redu.gnps2.org/dump`. A compact copy is saved in the browser (IndexedDB),
  so later visits don't need the file again.
- "Export compact copy": a gzipped 7-column ReDU file, a few MB instead of
  ~0.9 GB, that loads like the full dump.
- Biological sex and age as analysis variables. The Streamlit app requested
  columns it never loaded, so it offered only body part and disease.
- Shareable links that use the same URL parameters as the Streamlit app.
- Version number and changelog link in the footer.

### Changed (compared with the Streamlit app)

- Clustering no longer fails when a distance is undefined (for example,
  correlation with a constant row). Those distances are set to the largest
  finite distance, and the page says so.
- The dataset distribution chart shows the top 50 datasets by default
  (25 / 50 / 100 / all). The CSV download always includes every dataset.
- The merged table is paginated (25 rows per page) and sortable by any column,
  rather than showing only the first 100 rows.

### Known limitations

- FASST only accepts requests from pages served on `*.gnps2.org`, so searches
  fail when the file is opened locally or from localhost. The page explains
  this when it happens.
- MassIVE does not allow cross-site requests, so MassIVE datasets show links
  instead of title, summary and species. Set `MASSIVE_ALLOWS_CORS = true` if
  that changes.
- ReDU-normalized body-part values divide by ReDU's own category counts, before
  the skin and blood sub-sites are merged. This is kept deliberately so results
  match the Streamlit app.

[Unreleased]: https://github.com/wilhan-nunes/streamlit-Reverse-Metabolomics/commits/main/browser_app
[2026-10-14]: https://github.com/wilhan-nunes/streamlit-Reverse-Metabolomics/commits/main/browser_app
