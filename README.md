# Veridalis

Open-source web app for harmonizing biodiversity occurrence data, checking its quality, and screening it with standard diversity indices.

[Live app](https://eco-audit-ai-avci9d7appmnzeehkhfppwf.streamlit.app) · [Video walkthrough](https://youtu.be/pPBQ-fytl7c) · [MIT License](LICENSE) · [Citation](#citation)

[![Veridalis video walkthrough](https://img.youtube.com/vi/pPBQ-fytl7c/maxresdefault.jpg)](https://youtu.be/pPBQ-fytl7c)

Veridalis reads field survey tables and GBIF occurrence downloads, maps their columns to standard names (Darwin Core terms where they exist), runs quality checks, and reports richness, Shannon, Simpson, Pielou and Margalef indices per site and for the whole dataset, with an interactive map and downloadable tables. It is meant for students, researchers and field teams who want a quick, transparent look at a file before any further analysis.

> **Version 0.1.0. Exploratory research software; not a validated tool.** It does not replace expert ecological judgement, statutory environmental assessments or protected-species surveys. Read [Limitations](#limitations) and [Known issues](#known-issues) before interpreting any output.

The repository is named `eco-audit-ai`; the software is called Veridalis.

## Features

- **Flexible reading.** Comma-, tab-, semicolon- and pipe-separated text files (`.csv`, `.tsv`, `.txt`), UTF-8 or Latin-1, decimal commas, up to 200,000 rows. Zipped GBIF downloads must be unzipped first.
- **Column recognition.** Recognises common names such as `scientificName`, `decimalLatitude`, `decimalLongitude`, `individualCount` and `locality`, and shows the mapping it used. If a required column is not recognised, you choose it from a list.
- **Quality checks.**
  - *Excluded from analysis:* records with no taxon name, missing or unreadable coordinates, coordinates at exactly (0, 0), coordinates out of range, suspected swapped latitude and longitude, non-positive counts, and exact duplicate records (duplicate exclusion can be switched off in the sidebar).
  - *Warning only, record kept:* missing, unreadable or future dates.
  - For GBIF files, GBIF's own issue flags are listed separately for reference. They are not used to exclude records.
- **Diversity indices** per site and for the whole dataset: species richness (S), individuals (N), Shannon-Wiener H′ (natural log), Simpson 1 − D, Pielou J′ and Margalef richness.
- **Sites.** Site identifiers and names are used exactly as written. Records with neither a site identifier nor a site name are grouped into 0.1° grid cells named `cell_<lat>_<lon>`.
- **Unique species per habitat.** A table of unique species, records and sites for each habitat.
- **Screening.** Sites can be flagged for review with absolute thresholds (H′ and S) or relative to the dataset (lowest-H′ sites). The High / Moderate / Low class bands are adjustable. The rule in use is shown above the table and in the report.
- **Map.** Interactive map with markers coloured by diversity class, clustered markers and switchable basemaps.
- **Outputs.** Cleaned analysis table, flagged records and site metrics (CSV); a screening report with tables (Markdown, previewed in the app); an Excel workbook with one sheet per table (needs `openpyxl`); a Gmail draft button.
- **Single-species warning.** If a file contains one species, the app states that the diversity indices are not meaningful.

## Input

Header matching ignores case and punctuation. Only species, latitude and longitude are required.

| Field | Required | Recognised names |
|---|---|---|
| `species` | Yes | `species`, `scientificName`, `scientific_name`, `taxon_name`, `verbatimScientificName` |
| `latitude` | Yes | `latitude`, `decimalLatitude`, `lat`, `decimal_latitude`, `latitude_dd` |
| `longitude` | Yes | `longitude`, `decimalLongitude`, `lon`, `lng`, `long`, `decimal_longitude`, `longitude_dd` |
| `abundance` | No | `abundance`, `individualCount`, `count`, `n`, `individuals`, `organismQuantity` |
| `site_id` | No | `site_id`, `station`, `station_id`, `locationID`, `site` |
| `site_name` | No | `site_name`, `locality`, `verbatimLocality`, `location` |
| `date` | No | `date`, `eventDate`, `event_date`, `observed_on` |
| `habitat` | No | `habitat`, `verbatimHabitat` |
| `taxon` | No | `taxon`, `class` (a higher group such as Fish or Plant) |

A CSV template is available in the app sidebar.

## How counts and names are handled

- **Missing counts.** Many occurrence records have no count. A missing or unreadable count is treated as presence (n = 1) and is not flagged. The audit panel and the report state how many records were treated this way. Indices computed from such records reflect how often species were recorded, not how abundant they are. A count of zero or less excludes the record.
- **Passed.** "Passed" in the audit panel means a record has a taxon name, valid coordinates and a positive count. It does not mean the record is ecologically correct.
- **Names.** Only completely empty cells count as missing in site and habitat fields. Values such as `NA`, `None` or `unknown` are kept as written. Species names written as two words are standardised to *Genus epithet* capitalisation, and extra spaces are removed. If one site identifier has several names, they are shown joined with " / ".
- **Duplicates.** A duplicate is a record with the same species, location (5 decimal places), date and count as an earlier record.

## Methods

| Metric | Formula | Reference |
|---|---|---|
| Shannon-Wiener H′ | −Σ pᵢ ln pᵢ | Shannon (1948) |
| Simpson 1 − D | 1 − Σ pᵢ² | Simpson (1949) |
| Pielou evenness J′ | H′ / ln S | Pielou (1966) |
| Margalef richness | (S − 1) / ln N | Margalef (1958) |

pᵢ is the proportion of individuals belonging to species i. N is the sum of the counts, so it can exceed the number of records. H′ cannot exceed ln S, so species-poor sites score low whatever their habitat quality. J′ is undefined when S = 1.

**Screening flags.**
- *Absolute mode (default):* a site is marked "Review suggested" when H′ is below a threshold or S is below a minimum (defaults: H′ < 2.00 or S < 4).
- *Relative mode:* a site is marked when its H′ is at or below a chosen percentile of the site values in the dataset (default: 25th). If every site has the same H′, no sites are flagged.
- *Class bands:* defaults are High ≥ 3.00, Moderate ≥ 2.00, Low ≥ 1.00, otherwise Very low. They only set the class labels and map colours.

All thresholds and bands are rough screening guides, not ecological standards, and they vary by taxon and habitat. A flag means only that a site fell outside the rule you chose. It is not a legal determination and says nothing about environmental condition.

## Tested so far

| Dataset | Records | Species | Sites | Result |
|---|---|---|---|---|
| GBIF occurrence export, Slovenian freshwater records | 92 read, 90 analysed, 2 excluded (1 missing coordinates, 1 duplicate) | 1 (*Cyprinus carpio*) | 82 | 35 records without a date (warning only); 87 of the 90 analysed records had no count; single-species warning shown |
| GBIF occurrence export, Polish salmonids | 26 read, 26 analysed | 4 | 12 grid cells | Whole dataset: H′ = 0.939, 1 − D = 0.512, J′ = 0.677, Margalef = 0.921; no records flagged |

For the Polish file, the whole-dataset values match a hand calculation from the species counts (17, 6, 2 and 1 records). The Slovenian file is a single-species set, so it tests the checks, not the diversity indices.

The app was also run on small synthetic files with semicolon separators, decimal commas, unusual headers, blank counts, impossible coordinates, duplicates, and site names such as `NA`, `None` and `unknown`.

GBIF download DOIs for the two test files: add before release.

There is no formal test suite or benchmark yet.

## Limitations

- **Occurrence data are not surveys.** GBIF data are opportunistic records. Richness and evenness calculated from them are not survey estimates.
- **Sampling effort.** Indices are not standardised for sample size or unequal effort between sites.
- **Pooled taxa.** Abundances of different groups are pooled and are not comparable.
- **Sparse grid cells.** A 0.1° cell with one species always has S = 1, so in a sparse file most per-cell flags reflect the grid, not the habitat. Read them together with the whole-dataset values.
- **No uncertainty estimates.** No confidence intervals, spatial autocorrelation or detection corrections are applied.
- **No taxonomic matching.** Names are not reconciled with a taxonomic backbone, so synonyms and misspellings count as separate species.
- **Other record problems are not checked,** for example coordinate uncertainty, rounded coordinates, or records outside a species' known range. Duplicate detection is exact only.

## Known issues

- For sites with S = 1, J′ is shown as 0 instead of "n/a".
- In absolute mode, a file with one species still shows site-level flags.
- Long site names make tables wide.
- Basemaps come from third-party tile servers and need an internet connection.

## Run locally

You do not need to install anything to try Veridalis: use the [live app](https://eco-audit-ai-avci9d7appmnzeehkhfppwf.streamlit.app). Run it on your own computer if your data are sensitive, if you want to work offline (map tiles still need internet), or if you want to change the code.

You need [Python](https://www.python.org/downloads/) 3.11 or newer. On Windows, tick "Add python.exe to PATH" when installing. Veridalis was developed and tested on Windows with Python 3.14.

### Option 1: launcher (first run installs everything)

1. On GitHub choose **Code → Download ZIP** and extract it.
2. **Windows:** double-click `run_veridalis.bat`. **macOS or Linux:** run `bash run_veridalis.sh` in a terminal.

The launcher creates a private environment, installs the packages from `requirements.txt`, and opens the app at http://localhost:8501.

### Option 2: manual

```bash
git clone https://github.com/ananyastack/eco-audit-ai.git
cd eco-audit-ai
```

**Windows (PowerShell)**

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m streamlit run app.py
```

**macOS and Linux**

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m streamlit run app.py
```

These commands call the environment's own Python, so nothing depends on your PATH. Press Ctrl+C in the terminal to stop the app. On later runs only the last command is needed.

### Troubleshooting

| Message | Fix |
|---|---|
| `pip` or `streamlit` "is not recognized" | Use the `python -m ...` commands above instead of calling `pip` or `streamlit` directly |
| `Python was not found; run without arguments to install from the Microsoft Store` | Install Python from python.org, or turn off the Python entries under Settings > Apps > Advanced app settings > App execution aliases |
| `py` is not recognized | Install Python from python.org with the "py launcher" option ticked, or use the full path to `python.exe` |
| Long build error while installing pandas, numpy or pyarrow | Your Python version may be too new for the packages. Install Python 3.12 and create the environment with `py -3.12 -m venv .venv` |
| `Port 8501 is already in use` | Add `--server.port 8502` to the last command |
| The app says Excel export needs `openpyxl` | Run `.\.venv\Scripts\python.exe -m pip install openpyxl` (macOS and Linux: `.venv/bin/python -m pip install openpyxl`) |

If something else goes wrong, please [open an issue](https://github.com/ananyastack/eco-audit-ai/issues) with your operating system, your Python version and the last 15 lines of the error.

## Privacy

Veridalis does not save uploaded files; they stay in memory for your session. The hosted app runs on Streamlit Community Cloud. For sensitive data, such as locations of threatened species, run Veridalis on your own computer.

## Planned work

- Unit tests and a benchmark against an established reference implementation
- Show J′ as n/a and suppress site-level flags for single-species files
- Hill numbers (q = 0, 1, 2) and the Chao1 richness estimator
- Individual-based rarefaction curves
- Per-taxon analysis
- Taxonomic name reconciliation through the GBIF API

## AI assistance

AI coding assistants were used during development. The author tested the app on the datasets listed above and checked its outputs, including a hand calculation for the Polish file.

## Citation

If you use this software, please cite it using [`CITATION.cff`](CITATION.cff):

> Samadder, A. (2026). *Veridalis: biodiversity data harmonization, quality control and screening* (version 0.1.0) [Software]. https://github.com/ananyastack/eco-audit-ai

## License

MIT. See [LICENSE](LICENSE).                       
