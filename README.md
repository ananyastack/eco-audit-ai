# Veridalis

Open-source exploratory tool for harmonizing biodiversity occurrence data to Darwin Core, checking data quality, and screening assemblages with standard community diversity indices.

[Live app](https://eco-audit-ai-avci9d7appmnzeehkhfppwf.streamlit.app) · [Video walkthrough](https://youtu.be/pPBQ-fytl7c) · [MIT License](LICENSE) · [Citation](#citation)[cite: 1]

[![Veridalis video walkthrough](https://img.youtube.com/vi/pPBQ-fytl7c/maxresdefault.jpg)](https://youtu.be/pPBQ-fytl7c)[cite: 1]

Veridalis reads field survey tables and GBIF occurrence downloads, maps their columns to Darwin Core terms, runs quality checks, and reports richness, Shannon, Simpson, Pielou, and Margalef indices per site and for the whole dataset, alongside an interactive spatial preview. It is designed as an initial data hygiene and exploratory screening step for students, researchers, and field teams before undertaking formal statistical modeling.

> **Version 0.1.0 (Exploratory research software).** Veridalis is a first-pass data hygiene and screening tool[cite: 1]. It does not perform statutory ecological impact assessments, validate species taxonomy, or replace comprehensive community distribution modeling[cite: 1]. Review [Scope & Operational Realities](#scope--operational-realities) before interpreting any output.

The repository is named `eco-audit-ai`; the software is called Veridalis.

---

## Features

- **Flexible Ingestion:** Reads comma-, tab-, semicolon-, and pipe-delimited text files (`.csv`, `.tsv`, `.txt`), including GBIF occurrence exports, and handles decimal commas. If columns are not detected automatically, an interactive manual mapping interface is provided. Zipped GBIF downloads must be unzipped prior to upload.
- **Darwin Core Harmonization:** Maps recognized columns to standard Darwin Core terms (`scientificName`, `decimalLatitude`, `decimalLongitude`, `individualCount`, `locality`, etc.) and displays the active column mapping.
- **Data Quality Audits:**
  - *Excluded from analysis:* Records missing taxon names, missing or zero coordinates `(0, 0)`, coordinates outside valid geographic ranges, and exact duplicates (configurable in the sidebar).
  - *Flagged:* Suspected swapped latitude and longitude.
  - *Warning only (retained):* Missing sampling dates.
  - *External audit flags:* For GBIF exports, GBIF's internal issue tags are parsed and summarized separately for reference.
- **Community Diversity Indices:** Species richness ($S$), total abundance ($N$), Shannon–Wiener ($H'$), Gini–Simpson ($1 - D$), Pielou evenness ($J'$), and Margalef richness, calculated per site and dataset-wide.
- **Fallback Spatial Aggregation:** Records lacking a site identifier or locality name are aggregated into $0.1^\circ$ grid cells named `cell_<lat>_<lon>`.
- **Habitat Summary:** Summarizes unique species counts, total occurrence records, and sampling site representation across distinct habitat classes[cite: 1].
- **Interactive Mapping:** Spatial visualization with marker clustering, status-based coloration, and free basemap tile layers (Esri Light Gray, Esri Dark Gray, OpenStreetMap, and Esri World Imagery) that require no API keys[cite: 1].
- **Dual-Mode Screening & Multi-Format Exports:**
  - Screen sites using absolute cut-offs (e.g., $H' < 2.00$ or $S < 4$) or relative percentiles (e.g., lower quartile $Q1$ of the dataset's $H'$ distribution)[cite: 1].
  - Downloadable cleaned analysis tables (CSV) and flagged records (CSV)[cite: 1].
  - Downloadable Markdown screening audit report (`.md`) and direct email draft creation via Gmail[cite: 1].
  - Multi-tab Excel audit workbook (`.xlsx`) generated via `openpyxl`, organizing dataset summaries, site metrics, screening flags, habitat summaries, and cleaned records[cite: 1].
- **Single-Species Guard:** Detects datasets or assemblages with only one recorded species ($S = 1$) and displays an explicit notice that diversity indices cannot evaluate community variation[cite: 1].
- **Ephemeral Session Security:** Uploaded tables are processed purely in temporary session memory and are never written to disk or an external database[cite: 1].
- **CSV Template:** Downloadable standard template available directly from the sidebar[cite: 1].

---

## Input Formats

Veridalis accepts delimited text files (`.csv`, `.tsv`, `.txt`)[cite: 1]. Only three columns are strictly required[cite: 1]:

| Target Field | Required | Description & Accepted Aliases |
|---|---|---|
| `species` | Yes | Taxon binomial or scientific name (`scientificName`, `scientific_name`, `taxon_name`)[cite: 1] |
| `latitude` | Yes | Decimal latitude in WGS84, −90 to 90 (`lat`, `decimalLatitude`, `latitude_dd`)[cite: 1] |
| `longitude` | Yes | Decimal longitude in WGS84, −180 to 180 (`lon`, `lng`, `decimalLongitude`, `longitude_dd`)[cite: 1] |
| `abundance` | No | Count of individuals (`individualCount`, `count`, `n`). Unrecorded or non-numeric counts default to presence ($n = 1$)[cite: 1] |
| `site_id` | No | Unique sampling site or station code (`station`, `site`). If absent, $0.1^\circ$ grid cells are assigned[cite: 1] |
| `site_name` | No | Verbatim locality description (`locality`, `verbatimLocality`)[cite: 1] |
| `taxon` | No | Higher taxonomic rank or group (e.g., Fish, Plant, Insect)[cite: 1] |
| `habitat` | No | Habitat or vegetation type (`verbatimHabitat`)[cite: 1] |
| `date` | No | Observation or collection date (`eventDate`, `observed_on`)[cite: 1] |

---

## Handling Counts & Presences

Many historical and opportunistic occurrence records lack individual specimen counts[cite: 1]. Veridalis treats missing or unreadable counts as presence ($n = 1$) without discarding the record[cite: 1]. The audit interface and generated reports state the exact number of presences evaluated this way[cite: 1].

*Note:* Total individuals ($N$) represents the sum of abundances across records, which may exceed the total number of rows[cite: 1]. A "Passed" status in the audit panel indicates valid spatial and taxonomic coordinates; it does not constitute ecological verification[cite: 1].

---

## Methods & Formulas

| Metric | Formula | Reference |
|---|---|---|
| Shannon–Wiener $H'$ | $-\sum p_i \ln p_i$ | Shannon (1948)[cite: 1] |
| Gini–Simpson $1 - D$ | $1 - \sum p_i^2$ | Simpson (1949)[cite: 1] |
| Pielou Evenness $J'$ | $H' / \ln S$ | Pielou (1966)[cite: 1] |
| Margalef Richness | $(S - 1) / \ln N$ | Margalef (1958)[cite: 1] |

$p_i$ is the proportion of total individuals belonging to species $i$[cite: 1]. Because $H'$ cannot exceed $\ln(S)$, species-poor sites will score low regardless of habitat condition[cite: 1]. Pielou's $J'$ is mathematically undefined when $S = 1$[cite: 1].

**Screening Flags:** Sites are flagged for review based on either[cite: 1]:
1. *Absolute criteria:* $H'$ or $S$ falling below user-defined cut-offs (defaults: $H' < 2.00$ or $S < 4$)[cite: 1].
2. *Relative criteria:* Sites falling at or below a user-selected percentile (e.g., the lower quartile $Q1$) of the dataset's $H'$ distribution[cite: 1].

A flag denotes only that a site fell within the lower screening tail; it does not serve as a regulatory determination of ecological degradation[cite: 1].

---

## Scope & Operational Realities

To ensure appropriate scientific interpretation, users should observe the following domain considerations[cite: 1]:

- **Opportunistic Records vs. Standardized Surveys:** GBIF downloads represent presence-only collection events rather than effort-controlled transects or plots[cite: 1]. Metrics calculated from presence data reflect recording activity alongside community composition[cite: 1].
- **Single-Species Mathematical Limits ($S = 1$):** At sites with only one detected species, $H'$ is $0.0$ and $J' = \frac{0}{\ln(1)}$ is undefined[cite: 1]. Veridalis renders $J'$ as `0.0` to avoid division-by-zero application exceptions[cite: 1].
- **Grid Cell Artifacts:** In sparse datasets without locality names, an isolated record placed into a $0.1^\circ$ cell yields $S = 1$, triggering an absolute screening flag[cite: 1]. Site-level flags should always be cross-referenced against whole-dataset values[cite: 1].
- **No Taxonomic Backbone Scrubbing:** The software capitalizes two-part names and strips whitespace, but does not reconcile synonyms against external authorities (e.g., the GBIF Backbone Taxonomy)[cite: 1]. Synonyms or spelling variants are evaluated as independent species[cite: 1].
- **No Effort Standardization:** The application does not apply rarefaction curves, detectability adjustments, or spatial autocorrelation models[cite: 1].

---

## Tested Datasets

| Dataset | Records | Species | Sites | Verification Result |
|---|---|---|---|---|
| GBIF Slovenian Freshwater Records | 92 read, 90 analyzed, 2 excluded (1 missing coordinates, 1 duplicate)[cite: 1] | 1 (*Cyprinus carpio*)[cite: 1] | 82[cite: 1] | 35 missing dates flagged (non-critical); 87 records defaulted to presence ($n=1$); single-species guard triggered[cite: 1] |
| GBIF Polish Salmonids (TSV export)[cite: 1] | 26 read, 26 analyzed[cite: 1] | 4[cite: 1] | 12 grid cells[cite: 1] | $H' = 0.939$, $1 - D = 0.512$, $J' = 0.677$ dataset-wide; results verified against hand calculation[cite: 1] |

---

## Run Locally (Offline First)

For sensitive spatial records (such as protected or threatened species localities), Veridalis can be executed entirely offline on your local machine[cite: 1]:

### No Terminal Commands Required
1. Download this repository by clicking **Code** $\rightarrow$ **Download ZIP**, and extract the archive[cite: 1].
2. Ensure [Python 3.11+](https://www.python.org/downloads/) is installed (on Windows, ensure **"Add python.exe to PATH"** was checked during installation)[cite: 1].
3. Run the launcher script:
   - **Windows:** Double-click `run_veridalis.bat`
   - **macOS / Linux:** Run `bash run_veridalis.sh` in Terminal

The script will automatically set up the virtual environment (`.venv`), install dependencies from `requirements.txt`, and open the application in your browser at `http://localhost:8501`[cite: 1].

### Manual Setup

**Windows (PowerShell)**
```powershell
git clone [https://github.com/ananyastack/eco-audit-ai.git](https://github.com/ananyastack/eco-audit-ai.git)
cd eco-audit-ai
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m streamlit run app.py
