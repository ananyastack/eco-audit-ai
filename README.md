# Veridalis

**A free, open-source tool for automated Darwin Core harmonization, quality control, and biodiversity screening of field survey data.**

[Live app](https://eco-audit-ai-avci9d7appmnzeehkhfppwf.streamlit.app) · [MIT License](LICENSE)

Veridalis standardizes heterogeneous ecological datasets (including raw GBIF occurrence exports and custom field sheets) against Darwin Core (DwC) standards, executes automated spatial quality-control audits, and computes standard community-ecology diversity metrics paired with interactive GIS mapping. It is designed for students, researchers, and field teams who need rapid, transparent screening before downstream biostatistical modeling.

> **Exploratory and educational tool.** Veridalis is not a statutory regulatory instrument and does not replace statutory Environmental Impact Assessments (EIAs), protected-species legal surveys, or expert ecological judgement. See [Limitations](#limitations).

## Features

- **Automated Darwin Core Harmonization:** Ingests standard survey CSVs as well as raw, multi-column tab-separated GBIF occurrence archives. Dynamically resolves common Darwin Core fields (`scientificName`, `decimalLatitude`, `decimalLongitude`, `individualCount`, `locality`).
- **Automated QC & Audit Pipeline:** Detects and flags missing coordinates, inverted latitude/longitude, zero coordinates `(0, 0)`, non-numeric abundance counts, and duplicate headers with a live usability score.
- **Diversity Indices (Site & Landscape):** Computes species richness ($S$), individual abundance ($N$), Shannon-Wiener diversity ($H'$ using natural logarithms), Simpson diversity ($1 - D$), Pielou evenness ($J'$), and Margalef richness index.
- **Interactive Spatial Mapping:** Folium-powered GIS map featuring marker clustering, CartoDB Light/Dark, OpenStreetMap, and Esri Satellite basemaps. Markers are classified and color-coded by diversity status with site-level telemetry on click.
- **Screening & Audit Reporting:** Configurable Shannon-Wiener and richness screening thresholds that flag simplified assemblages. Generates downloadable audit reports (`.txt`) and direct draft email dispatch via Gmail.
- **Built-in Demonstration Data:** Ships with synthetic riparian and woodland benchmark datasets alongside immediate GBIF file drop compatibility.

## Input formats

Veridalis accepts both comma-separated (`.csv`) and tab-separated (`.tsv`, `.txt`) text files.

### 1. Standard Survey CSV
| Column | Type | Required | Description |
|---|---|---|---|
| `site_id` | String | Yes | Unique site or sampling station identifier (aliases: `station`, `site`) |
| `latitude` | Float | Yes | WGS84 decimal latitude ($-90$ to $+90$; aliases: `lat`, `decimalLatitude`) |
| `longitude` | Float | Yes | WGS84 decimal longitude ($-180$ to $+180$; aliases: `lon`, `lng`, `decimalLongitude`) |
| `species` | String | Yes | Binomial taxon name (aliases: `scientific_name`, `taxon_name`, `scientificName`) |
| `abundance` | Numeric | Yes | Count or individual tally (aliases: `count`, `n`, `individualCount`) |
| `site_name` | String | Optional | Descriptive site name |
| `taxon` | String | Optional | Higher taxonomic grouping (e.g., Plant, Insect, Fish) |
| `habitat` | String | Optional | Habitat type or vegetation class |
| `date` | Date/String | Optional | Sampling date (`YYYY-MM-DD`) |

### 2. GBIF Occurrence Exports
Raw occurrence dumps downloaded directly from GBIF can be uploaded as-is. The embedded `harmonize.py` engine strips duplicate headers, extracts coordinates, imputes default abundances for occurrence-only records, and reports flagged anomalous records without crashing.

## Methods

| Metric | Formula | Reference |
|---|---|---|
| Shannon-Wiener $H'$ | $-\sum p_i \ln p_i$ | Shannon (1948) |
| Simpson $1 - D$ | $1 - \sum p_i^2$ | Simpson (1949) |
| Pielou evenness $J'$ | $H' / \ln S$ | Pielou (1966) |
| Margalef richness | $(S - 1) / \ln N$ | Margalef (1958) |

## Quality control & Screening flags

- **Quality Flags:** Records with missing coordinates, non-numeric values, or coordinates outside $[-90, 90]$ and $[-180, 180]$ are quarantined by the harmonization engine to preserve audit trail integrity while maintaining app stability.
- **Screening Flags:** A screening flag indicates that a site fell below user-defined $H'$ or $S$ thresholds. It highlights simplified communities or potential sampling deficits; it is not an indicator of legal non-compliance or environmental degradation.

## Limitations

- **Sampling Effort:** Diversity metrics are unstandardized for sample size or unequal sampling effort across sites.
- **Pooled Taxa:** Abundances across differing taxonomic groups (e.g., plant stem counts vs. avian call tallies) are pooled.
- **Single-Taxon Datasets:** Single-species occurrence archives (e.g., targeted taxon exports) naturally produce $H' = 0$ and $J' = 0$, triggering low-diversity screening flags as designed.
- **Spatial Autocorrelation:** No spatial distance weights or imperfect detection corrections are applied.

## Roadmap

- Hill numbers ($q = 0, 1, 2$) and Chao1 richness estimator
- Individual-based sample rarefaction curves
- Per-taxon isolated stratification
- Automated GBIF Taxonomic Backbone species reconciliation via API

## Run locally

```bash
git clone [https://github.com/ananyastack/eco-audit-ai.git](https://github.com/ananyastack/eco-audit-ai.git)
cd eco-audit-ai
pip install -r requirements.txt
streamlit run app.py
