# EcoAudit AI (Veridalis)

**A free, open-source tool for quick biodiversity screening of field survey data.**

[Live app](https://eco-audit-ai-avci9d7appmnzeehkhfppwf.streamlit.app) · [MIT License](LICENSE)
EcoAudit AI takes a species-by-site survey CSV and returns standard community-ecology diversity metrics, an interactive map of sampling sites, and a downloadable screening report. It is designed for students, early-career researchers and field teams who want a transparent first look at their data before a full analysis.

> **Exploratory and educational tool.** It is not a regulatory instrument and does not replace statutory Environmental Impact Assessments, protected-species surveys, or expert ecological judgement. See [Limitations](#limitations).

## Features

- **Diversity indices per site and for the whole landscape:** species richness (S), abundance (N), Shannon-Wiener H′ (natural log), Simpson 1−D, Pielou evenness J′, and Margalef richness.
- **Interactive map** (Folium) with street, light, dark and satellite layers. Site markers are clustered and colour-coded by diversity class, with metrics shown on click.
- **Adjustable screening thresholds** for H′ and minimum species richness, to flag sites that may warrant a follow-up survey.
- **Screening report (.txt)** summarising the settings, site metrics and flags.
- **Built-in demonstration data** (synthetic) and a downloadable CSV template.

## Input format

Required columns: `site_id`, `latitude`, `longitude`, `species`, `abundance`
Optional columns: `site_name`, `taxon`, `habitat`, `date`, `notes`

Common aliases (`lat`, `lon`, `count`, `scientific_name`, etc.) are recognised.

## Methods

| Metric | Formula | Reference |
|---|---|---|
| Shannon-Wiener H′ | −Σ pᵢ ln pᵢ | Shannon (1948) |
| Simpson 1−D | 1 − Σ pᵢ² | Simpson (1949) |
| Pielou evenness J′ | H′ / ln S | Pielou (1966) |
| Margalef richness | (S − 1) / ln N | Margalef (1958) |

## Interpreting screening flags

A flag means only that a site falls below the thresholds you set, and that a follow-up survey may be worth considering. It is not a judgement of legal compliance or habitat condition. Note that H′ cannot exceed ln S, so a site with few species will always score low on H′ regardless of habitat quality. The demonstration data is synthetic and should not be interpreted ecologically.

## Limitations

- Indices depend on sampling effort, season, detectability and taxonomic resolution. The current version does not standardise for unequal effort.
- Taxa are pooled when computing indices, so abundance units across groups (for example plant counts and bird counts) are not comparable.
- Commonly quoted H′ bands vary by taxon and habitat and should be treated as rough guides.
- No correction for imperfect detection or spatial autocorrelation.
- Outputs must not be the sole basis for regulatory, legal or land-management decisions.

## Roadmap

- Hill numbers and the Chao1 estimator
- Individual-based rarefaction for comparing sites with unequal effort
- Per-taxon analysis
- Data cleaning log and species-name normalisation
- Automated tests validated against scikit-bio

## Run locally

```bash
git clone https://github.com/ananyastack/eco-audit-ai.git
cd eco-audit-ai
pip install -r requirements.txt
streamlit run app.py
```

## Tech stack

Python · Streamlit · pandas · NumPy · Folium / streamlit-folium

## References

- Margalef, R. (1958). Information theory in ecology. *General Systems*, 3, 36–71.
- Pielou, E. C. (1966). The measurement of diversity in different types of environments. *Journal of Theoretical Biology*, 13, 131–144.
- Shannon, C. E. (1948). A mathematical theory of communication. *Bell System Technical Journal*, 27, 379–423.
- Simpson, E. H. (1949). Measurement of diversity. *Nature*, 163, 688.

## License

Released under the [MIT License](LICENSE). © 2026 Ananya Samadder.
