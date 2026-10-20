# eco-audit-ai
https://eco-audit-ai-...streamlit.app

EcoAudit AI: Project Description
Ananya Samadder · October 2026
Project description
EcoAudit AI: an open-source web application for biodiversity screening and spatial analysis of field survey data
EcoAudit AI is a web-based tool, developed in Python and Streamlit, that streamlines the first stage of analysing biodiversity survey data. It accepts species-by-site survey records with geographic coordinates and produces a standardised set of community-ecology metrics, effort-aware site comparisons, an interactive map of sampling locations, and a structured summary report.
Analytical capabilities:
For each site and for the landscape as a whole, the application computes species richness, abundance, Shannon-Wiener diversity (H′), Simpson diversity (plug-in and unbiased estimators), Pielou's evenness (J′), Margalef richness, Hill numbers (q = 0, 1, 2) and the Chao1 richness estimator. Because observed diversity is strongly dependent on sampling effort, the tool also provides individual-based rarefaction curves and richness standardised to a common sample size, allowing sites with unequal sampling effort to be compared on a defensible basis. Analyses are run per taxon group, and the user is warned when taxa with non-comparable abundance units are pooled.
Data integrity and transparency:
Input data are validated and cleaned automatically, including normalisation of species names, handling of missing or invalid values, and coordinate checks. Every modification is recorded in a cleaning log that is included in the output. All formulas are documented with primary references, and the implementations are verified through automated unit tests against analytical results and the scikit-bio library.
Spatial and screening functions:
Sampling sites are displayed on an interactive Folium map with multiple base layers and colour-coded by screening status. Users can define their own diversity thresholds, including habitat-specific values, to identify sites that may warrant follow-up survey effort. The application flags thresholds that are mathematically unreachable given a site's richness (since H′ ≤ ln S), so that screening outputs are not mistaken for assessments of habitat condition.
Outputs
The application generates a downloadable report documenting the settings used, taxon subset analysed, data cleaning log, site-level metrics, methodology and limitations, together with a CSV export of the metrics table.
Scope and limitations:
EcoAudit AI is intended as an educational and exploratory aid for students, researchers and field teams. It is not a regulatory instrument and does not replace statutory Environmental Impact Assessments, protected-species surveys or expert ecological judgement. It does not correct for imperfect detection, seasonal variation or taxonomic uncertainty.
Technical stack and availability:
Python, Streamlit, pandas, NumPy, SciPy, Folium and Plotly, with continuous integration via GitHub Actions. The source code is released under the MIT License at github.com/ananyastack/eco-audit-ai.
