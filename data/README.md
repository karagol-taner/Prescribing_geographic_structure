# Data

All files are aggregates of national open data; no patient-level data are present. Areas are the 42 Integrated Care Boards (ICBs), identified by their ODS code (for example QWO, West Yorkshire); `icb_order.csv` gives the ONS codes and names.

## Inputs
- `icb_area_substance_items_full.csv` - columns AREA_CODE, AREA_NAME, DRUG, SECTION, ITEMS: dispensed items by area, area name (each area appears under its STP name and its ICB name), BNF chemical substance or product and BNF section, January 2021 to December 2025. Made from the raw monthly files by `code/01_aggregate_pca_data.ipynb`. The 42-area profile sums over area name and section for each area code; six item names carry a leading or trailing space in the source data, which `code/analyze_af.py` removes (2,161 names, 2,160 items).
- `icb_substance_items_by_year.csv` - columns AREA_CODE, DRUG, YEAR, ITEMS: dispensed items by area, item and calendar year. Made by `code/01_aggregate_pca_data.ipynb`.
- `icb_boundaries_small.geojson` - the 42 ICB polygons (ONS Integrated Care Boards (April 2023) EN BFE, WGS84), simplified by `code/simplify_boundaries.py` (coordinates rounded to 3 decimal places) from the file downloaded by `code/02_fetch_boundaries_and_maps.ipynb`.
- `icb_covariates.csv` - need covariates for the 42 ICBs from OHID Fingertips (area type 221, persons): ods and ons codes; % of the registered population in seven age bands (age_0_4, age_5_14, age_u18, age_15plus, age_65plus, age_75plus, age_85plus; mean of 2023 to 2025, for age_15plus of 2023 and 2025); IMD 2025 score (imd2025); IMD 2019 score (imd2019, 29 ICBs only; reported in Table S1 and compared with IMD 2025, not used in the models); and 21 QOF recorded prevalences (qof_*; mean of 2021/22 to 2024/25 for most). Extracted with `code/fetch_covariates.js`.
- `icb_populations.csv` - registered population by ICB and year (pop2021 to pop2025, denominator of Fingertips indicator 93468; the analysis uses the mean of 2023 to 2025) and QOF diabetes register size by financial year (dm2021_22 to dm2024_25; not used in the published analyses). Extracted with `code/fetch_populations.js`.
- `temporal_data.csv` - columns Period, STP_CODE, Total_Items, Total_Cost, Unique_Chemicals: quarterly items, net ingredient cost (GBP) and number of distinct chemical substances (largest monthly count in the quarter) by area. Made from the raw monthly files by `code/quarterly_totals.py`; its item totals agree exactly with `icb_substance_items_by_year.csv` for every area and year.

## Results (produced by the code)
- `mantel_af.csv` - per-item Mantel r, 1,000-permutation P, national items and top-area share.
- `af_summary.json` - headline numbers of the per-item Mantel test on the 42 areas.
- `temporal_af.csv` - per-item pattern stability (2021 vs 2025), between-area CV in 2021 and 2025 and their ratio, and national growth.
- `mantel_af_fdr.csv` - per-item Mantel results with P-values refined by 99,999 permutations for nominally significant items, Benjamini-Hochberg q-values and the appliance flag.
- `loo_robustness.csv` - Mantel r recomputed against leave-one-out distance matrices for ten named agents and the 20 items with the highest mean within-area share.
- `icb_order.csv` - ICB codes, ONS codes and names in analysis order.
- `moran_af.csv` - per-item Moran's I (contiguity and k-nearest-neighbour weights), permutation P, q-value, nominal and clustered flags (positive I), Mantel results, appliance flag, number of areas with any use, and the P90/P10 amplitude ratio.
- `geo_summary.json` - Moran summary counts, including the k-nearest-neighbour comparison of Table S2 (853 nominal items, rank correlation with the contiguity analysis 0.93), and an unadjusted test of overall prescribing dissimilarity against distance for all items (the similarity analyses reported in the paper are in `need_similarity.json`).
- `moran_by_year.csv` - Moran's I and P by calendar year for named agents.
- `extras_summary.json` - Moran's I screen by calendar year, P90/P10 amplitude summary and the coding-stable sensitivity analysis.
- `access_choice_unadjusted.json` - means, ranges and unadjusted Moran's I of the access-versus-choice measures by year (key `by_year`) and the Moran values of Table S3 (key `tableS3_moran`).
- `trends_summary.json` - national yearly totals and trend tests (Table S5).
- `summary_numbers.txt` - printed output of `code/summary_numbers.py`.
- `need_adjusted_moran.csv` - per-item R2 of the need models, need-adjusted Moran's I, Freedman-Lane P, q-value and nominal and clustered flags for the primary model, the age and deprivation model and the six-component model, plus the primary model with k-nearest-neighbour weights.
- `need_summary.json` - summary of the need-adjusted item-level analysis (variance explained by the QOF components, expected I, counts).
- `need_similarity.json` - similarity of neighbouring areas and rank correlation with distance, unadjusted and need-adjusted, excluding dispensing appliances (primary), for all items, and for items dispensed in every area.
- `placebo_summary.json`, `placebo_runs.csv` - calibration with placebo covariates (summary and one row per placebo set).
- `signflip_check.json`, `signflip_check.csv` - sign-flip check of the need-adjusted item-level tests.
- `class_need.json` - access versus choice by year, unadjusted and need-adjusted (Table 2, Figure 2).
- `class_robustness.json` - robustness of the need-adjusted clustering of dapagliflozin's share (Table S8), standardised coefficients and q-values.
- `section_access.csv`, `section_choice.csv`, `section_choice_summary.json` - BNF section totals and item shares within sections, unadjusted and need-adjusted.
- `lisa_2025.csv`, `lisa_summary.json` - local Moran's I for SGLT2 inhibitor use and dapagliflozin's share in 2025 (Figure 3).
- `tableS6_items.csv`, `tableS7_sections.csv` - the items and sections clustered after adjustment (Tables S6 and S7).
- `need_summary_numbers.txt`, `need_summary_numbers.json` - other numbers quoted in the paper (printed output of `code/need_summary_numbers.py` and the numbers it saves).
- `tableS1.md`, `tableS6.md`, `tableS7.md`, `tableS8.md` - supplementary tables as written by `code/supplementary_tables.py`.
- `mantel_84.csv` - per-item Mantel r and 1,000-permutation P on the 84 area-period profiles (`code/mantel_84.py`; item names as in the source data, 2,161 items).
- `tableS3_items.csv`, `tableS4_agents.csv` - numeric columns of Tables S3 and S4 (`code/phylogeography_tables.py`).
- `spatial_weights_contiguity.csv`, `spatial_weights_knn4.csv` - the contiguity (104 neighbouring pairs) and k-nearest-neighbour (k=4) weights as binary 42 x 42 tables labelled with ODS codes (row i of the k-nearest-neighbour table marks the four nearest areas to area i); `code/spatial_weights.py` row-standardises them for the analysis.

## Sources and licences
See also `LICENSE.md` in this folder.

Community dispensing: NHS Business Services Authority, Prescription Cost Analysis monthly data, https://opendata.nhsbsa.net/dataset/prescription-cost-analysis-pca-monthly-data (Open Government Licence v3.0). Each item is attributed to the ICB of the dispensing contractor. Need covariates and registered populations: Office for Health Improvement and Disparities, Fingertips, https://fingertips.phe.org.uk; contains public sector information licensed under the Open Government Licence v3.0. ICB boundaries: Source: Office for National Statistics licensed under the Open Government Licence v.3.0. Contains OS data © Crown copyright and database right 2023.
