# Geographic structure in English NHS prescribing (2021-2025)

Whole-formulary spatial analysis of national NHS community dispensing data for England, with adjustment for population need, separation of access to a therapeutic class from the choice of agent within it, local indicators of spatial association, and a comparison with a phylogeographic signal method.

This repository accompanies the manuscript *Geographic structure in English NHS prescribing (2021-2025): a need-adjusted spatial analysis of 6 billion prescription items* (Karagol et al.).

## Summary

We analyse every item dispensed in the community in England from January 2021 to December 2025 (NHS Business Services Authority Prescription Cost Analysis, 6.04 billion items), aggregated by Integrated Care Board (ICB) area (n=42) and British National Formulary item (n=2,160). Each item's value in an area is its share of all items dispensed there. Spatial clustering is measured with global Moran's I on ICB contiguity, before and after regressing each item's share on ICB-level need covariates (age structure, Index of Multiple Deprivation 2025, and Quality and Outcomes Framework prevalence), with Freedman-Lane permutation tests for the regression residuals.

Key findings:

- Excluding dispensing appliances, neighbouring areas prescribe more alike than other areas (Jensen-Shannon distance 16.1% lower). Adjustment for age and deprivation reduces this to 7.2%, and further adjustment for recorded morbidity to 3.6%, which remains significant.
- 577 of 2,160 items (26.7%) are spatially clustered before adjustment (positive Moran's I, false-discovery-rate q<0.05), covering 61% of all items dispensed; 120 (5.6%) remain clustered after adjustment, covering 10.5%. Placebo covariates with the same spatial structure as the need covariates leave a mean of 284 items clustered (with the same procedure, the real covariates leave 110).
- The increasingly regional pattern in the share of dispensing taken by SGLT2 inhibitors (Moran's I 0.17 in 2021, 0.40 in 2025) is accounted for by the need covariates, mainly age. The choice between SGLT2 inhibitors is not: dapagliflozin's share of the class ranged from 37% to 83% between areas in 2025 and remained spatially autocorrelated after adjustment (Moran's I 0.18, P=0.008; q=0.10 across the 30 need-adjusted class-year tests), consistently from 2023 to 2025, with a low-share block across Yorkshire, Greater Manchester and the East Midlands.
- The phylogeographic (Mantel) signal method, which does not use geography, largely misses this structure; its strongest signals are dispensing appliances and very rarely dispensed items.

## Repository layout

    .
    ├── run_all.sh                        runs steps 1 to 24 in order
    ├── code/
    │   ├── method.py                     Jensen-Shannon distance and Mantel test
    │   ├── need_adjustment_core.py       shared functions for the need-adjusted analyses
    │   ├── analyze_af.py                 step 1: 42-area profiles, distance matrix, per-item Mantel test
    │   ├── temporal_af.py                step 2: year-resolved shares, pattern stability, dispersion
    │   ├── fdr_refine.py                 step 3: 99,999-permutation Mantel P-values, FDR, appliance flag
    │   ├── robustness_loo.py             step 4: leave-one-out distance-matrix check
    │   ├── mantel_84.py                  step 5: Mantel test on the 84 area-period profiles
    │   ├── spatial_weights.py            step 6: contiguity and k-nearest-neighbour weights, centroid distances
    │   ├── spatial_analysis.py           step 7: unadjusted Moran's I for every item
    │   ├── spatial_sensitivity.py        step 8: Moran's I by year, amplitude, coding-stable items
    │   ├── trends.py                     step 9: national trends (Table S5)
    │   ├── make_geography_figures.py     step 10: Figure S1; access-versus-choice means and Table S3 Moran values
    │   ├── make_figures.py               step 11: Figures S3 to S7
    │   ├── summary_numbers.py            step 12: numbers of the unadjusted spatial and phylogeographic analyses
    │   ├── phylogeography_tables.py      step 13: numeric columns of Tables S3 and S4
    │   ├── need_adjustment.py            step 14: need-adjusted Moran's I for every item (Table 1)
    │   ├── need_similarity.py            step 15: similarity of neighbouring areas before and after adjustment
    │   ├── placebo_covariates.py         step 16: calibration with placebo covariates
    │   ├── signflip_check.py             step 17: sign-flip check of the clustered items
    │   ├── class_need.py                 step 18: access versus choice by year, unadjusted and need-adjusted (Table 2)
    │   ├── class_robustness.py           step 19: robustness of the choice result, multiplicity (Table S8)
    │   ├── section_choice.py             step 20: BNF sections and within-section shares
    │   ├── lisa.py                       step 21: local Moran's I for SGLT2 inhibitors (Figure 3)
    │   ├── need_summary_numbers.py       step 22: other numbers quoted in the paper
    │   ├── supplementary_tables.py       step 23: Tables S1, S6, S7 and S8 (markdown)
    │   ├── make_need_figures.py          step 24: Figures 1 to 3, S2 and S8
    │   ├── 01_aggregate_pca_data.ipynb   Colab: raw monthly files to the two area aggregates
    │   ├── quarterly_totals.py           raw monthly files to data/temporal_data.csv
    │   ├── 02_fetch_boundaries_and_maps.ipynb  Colab: ICB boundaries
    │   ├── simplify_boundaries.py        full ONS boundaries to data/icb_boundaries_small.geojson
    │   ├── fetch_covariates.js           browser script that extracts data/icb_covariates.csv from Fingertips
    │   └── fetch_populations.js          browser script that extracts data/icb_populations.csv from Fingertips
    ├── data/                             inputs and results (see data/README.md and data/LICENSE.md)
    ├── figures/                          Figures 1 to 3 and S1 to S8 (PNG)
    ├── requirements.txt                  minimum versions
    ├── requirements-pinned.txt           exact versions
    ├── LICENSE
    └── README.md

## Reproducing the analysis

The aggregated data and covariates needed to reproduce every number, table and figure are included in `data/`. Install the dependencies with `pip install -r requirements.txt` and run

    bash run_all.sh

or the steps one by one:

    python code/analyze_af.py
    python code/temporal_af.py
    python code/fdr_refine.py            # about 7 minutes
    python code/robustness_loo.py
    python code/mantel_84.py             # about 2 minutes
    python code/spatial_weights.py
    python code/spatial_analysis.py      # about 2 minutes
    python code/spatial_sensitivity.py
    python code/trends.py
    python code/make_geography_figures.py
    python code/make_figures.py
    python code/summary_numbers.py > data/summary_numbers.txt
    python code/phylogeography_tables.py
    python code/need_adjustment.py       # about 4 minutes
    python code/need_similarity.py       # about 30 minutes
    python code/placebo_covariates.py    # 75 to 100 minutes
    python code/signflip_check.py
    python code/class_need.py
    python code/class_robustness.py
    python code/section_choice.py        # about 3 minutes
    python code/lisa.py
    python code/need_summary_numbers.py > data/need_summary_numbers.txt
    python code/supplementary_tables.py
    python code/make_need_figures.py

The whole run takes about 2 to 2.5 hours on one core. Each script sets its own working directory to `data/`, so it can be run from anywhere; figures are written to `figures/`. Steps 1, 2 and 6 write intermediate arrays (`.npy`, `.pkl`) into `data/` that later steps use, so run the steps in order (step 15 also writes the three distance matrices behind Figure 1A as `.npy` files). All permutation tests use fixed seeds, so a rerun reproduces the results in `data/` exactly. For a quicker run (about 50 minutes), `PLACEBO_ARGS="10 5 5" bash run_all.sh` uses 10, 5 and 5 placebo sets instead of 200, 100 and 100; this overwrites the placebo results, the placebo numbers in `data/need_summary_numbers.txt` and Figure S8, which `git checkout -- data/placebo_summary.json data/placebo_runs.csv data/need_summary_numbers.txt figures/figS8_placebo.png` restores.

Tested with Python 3.11 and numpy 2.4, pandas 3.0, scipy 1.17, matplotlib 3.10, geopandas 1.1 (with pyogrio 0.13) and shapely 2.1; `requirements-pinned.txt` lists the exact versions.

## Where each table and figure comes from

| Item | Scripts | Result files in `data/` |
|---|---|---|
| Table 1 | spatial_analysis.py, need_adjustment.py, need_similarity.py, placebo_covariates.py, section_choice.py, need_summary_numbers.py | moran_af.csv, need_adjusted_moran.csv, need_summary.json, need_similarity.json, placebo_summary.json, section_choice_summary.json, need_summary_numbers.txt |
| Table 2, Figure 2 | class_need.py, make_need_figures.py | class_need.json |
| Figure 1 | need_similarity.py, need_adjustment.py, make_need_figures.py | need_similarity.json, need_adjusted_moran.csv |
| Figure 3 | lisa.py, make_need_figures.py | lisa_2025.csv, lisa_summary.json |
| Table S1 | supplementary_tables.py | tableS1.md |
| Table S2 | analyze_af.py, fdr_refine.py, mantel_84.py, spatial_analysis.py, spatial_sensitivity.py, summary_numbers.py | mantel_af_fdr.csv, mantel_84.csv, moran_af.csv, geo_summary.json, extras_summary.json, summary_numbers.txt |
| Table S3 | phylogeography_tables.py | tableS3_items.csv |
| Table S4 | phylogeography_tables.py | tableS4_agents.csv |
| Table S5, Figure S7 | trends.py, make_figures.py | trends_summary.json |
| Table S6 | need_adjustment.py, need_summary_numbers.py, supplementary_tables.py | tableS6_items.csv, tableS6.md |
| Table S7 | section_choice.py, need_summary_numbers.py, supplementary_tables.py | section_access.csv, tableS7_sections.csv, tableS7.md |
| Table S8 | class_need.py, class_robustness.py, supplementary_tables.py | class_need.json, class_robustness.json, tableS8.md |
| Figure S1 | make_geography_figures.py | moran_af.csv, icb_boundaries_small.geojson and the arrays written by steps 1, 2 and 6 |
| Figures S2 and S8 | make_need_figures.py | moran_af.csv, placebo_runs.csv, placebo_summary.json |
| Figures S3 to S6 | make_figures.py | mantel_af.csv, mantel_84.csv and the arrays written by steps 1 and 2 |

## Where the input files come from

- `icb_area_substance_items_full.csv` and `icb_substance_items_by_year.csv`: the raw monthly Prescription Cost Analysis files (about 15 GB, not included) aggregated by `code/01_aggregate_pca_data.ipynb`.
- `temporal_data.csv`: quarterly items, net ingredient cost and number of chemical substances by area, made from the same raw files by `code/quarterly_totals.py`. Its item totals agree exactly, area by area and year by year, with `icb_substance_items_by_year.csv`.
- `icb_boundaries_small.geojson`: the ONS Integrated Care Boards (April 2023) boundaries downloaded by `code/02_fetch_boundaries_and_maps.ipynb` (SHA-256 86b8ee43371839c2af75be45bf2418d867ef10831d1a855c79aa8de3877c422a) and simplified by `code/simplify_boundaries.py`, which gives this file byte for byte.
- `icb_covariates.csv` and `icb_populations.csv`: extracted from OHID Fingertips on 1 October 2026 with the two browser scripts in `code/`, which give these files byte for byte (SHA-256 fe78a076dbf3c5ad440502b32f64265fee96bbaf206cb4f7caf9c548fca2baca and 5dfe6d1f875c903201d89c802ad4f7f87f3eb45c1a337a3b49b83f3c70c6bb78; OHID may revise the source data).

## Method

Area profiles are normalised item distributions over 2,160 BNF chemical and product items. Pairwise distance is the square root of the Jensen-Shannon divergence (natural logarithm, smoothing constant 1e-10). Comparisons of overall profiles exclude the 49 dispensing appliance items, whose dispensing is attributed to the few areas that host dispensing appliance contractors.

Spatial clustering: ICBs are neighbours if their boundaries touch (queen contiguity, April 2023 boundaries on the British National Grid, 100 m tolerance; 104 neighbouring pairs), with row-standardised weights and k-nearest-neighbour weights (k=4) as a sensitivity analysis. For each item, global Moran's I of the within-area share is tested one-sided with 9,999 permutations, refined to 99,999 for items with P<0.05, P = (b+1)/(N+1), and the Benjamini-Hochberg false discovery rate across all 2,160 items. An item is classed as spatially clustered when I is positive and q<0.05. Permuted statistics within 1e-12 of the observed value count as ties.

Need adjustment: each item's share is regressed on ICB need covariates by ordinary least squares (primary model: % aged 65+, % aged under 18, IMD 2025 score and the first two principal components of 21 QOF prevalences). Moran's I of the residuals is tested with the Freedman-Lane scheme (residuals permuted and re-projected through the residual-maker matrix). Need-adjusted overall profiles are the mean profile plus each area's residuals (negative values set to zero, rows renormalised), tested with a Freedman-Lane-type permutation of the residual rows. Placebo covariates are generated by Moran spectral randomisation (pair method), preserving the covariance and approximately the spatial autocorrelation of the real covariates.

Access versus choice: for SGLT2 inhibitors, GLP-1 and GIP/GLP-1 receptor agonists, and thiazide-type diuretics, the class's total within-area share (use of the class) is separated from one agent's share of the class (choice of agent), each tested with Moran's I in every year, unadjusted and after adjustment with a class-specific need model. Local Moran's I (conditional permutation, folded pseudo P-values) locates the regional blocks.

Phylogeographic signal: each item's pairwise absolute-difference matrix is correlated with the Jensen-Shannon matrix (Mantel r), with 1,000 area-label permutations for the screen and 99,999 for items with P<0.05, followed by the Benjamini-Hochberg procedure. The approach does not use the locations of the areas.

## Data sources and licences

- Community dispensing: NHS Business Services Authority, Prescription Cost Analysis (PCA) monthly data, https://opendata.nhsbsa.net/dataset/prescription-cost-analysis-pca-monthly-data. Contains public sector information licensed under the Open Government Licence v3.0. Each item is attributed to the ICB of the dispensing contractor, not the prescriber.
- Need covariates and registered populations: Office for Health Improvement and Disparities, public health profiles (Fingertips), https://fingertips.phe.org.uk, area type 221. Contains public sector information licensed under the Open Government Licence v3.0.
- ICB boundaries: Office for National Statistics Open Geography Portal, Integrated Care Boards (April 2023) EN BFE. Source: Office for National Statistics licensed under the Open Government Licence v.3.0. Contains OS data © Crown copyright and database right 2023.

No patient-level data are used.

## Citation

Citation to be added on publication.

## License

Code is released under the MIT License (see `LICENSE`). The data in `data/` are derived from the sources above and remain under the Open Government Licence v3.0 (see `data/LICENSE.md`).
