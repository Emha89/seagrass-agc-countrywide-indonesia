# Country-Wide Seagrass Above-Ground Carbon Estimation

R and Google Earth Engine scripts for:

Hafizt, M., Roelfsema, C., Phinn, S., Lyons, M., Wicaksono, P., Salsabila, H.N.,
McMahon, K., Hernawan, U.E., Irawan, A., Adi, N.S., & Rahmawati, S. (in
preparation). **From Area-Based to Pixel-Based: A Country-Wide Approach for
Seagrass Above-Ground Carbon Estimation in Tropical Coastal Ecosystems**.

Companion repositories for the earlier chapters in this thesis:
- Chapter 3 (field data harmonization):
  [seagrass-agc-fieldharmonization-indonesia](https://github.com/Emha89/seagrass-agc-fieldharmonization-indonesia)
- Chapter 4, Study 1 (GSE spectral percent-cover):
  [seagrass-spc-gse-indonesia](https://github.com/Emha89/seagrass-spc-gse-indonesia)
- Chapter 4, Study 2 (proxy-based AGC, 6 regions):
  [seagrass-agc-proxy-indonesia](https://github.com/Emha89/seagrass-agc-proxy-indonesia)

## Overview

This repository extends the proxy-based AGC framework from Chapter 4 (Study
2) to country-wide coverage across Indonesia, and compares two ways of
estimating total AGC: an **area-based** approach (persistence area x mean
density) and a **pixel-based** approach (summing per-pixel model
predictions). Uncertainty is quantified using a cluster design effect (DEFF)
/ intraclass correlation (ICC) framework, reflecting the compiled field data's
survey-cluster structure -- a different approach from Study 2, which retains
a simpler patch-size-based correction (L=30m) at regional scale.

STATUS: R-side script review complete (25 files). GEE-side review in
progress (11 files so far, covering the core national deployment and
the Table S3 per-ROI comparison).

## Repository structure

```
seagrass-agc-countrywide-indonesia/
├── R/
│   ├── func_reg_rf.R                    shared regression utilities
│   ├── func_class_rf.R                  shared classification utilities
│   ├── 01_build_master_raw.R            builds master_raw_<year>.csv
│   ├── A1_rf_model_PA.R                 PA model
│   ├── A2_rf_applyModel_PA.R
│   ├── B1_rf_model_leafMorpho.R         morphology model
│   ├── B2_rf_applyModel_leafMorpho.R
│   ├── B3_rf_evalModel_leafMorpho.R
│   ├── C1_rf_model_cIndex.R             carbon index model
│   ├── C2_rf_applyModel_cIndex.R
│   ├── C3_rf_evalModel_cIndex.R
│   ├── D1_rf_model_SPC.R                percent cover model
│   ├── D2_rf_applyModel_SPC.R
│   ├── D3_rf_evalModel_SPC.R
│   ├── 04_build_master_AGB_AGC.R        builds master_AGB/AGC_<year>.csv
│   │                                     (run twice -- see execution order)
│   ├── E1_rf_model_AGB.R                AGB model
│   ├── E2_rf_applyModel_AGB.R
│   ├── E3_rf_evalModel_AGB.R
│   ├── F1_rf_model_AGC.R                final AGC model
│   ├── F2_rf_applyModel_AGC.R
│   ├── F3_rf_evalModel_AGC.R
│   ├── 22_compile_training_GEE.R        compiles the full GEE training asset
│   ├── calc_SDdata_AGC.R                field-uncertainty scalar + Method B
│   ├── recompute_N_chapter5_DEFF.R      adopted uncertainty method (DEFF/ICC)
│   └── check_DEFF_by_year_chapter5.R    validity check for pooled DEFF
├── gee/
│   ├── 01_trainingData_prep_indo.js       extracts GSE+depth at training points
│   ├── 02_PROB_rf_modelDev_indo.js        PA model
│   ├── 02b_seagrassExtent_indo.js         persistence mask (from 02's output)
│   ├── 03_MORPHO_rf_modelDev_indo.js      morphology model
│   ├── 04_COVER_rf_modelDev_indo.js       SPC model (two-step workflow)
│   ├── 05_AGB_rf_modelDev_indo.js         AGB model
│   ├── 06_cIndex_rf_modelDev_indo.js      carbon index model
│   ├── 07_AGC_rf_modelDev_indo.js         final AGC model (two-step workflow)
│   ├── 08_AGC_rf_modelUncer_indo.js       national total AGC + DEFF/ICC CI
│   ├── export_ROI_collection.js           builds the named-ROI asset
│   ├── 09_AGC_rf_compare_indo.js          per-ROI Method A vs B (Table S3)
│   └── GEE_Apps/
│       └── 09_AGC_timeSeries_forecast_app.js   live App (see below)
└── data/
    └── field_data_template.xlsx   column headers only, see Data availability below
```

STATUS: GEE-side review complete for the core national deployment
(01-08), the Table S3 per-ROI comparison, and the App. Remaining
scripts from the original inventory not reviewed: 07_AGC_rf_calculationC
(a possible v9/v10 variant) and any indoMPA/indoROI-specific variants of
the model-uncertainty script beyond the national one already included.

## Suggested execution order

Similar to the Study 2 repository, several stages here depend on
outputs from later-numbered scripts, and one script needs to run
twice:

1. `01_build_master_raw.R` -- produces `master_raw_<year>.csv`
2. `A1_rf_model_PA.R` then `A2_rf_applyModel_PA.R` -- produces
   `predicted_PA_<year>.csv`
3. `B1_rf_model_leafMorpho.R` then `B2_rf_applyModel_leafMorpho.R` --
   produces `predicted_MORPH3_probs_<year>.csv` (needs step 2's output)
4. `C1_rf_model_cIndex.R` then `C2_rf_applyModel_cIndex.R` -- needs
   only step 2's output; produces `predicted_CINDEX_<year>.csv`
5. `D1_rf_model_SPC.R` then `D2_rf_applyModel_SPC.R` -- needs step 3's
   morphology probabilities; produces `predicted_tSPC_<year>.csv`
6. `04_build_master_AGB_AGC.R` -- **pass 1**: produces
   `master_AGB_<year>.csv` (needs step 5's output). Also writes
   `master_AGC_<year>.csv` at this point, but with `tAGB_pred` all NA
   since step 7 hasn't run yet -- not usable for training until pass 2.
7. `E1_rf_model_AGB.R` then `E2_rf_applyModel_AGB.R` -- produces
   `predicted_tAGB_<year>.csv`
8. `04_build_master_AGB_AGC.R` again -- **pass 2**: re-run so that
   `master_AGC_<year>.csv` picks up step 7's `tAGB_pred` values
9. `F1_rf_model_AGC.R`, `F2_rf_applyModel_AGC.R`, `F3_rf_evalModel_AGC.R`
10. Any `*_evalModel_*.R` script (B3, C3, D3, E3) can run any time
    after its matching apply-model step above.
11. `22_compile_training_GEE.R` -- needs every stage's `training_*.csv`
    to exist (i.e. run after step 9)
12. `calc_SDdata_AGC.R`, `recompute_N_chapter5_DEFF.R`, and
    `check_DEFF_by_year_chapter5.R` -- all need `training_AGC.csv`
    (from step 9); the latter two are the adopted uncertainty method
    and its per-year validity check respectively.

**Confirmed predictor set for the final AGC model** (a meaningful
difference from Study 2): `F1_rf_model_AGC.R` uses only the 64 GSE
bands plus `tAGB_pred`, `carbon_index`, and `tSPC_pred` as predictors --
no depth, no PA_prob, no morphology probabilities. See that script's
own header for details.

## Data availability

This repository contains **analysis code only**. Field survey data compiled
from third-party sources are subject to the data-sharing policies of the
original providers and are not included here (see the paper's Data
Availability statement once published).

**Field data template**: `data/field_data_template.xlsx` documents the
exact column headers the R and GEE scripts expect, across three sheets:
`GT_Label_Data` (identity, PA, tSPC, morphology, AGB, AGC -- read by
`01_build_master_raw.R`), `Species_Composition` (per-species percent
cover, joined on `compositio` -- read by the same script), and
`GSE_Training_PerYear` (reference only: the GSE band + depth columns
`01_trainingData_prep_indo.js` extracts automatically, not filled in by
hand). Raw data is not included; only the template.

**GEE asset paths**: all Earth Engine asset paths in `gee/*.js` are
placeholders (`ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER'`).
The underlying assets themselves (training points, bathymetry, model
outputs, etc.) are private to the corresponding author's Earth Engine
account and are not made public -- this mirrors how field survey data is
handled on the R side. The code documents the full workflow and is
reusable with an equivalent set of assets under your own GEE project;
it is not runnable as-is without them.

## Exploratory work not included

Two supplementary analyses were carried out locally but are not part of
this repository, because neither directly produces a result reported in
the paper:

- **Spatial correlation length (L) via variogram**: an earlier attempt
  at the national uncertainty correction, producing unstable/implausible
  fits (spurious "hump" shapes, fitted ranges of hundreds of km) because
  the training data is a compilation of discrete field campaigns, not a
  spatially continuous survey. Replaced entirely by the cluster design
  effect (DEFF/ICC) approach in `recompute_N_chapter5_DEFF.R`.
- **Pixel-based vs area-based ROI selection**: exploratory analysis used
  to choose which ROIs to feature when comparing the two estimation
  methods for the paper -- an internal decision-making step, not itself
  a reported result.
- **Interactive AGC/Sentinel-2 validation viewer**: a personal visual
  QA tool (adjustable PA-probability slider over a Sentinel-2 true-colour
  composite) used during development, not a script that produces a
  reported figure or table.
- **Earlier per-ROI Method A vs B comparison script**: superseded by
  09_AGC_rf_compare_indo.js. The earlier version used the pre-DEFF/ICC,
  L=30-based uncertainty formula and didn't explicitly target Table S3;
  the current version removes that (locally non-representative)
  uncertainty calculation and formats its output directly as Table S3
  rows.
- **Alternate GEE Apps, superseded or not in current use**: an earlier
  App revision (v8r2, superseded by the App included here), an
  MPA-dropdown App with live Monte Carlo ensemble training (a different,
  more expensive design than the pre-computed-asset Apps kept here), and
  a separate ROI-dropdown comparison App -- none of these are the
  currently published Apps.

## Interactive AGC Viewer

A live App lets you draw a polygon anywhere in Indonesia and see annual
AGC totals (2017-2024, with DEFF/ICC confidence intervals) plus a
forecast to 2027, at:
https://muhammadhafizt.users.earthengine.app/view/seagrassagcmpaindonsesia

Source: `gee/GEE_Apps/09_AGC_timeSeries_forecast_app.js`.

A second script, used to view every proxy layer (PA, morphology, SPC,
AGB, carbon index, AGC) across all of Indonesia and its 6 study
locations, by year, is available directly in the GEE Code Editor
(requires a Google/GEE account to open):
https://code.earthengine.google.com/308cd0b53b30d84902fe450210c17236

## Citation

_(to be added once the paper is accepted)_

## License

This project is licensed under the MIT License -- see the
[LICENSE](LICENSE) file for details. Field data used to produce the
results are not covered by this license -- see Data availability above.

## Contact

Muhammad Hafizt
School of the Environment, The University of Queensland, Brisbane, Australia
National Research and Innovation Agency of Indonesia (BRIN)
m.hafizt@uq.edu.au