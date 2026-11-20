# What pixel-based maps add to area-based estimates of seagrass above-ground carbon across Indonesia

R and Google Earth Engine (GEE) scripts for:

Hafizt, M., McMahon, K., Roelfsema, C., Hernawan, U.E., Adi, N.S., Lyons, M.,
Phinn, S. (in preparation). **What pixel-based maps add to area-based estimates
of seagrass above-ground carbon: a national test across Indonesia.**

The scripts map seagrass above-ground carbon (AGC, gC m⁻²) at 10 m across
Indonesia for 2017–2024 with a random forest on Google Satellite Embeddings
(GSE, 64 bands). On one fixed seagrass extent they compare the national and
ecoregion totals of three methods:

- **Pixel-based GSE**: the sum of the mapped AGC over all extent pixels.
- **Tier 1**: extent area × the national field mean AGC.
- **Tier 2**: area of each marine ecoregion (MEOW) × its field mean AGC, or the
  national mean where an ecoregion has fewer than 30 samples.

The methods are also tested on the same held-out samples, with two benchmarks:
the **local mean** of the sampled meadow, and a **Sentinel-2** model in the same
annual format as GSE. Tests run at the sample level under five validation
designs, and at the inventory level for held-out sites, meadows and
ecoregions:

| Design | Held-out unit | Distance to the nearest training sample |
|---|---|---|
| `pixel` | 10 m pixel | ≤ 10 m |
| `b100` | 100 m block (site) | 10–100 m |
| `b1000` | 1 km block | 100 m–1 km |
| `cluster` | meadow (DBSCAN cluster, 5 km) | 1–20 km |
| `ecoregion` | marine ecoregion (MEOW) | > 20 km |

## Repository structure

```
seagrass-agc-countrywide-indonesia/
├── R/
│   ├── P4_00_config_functions.R    configuration (edit the CONFIG block) and shared functions
│   ├── P4_01_build_training.R      sample table, validation groups, ecoregions, DBSCAN, Sentinel-2 join
│   ├── P4_01a_link_measured_agc.R  optional diagnostic: links field-measured AGC (not used for the results)
│   ├── P4_01b_trace_actual_agb.R   optional diagnostic: traces the measured AGB (not used for the results)
│   ├── P4_02_tuning.R              random forest settings (pixel-grouped CV)
│   ├── P4_03_validation.R          out-of-fold predictions, five designs, all methods (long run)
│   ├── P4_04_metrics.R             sample-level and inventory-level accuracy
│   ├── P4_05_tier_uncertainty.R    design effect, effective sample size, Tier 1 and Tier 2
│   ├── P4_06_gee_prepare.R         training table and ecoregion shapefile for the GEE upload
│   ├── P4_07_gee_results.R         national and ecoregion totals from the GEE exports
│   ├── P4_08_figures.R             figures of the paper and the Supplement
│   ├── P4_08b_qgis_layers.R        GIS layers for drawing the study area map in QGIS
│   ├── P4_09_tables.R              tables of the paper and the Supplement
│   └── P4_11_zenodo_package.R      open data deposit from processed results only
├── GEE/
│   ├── README.md                   role, inputs and outputs of every GEE script
│   ├── P4_G00_extract_gse.js       GSE values at the samples (current dataset version)
│   ├── P4_G01_train_apply.js       AGC and cover models (cover unused in v3), core extent, annual maps
│   ├── P4_G02_totals.js            area and AGC sums by year and ecoregion
│   ├── P4_G03_extent_coverage.js   samples inside the extent; distance to the nearest AGC sample
│   ├── P4_G04_temporal.js          between-year CV of AGC; GSE similarity between years
│   ├── P4_G05_extent_sensitivity.js  extent under three geomorphic mask options
│   ├── P4_G06_addon_maps.js        add-on extent and maps (see below)
│   ├── P4_G07_addon_tables.js      sums of G02, G03 and G04 for the add-on extent
│   ├── P4_G08_figure_grid.js       ~500 m cell sums for the national map figures
│   ├── P4_G08b_zoom_windows.js     10 m zoom windows (Rote, Derawan, Seribu Islands, Selayar)
│   ├── P4_G09_pa_cv.js             cross-validation of the presence model at three thresholds
│   ├── P4_G10_export_zenodo.js     GeoTIFF export of the extent and the AGC maps (open data)
│   └── P4_G11_extract_s2.js        Sentinel-2 L2A predictors at the samples (benchmark model)
└── data/
    ├── README.md                   input format and folders
    └── master_raw_template.csv     column headers of the field input table (no data)
```

## Requirements

- R 4.3 or later with `dplyr`, `readr`, `tidyr`, `purrr`, `ggplot2`, `scales`,
  `tibble`, `ranger`, `sf`, `dbscan` and `FNN`; the figures also need `terra`,
  `patchwork` (1.1 or later) and `maps`; `P4_11` uses `zip` if installed.
- A Google Earth Engine account and a Cloud project for the GEE scripts.

## Setup

1. Open the `.Rproj` file in RStudio. Every script runs from the repository
   root, for example `source("R/P4_01_build_training.R")`, and every path in
   `CFG` is relative to the root.
2. Edit the CONFIG block of `R/P4_00_config_functions.R` if your folders differ.
3. In every GEE script, set `ASSET_ROOT` to your own asset folder
   (placeholder: `projects/YOUR-GEE-PROJECT/assets`).
4. Place the inputs listed in `data/README.md`. Field data are not shared.

## Run order

R and GEE steps alternate. Tables exported by GEE to Google Drive are copied
to `data/gee/`.

| Step | Script | What it does |
|---|---|---|
| 1 | `R/P4_01_build_training.R`, `R/P4_06_gee_prepare.R` | Sample table and validation groups; GEE upload table and ecoregion shapefile |
| 2 | `GEE/P4_G00_extract_gse.js` | Re-extracts the 64 GSE values at every sample from the current dataset version (set `gse_path` in `P4_00` and rerun step 1) |
| 3 | `GEE/P4_G11_extract_s2.js` | Sentinel-2 L2A percentile predictors at every AGC sample, built like the GSE annual layer (rerun `P4_01` to join them) |
| 4 | `R/P4_02_tuning.R` | Random forest settings; copy them into `P4_G01` |
| 5 | `GEE/P4_G01_train_apply.js`, then `P4_G06` and `P4_G07` | Annual AGC maps on the core extent, then the add-on extent |
| 6 | `GEE/P4_G02` to `P4_G05`, `P4_G08`, `P4_G08b`, `P4_G09` | Sums, distance to samples, temporal checks, mask options, figure grids, presence model CV |
| 7 | `R/P4_03_validation.R` | Leak-free cross-validation of every method (2–3.5 h) |
| 8 | `R/P4_04_metrics.R`, `R/P4_05_tier_uncertainty.R`, `R/P4_07_gee_results.R` | Accuracy, uncertainty, national and ecoregion totals |
| 9 | `R/P4_08_figures.R`, `R/P4_08b_qgis_layers.R`, `R/P4_09_tables.R` | Figures, GIS layers and tables |
| 10 | `GEE/P4_G10_export_zenodo.js`, then `R/P4_11_zenodo_package.R` | Open data deposit (see Data availability) |

`R/P4_01a` and `R/P4_01b` are optional diagnostics and are not used for the
results. Each script prints its version and lists what it reads and writes in
its header. All outputs carry the prefix `P4_`. Figure and table numbers are
set in `FIGNUM` (`P4_08`) and `TABNUM` (`P4_09`).

### Core and add-on extent

The seagrass extent is shallow water (Allen Coral Atlas depth 0–5 m) with a
maximum annual presence probability ≥ 0.7, minus the Allen Coral Atlas
geomorphic class Terrestrial Reef Flat. In `P4_G01` the geomorphic test also
removed shallow pixels that have no geomorphic data. `P4_G06` maps these pixels
as an add-on with the same model, settings and grid, without changing the
`P4_G01` assets. The final extent and maps are core plus add-on, and all
tables are sums of both parts.

## Versions

- **Current (design v3).** Tier 2 uses ecoregion strata, the cover model is
  removed, a Sentinel-2 benchmark and a local mean benchmark are added, and
  accuracy is also tested for held-out sites, meadows and ecoregions. The
  sample build, folds and seeds are unchanged, so the GSE predictions and maps
  are those of the previous version.
- The thesis version is kept on branch `thesis-chapter5` and tag `thesis-v1`.

## Data availability

Field data are not shared, under Indonesian data-sharing policy. They are
available from the corresponding author on reasonable request, with
permission from the data custodians. The processed results are on Zenodo
(DOI to be added): the seagrass extent, the annual AGC maps for 2017–2024 and
their mean, the tables of the paper and the Supplement, and the extent by
distance class per 0.1° cell. `P4_11_zenodo_package.R` builds this deposit and
stops if a table holds sample IDs, coordinates or sample-level rows.

## Licence

MIT (see `LICENSE`).
