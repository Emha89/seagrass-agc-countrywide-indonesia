# Seagrass above-ground carbon across Indonesia: pixel-based and area-based estimates

R and Google Earth Engine (GEE) scripts for:

Hafizt, M., McMahon, K., Roelfsema, C., Hernawan, U.E., Adi, N.S., Lyons, M.,
Phinn, S. (in preparation). **Pixel-based and area-based estimates of seagrass
above-ground carbon across Indonesia tested on the same held-out samples.**

The scripts map seagrass above-ground carbon (AGC, gC m⁻²) and percent cover at
10 m across Indonesia for 2017–2024 with random forest models on Google
Satellite Embeddings (GSE, 64 bands). On one fixed seagrass extent they compare
three national estimates:

- **Pixel-based**: the sum of the mapped AGC over all extent pixels.
- **Tier 1**: extent area × one field mean AGC.
- **Tier 2**: area of each cover class × the field mean AGC of that class
  (classes < 30 %, 30–59.9 % and ≥ 60 %).

All three methods are tested on the same held-out samples under five
validation designs, from held-out pixels to held-out ecoregions:

| Design | Held-out unit | Distance to the nearest training sample |
|---|---|---|
| `pixel` | 10 m pixel | ≤ 10 m |
| `b100` | 100 m block | 10–100 m |
| `b1000` | 1 km block | 100 m–1 km |
| `cluster` | meadow cluster (DBSCAN, 5 km) | 1–20 km |
| `ecoregion` | marine ecoregion (MEOW) | > 20 km |

## Repository structure

```
seagrass-agc-countrywide-indonesia/
├── R/
│   ├── P4_00_config_functions.R   configuration (edit the CONFIG block) and shared functions
│   ├── P4_01_build_training.R     training tables, validation groups, ecoregions, DBSCAN clusters
│   ├── P4_02_tuning.R             random forest settings (pixel-grouped CV)
│   ├── P4_03_validation.R         out-of-fold predictions for the five designs (long run)
│   ├── P4_04_metrics.R            RMSE, MAE, R², MEC and bias; accuracy against distance
│   ├── P4_05_tier_uncertainty.R   Tier 1 and Tier 2 means, design effect, effective sample size
│   ├── P4_06_gee_prepare.R        training table and ecoregion shapefile for the GEE upload
│   ├── P4_07_gee_results.R        national, ecoregion and ROI tables from the GEE exports
│   ├── P4_08_figures.R            figures of the paper and the supplement
│   └── P4_09_tables.R             tables of the paper and the supplement
├── gee/
│   ├── P4_G00_extract_gse.js      GSE values at the training samples (current dataset version)
│   ├── P4_G01_train_apply.js      AGC and cover models, core extent, annual maps 2017–2024
│   ├── P4_G02_totals.js           area and AGC sums by year, ecoregion, cover class and ROI
│   ├── P4_G03_extent_coverage.js  samples inside the extent; distance to the nearest AGC sample
│   ├── P4_G04_temporal.js         between-year CV of AGC; GSE similarity between years
│   ├── P4_G05_extent_sensitivity.js  extent under three geomorphic mask options
│   ├── P4_G06_addon_maps.js       add-on extent and maps (see below)
│   ├── P4_G07_addon_tables.js     sums of G02, G03 and G04 for the add-on extent
│   ├── P4_G08_figure_grid.js      ~500 m cell sums for the national map figures
│   ├── P4_G08b_zoom_windows.js    10 m zoom windows (Rote, Derawan, Seribu Islands, Selayar)
│   ├── P4_G09_pa_cv.js            cross-validation of the presence model at three thresholds
│   └── P4_G10_export_zenodo.js    GeoTIFF export of the extent and annual maps (open data)
└── data/
    ├── README.md                  input format
    └── master_raw_template.csv    column headers of the field input table (no data)
```

## Requirements

- R ≥ 4.1 with `dplyr`, `readr`, `tidyr`, `purrr`, `ggplot2`, `scales`, `tibble`,
  `ranger`, `sf`, `dbscan` and `FNN`; for the figures also `terra`, `patchwork`,
  `ggrepel` and `maps`.
- A Google Earth Engine account and a Cloud project for the GEE scripts.

## Setup

1. Set the working directory to the `R/` folder (each script has a commented
   `setwd()` line at the top).
2. Edit the CONFIG block of `P4_00_config_functions.R`. The paths are relative
   to `R/` (`../result`, `../result/P4v2`, `../data/MEOW/...`). You can also
   define a list `CFG_USER` with your own values before sourcing the file; it
   overrides the defaults.
3. In every GEE script, set `ASSET_ROOT` to your own asset folder
   (placeholder: `projects/YOUR-GEE-PROJECT/assets`).

## Run order

R and GEE steps alternate. Files exported by GEE to Google Drive are copied to
`<out_dir>/gee/` before the next R step.

1. **Training tables.** `P4_01_build_training.R`, then `P4_06_gee_prepare.R`.
   Upload `P4_gee_training.csv` and `P4_meow_indo.shp` to GEE.
2. **GSE version (optional).** `P4_G00_extract_gse.js` re-extracts the 64 GSE
   values at every sample from the current dataset version. Set `gse_path` in
   `P4_00` to the exported CSV and rerun `P4_01` and `P4_06`.
3. **Tuning.** `P4_02_tuning.R`. Copy the values in `P4_best_params.csv` to
   `RF_AGC` and `RF_COV` in `P4_G01` (`mtry` → `variablesPerSplit`,
   `min_node_size` → `minLeafPopulation`, `sample_fraction` → `bagFraction`).
4. **Maps.** `P4_G01_train_apply.js`. The console prints a fingerprint of the
   training table, which must match the one printed by `P4_06`.
5. **Validation (while GEE runs).** `P4_03_validation.R` (1–2 h), then
   `P4_05_tier_uncertainty.R`.
6. **Add-on extent.** `P4_G06_addon_maps.js`, then `P4_G07_addon_tables.js`.
7. **Sums and checks.** `P4_G02_totals.js`, `P4_G03_extent_coverage.js`,
   `P4_G04_temporal.js` and `P4_G05_extent_sensitivity.js`. Set `roi_path` in
   `P4_00` to `<out_dir>/gee/P4_roi_polygons.geojson`.
8. **Results.** `P4_04_metrics.R`, then `P4_07_gee_results.R`.
9. **Figure inputs.** `P4_G08_figure_grid.js`, `P4_G08b_zoom_windows.js` and
   `P4_G09_pa_cv.js`.
10. **Figures and tables.** `P4_08_figures.R` and `P4_09_tables.R`. Figures are
    printed; uncomment the `ggsave()` lines to write the files.
11. **Open data.** `P4_G10_export_zenodo.js` writes the extent and the annual
    AGC and cover maps as GeoTIFFs.

Each script prints its version and lists what it reads and writes in its
header. All outputs carry the prefix `P4_`.

### Core and add-on extent

The seagrass extent is shallow water (Allen Coral Atlas depth 0–5 m) with a
maximum annual presence probability ≥ 0.7, minus the Allen Coral Atlas
geomorphic class Terrestrial Reef Flat. In `P4_G01` the geomorphic test also
removed shallow pixels that have no geomorphic data. `P4_G06` maps these pixels
as an add-on with the same models, settings and grid, without changing the
`P4_G01` assets. The final extent and maps are core plus add-on; `P4_07` and
`P4_09` add the two sets of sums.

## Inputs

- **Field data** (not shared): one table per survey year,
  `master_raw_<year>.csv`, with the columns described in `data/README.md`.
- **Marine Ecoregions of the World** (Spalding et al., 2007), NoCoast
  version, from UNEP-WCMC.
- **GEE assets** (private to the corresponding author's account): the
  maximum presence probability mask, the Allen Coral Atlas bathymetry and
  geomorphic layers clipped to Indonesia, the presence training points, the ROI
  polygons and the ecoregion polygons. Google Satellite Embeddings
  (`GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL`) and Sentinel-2 are read from the
  public Earth Engine catalog.

The code documents the full workflow and runs with an equivalent set of inputs
under your own GEE project; it does not run as-is without them.

## Data availability

The field data are owned by several Indonesian institutions and were collected
through national survey programmes. Under Indonesian data-sharing policy they
cannot be redistributed publicly. They are available from the corresponding
author on reasonable request and with permission from the data custodians.

The seagrass extent and the annual AGC and cover maps for 2017–2024 are
available at Zenodo (DOI to be added). Format of the files written by
`P4_G10_export_zenodo.js`:

| File | Bands | Type | Values | NoData |
|---|---|---|---|---|
| `P4_seagrass_extent*.tif` | `extent` | UInt8 | 1 = core, 2 = add-on | 0 |
| `P4_seagrass_AGC_cover_<year>*.tif` | `AGC_x10`, `cover_x10` | UInt16 | AGC (gC m⁻²) × 10, cover (%) × 10 | 65535 |

Grid: EPSG:4326, native map grid (about 10 m). Large products are split into
tiles by GEE (suffix `-0000000000-0000000000`).

## Citation

To be added on publication.

## License

MIT License (see [LICENSE](LICENSE)). The field data used to produce the
results are not covered by this license.

## Contact

Muhammad Hafizt, School of the Environment, The University of Queensland,
Australia, and Research Centre for Oceanography, National Research and
Innovation Agency (BRIN), Indonesia. m.hafizt@uq.edu.au
