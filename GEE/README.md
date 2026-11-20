# Paper 4: Google Earth Engine scripts (design v3)

These scripts map seagrass above-ground carbon (AGC) in Indonesia from Google Satellite Embeddings (GSE, 64 bands only), for 2017–2024. They also compute the national and ecoregion sums and the coverage and temporal checks, and extract the Sentinel-2 predictors of the benchmark model.

Every R step that comes before or after GEE is in the `R/` folder (`P4_00` to `P4_07`).

## What changed in v3

- `P4_G00`–`P4_G09` are the P4v2 scripts. Their assets and Drive tables are reused as they are, because the GSE model and the samples did not change. They are not run again.
- The cover layers of `P4_G01`/`P4_G06`, the cover-class split of `P4_G02` and the ROI tables are no longer used by R (Tier 2 is now based on ecoregion strata, and the ROIs were dropped).
- New: `P4_G11_extract_s2.js` extracts the Sentinel-2 percentile predictors at every AGC sample. It writes to the Drive folder `P4v3_GEE`.

## Setup

Set `ASSET_ROOT` in every script to your own asset folder (placeholder: `projects/YOUR-GEE-PROJECT/assets`).

## Scripts in run order

| Script | Role | Main outputs |
|---|---|---|
| `P4_G00_extract_gse.js` | Re-extracts the 64 GSE values at every training sample from the current dataset (version 1.1, regenerated 2017 layer), so the training values and the maps use the same version | `P4_gee_gse_extract.csv` (read by `P4_01` through `gse_path`) |
| `P4_G01_train_apply.js` | Trains the GSE-only AGC and cover random forests (settings from `P4_02`) and maps 2017–2024 inside the **core** extent | `P4_extent_mask`, `P4_AGC_COV_YYYY`, `P4_gee_points_fitted.csv` |
| `P4_G02_totals.js` | Area and AGC sums by year × ecoregion × cover class, and per ROI (core extent) | `P4_gee_totals_YYYY.csv`, `P4_gee_roi_totals.csv`, `P4_roi_polygons.geojson` |
| `P4_G03_extent_coverage.js` | Samples inside the extent, PA probability at the PA points, and distance to the nearest AGC sample (core extent) | `P4_gee_samples_extent.csv`, `P4_gee_pa_points.csv`, `P4_gee_distance_area.csv` |
| `P4_G04_temporal.js` | Per-pixel AGC CV 2017–2024 and GSE similarity between years (core extent) | `P4_AGC_CV_2017_2024`, `P4_gee_cv_summary.csv`, `P4_gee_gse_sim_*.csv` |
| `P4_G05_extent_sensitivity.js` | Extent area under three geomorphic masks, and a check of the old 100 m area method | `P4_gee_extent_sensitivity_30m.csv` |
| `P4_G06_addon_maps.js` | **Add-on.** Maps AGC and cover for the shallow pixels that `P4_G01` dropped by mistake (see below) | `P4_extent_addon_mask`, `P4_AGC_COV_ADDON_YYYY` |
| `P4_G07_addon_tables.js` | **Add-on.** Computes the same sums as `G02`, `G03` (distance) and `G04` (CV) for the add-on area | `*_addon*.csv`, `P4_AGC_CV_ADDON_2017_2024` |
| `P4_G08_figure_grid.js` | Sums the 10 m pixels of the final extent (core + add-on) into ~500 m cells for the national map figures: area, AGC per year, area per distance class, AGC CV. Drive tables only | `P4_grid500m_agcmean.csv`, `P4_grid500m_agc_YYYY.csv`, `P4_grid500m_dist.csv`, `P4_grid500m_cv.csv` |
| `P4_G08b_zoom_windows.js` | 10 m zoom windows (0.06°) at Rote, Derawan, Kepulauan Seribu and Selayar: annual and mean AGC, mean cover and cover class, distance class, core/add-on part, which extent mask removes each pixel, and a Sentinel-2 true-colour background (Cloud Score+ masked median 2020–2024) | `P4_zoom_<site>_layers.tif`, `P4_zoom_<site>_s2rgb.tif`, `P4_zoom_<site>_agc_annual.tif`, `P4_zoom_<site>_maskdiag.tif` |
| `P4_G09_pa_cv.js` | Repeats the 5-fold cross-validation of the thesis presence (PA) model (training points 2018–2023, 64 GSE bands + depth, 500 trees, 10 variables per split, min leaf 3, bag fraction 0.8, seed 42, same folds) at thresholds 0.5, 0.6 and 0.7. Reads thesis assets only | `P4_gee_pa_cv.csv`, `P4_gee_pa_cv_counts.csv` |
| `P4_G10_export_zenodo.js` | **v3 (20 Nov).** GeoTIFF export for the Zenodo deposit: the final extent (1 = core, 2 = add-on), the annual AGC maps 2017–2024 and their mean, as AGC × 10 (UInt16). The cover band is no longer exported. No field data and no distance layer leave Earth Engine | Drive folder `P4v3_ZENODO` (zipped by `R/P4_11_zenodo_package.R`) |
| `P4_G11_extract_s2.js` | **v3.** Sentinel-2 benchmark predictors at every AGC sample from L2A surface reflectance (`COPERNICUS/S2_SR_HARMONIZED`). The annual mosaic is built like the GSE annual layer: calendar year of the survey, every scene with a Cloud Score+ pixel mask (cs_cdf ≥ 0.60), one observation per satellite pass, bands B2, B3, B4, B8 and B11, percentiles 10, 25, 50, 75 and 90 (25 predictors). The console prints the L2A and L1C scenes per survey year (Earth Engine L2A is not yet global in 2017–2018) | `P4_gee_s2_extract_YYYY.csv` (read by `P4_01` from `CFG$s2_dir`) |

## Core + add-on: why there are two sets of maps

The extent is defined as:
- shallow water (ACA depth 0–500),
- with maximum presence probability ≥ 0.7,
- minus ACA geomorphic class 16, which is removed to limit false seagrass in turbid water.

`P4_G01` reused the thesis mask `geomorphic.neq(16)` together with `updateMask`. In GEE, `neq(16)` on a pixel without ACA geomorphic data returns no data, so `updateMask` drops that pixel. As a result, `P4_G01` also removed every shallow pixel outside the ACA geomorphic layer, about 34% of the intended extent. This included all of Derawan and most of the Malacca Strait and Papua.

`P4_G06` maps only those missing pixels. It uses the same models, settings and grid, and does not change any `P4_G01` asset. `P4_G07` sums them.

**Final product:**
- **Extent:** `P4_extent_mask` OR `P4_extent_addon_mask`
- **Maps:** `ee.ImageCollection([P4_AGC_COV_YYYY, P4_AGC_COV_ADDON_YYYY]).mosaic()`
- **Tables:** `P4_07_gee_results.R` adds the core and add-on sums. All tables are sums, so core + add-on = total.

## Conventions

- **Asset locations:** all new assets are in `chapter_3/output/P4v2` with the prefix `P4_`. Assets from the thesis version are only read, never overwritten.
- **Drive exports:** all tables go to the Drive folder `P4v2_GEE` (`P4_G11`: `P4v3_GEE`). Copy them to `data/gee` of the R project.
- **Grid:** maps are exported on the grid of the MaxProb mask (EPSG:4326, 10 m).
- **Areas:** computed from `ee.Image.pixelArea()`.
- **Training points:** rebuilt from `xcoord`/`ycoord`, because the CSV upload leaves the geometry empty.
