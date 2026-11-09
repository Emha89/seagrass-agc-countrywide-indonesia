# Input format

The field data are not shared (see Data availability in the main README).
This folder documents the input format only.

## `master_raw_<year>.csv` (one file per survey year, 2017–2023)

Read by `read_master()` in `R/P4_00_config_functions.R` from `CFG$result_dir`.
`master_raw_template.csv` holds the column headers.

| Column | Content |
|---|---|
| `gee_id` | unique sample ID (text) |
| `xcoord`, `ycoord` | longitude and latitude (decimal degrees, WGS 84) |
| `PA` | seagrass presence (1) or absence (0) |
| `tSPC` | total seagrass percent cover (%); values above `CFG$cover_cap` are set to it |
| `AGC_pred` | field above-ground carbon (gC m⁻²), point estimate |
| `AGC_low`, `AGC_up` | lower and upper limits of the 95% prediction interval of `AGC_pred` |
| `GSE_A00` … `GSE_A63` | Google Satellite Embeddings of the survey year at the sample pixel |
| `depth`, `distToLand` | optional, diagnostics only (not used as predictors) |

A sample enters the AGC model when `AGC_pred > 0` and its interval is valid, and
the cover model when `tSPC > 0` (`CFG$cover_train_set = "present"`). The value
`-9999` in a predictor column is read as missing. When `CFG$gse_path` is set,
the GSE values are replaced by those exported by `gee/P4_G00_extract_gse.js`.

## Other inputs

- `MEOW/`: Marine Ecoregions of the World, NoCoast version (UNEP-WCMC), path in
  `CFG$meow_path`.
- Field label file with `gee_id` and a location column (`TAB$location_path` in
  `R/P4_09_tables.R`): used only for the site names of the ROI table.
