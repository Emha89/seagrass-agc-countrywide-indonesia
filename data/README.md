# Input data

The field data are not shared (see Data availability in the main README).
This folder documents the inputs and their format only. Place them as follows
before running the R scripts from the repository root.

| Folder | Files | Source |
|---|---|---|
| `data/master_raw/` | `master_raw_2017.csv` to `master_raw_2023.csv` | field data, one file per survey year (not shared) |
| `data/MEOW/` | `WCMC-036-MEOW-PPOW-2007-2012-NoCoast.shp` and its sidecar files | UNEP-WCMC Marine Ecoregions of the World |
| `data/gee/` | CSV exports of the GEE scripts | `GEE/P4_G00`–`P4_G09` and `P4_G11` (Google Drive) |
| `data/cover_agb/` | input, training table and output of the cover-AGB model | only for the optional diagnostics `P4_01a` and `P4_01b` (not shared) |
| `data/zenodo_maps/` | GeoTIFFs of `GEE/P4_G10_export_zenodo.js` | read by `R/P4_11_zenodo_package.R` |

R outputs are written to `output/P4v3/` (created by `P4_00`).

## `master_raw_<year>.csv`

Read by `read_master()` in `R/P4_00_config_functions.R` from `CFG$result_dir`.
`master_raw_template.csv` holds the column headers.

| Column | Content |
|---|---|
| `gee_id` | unique sample ID (text) |
| `xcoord`, `ycoord` | longitude and latitude (decimal degrees, WGS 84) |
| `PA` | seagrass presence (1) or absence (0) |
| `tSPC` | total seagrass percent cover (%); values above `CFG$cover_cap` are set to it |
| `AGC_pred` | field above-ground carbon (gC m⁻²) predicted from species cover, point estimate |
| `AGC_low`, `AGC_up` | lower and upper limits of the 95% prediction interval of `AGC_pred` |
| `GSE_A00` … `GSE_A63` | Google Satellite Embeddings of the survey year at the sample pixel |
| `depth`, `distToLand` | optional, diagnostics only (not used as predictors) |

A sample enters the AGC model when `AGC_pred > 0` and its interval is valid.
The value `-9999` in a predictor column is read as missing. When
`CFG$gse_path` is set, the GSE values are replaced by those exported by
`GEE/P4_G00_extract_gse.js`.
