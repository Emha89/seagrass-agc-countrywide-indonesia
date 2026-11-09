# =============================================================================
# P4_06_gee_prepare.R
# Paper 4 revised design: the two files to upload to GEE before P4_G01-G04.
#
# Reads : P4_samples.csv (P4_01), MEOW shapefile (CFG$meow_path)
# Writes: <out_dir>/gee_upload/P4_gee_training.csv
#             one row per sample (AGC and/or cover), the 64 GSE values used in R,
#             rid = row id used to join every GEE export back to the samples,
#             missing responses = -9999.
#             Upload: Assets > New > CSV file; X column = xcoord, Y column = ycoord
#             -> chapter_3/output/P4v2/P4_gee_training
#         <out_dir>/gee_upload/P4_meow_indo.shp (+ .shx .dbf .prj .cpg)
#             MEOW ecoregions (TYPE == "MEOW") cropped to Indonesia, integer eco_id.
#             Upload: Assets > New > Shape files (select all five files)
#             -> chapter_3/output/P4v2/P4_meow_indo
# All new GEE assets (inputs and outputs) live in chapter_3/output/P4v2 with a
# P4_ prefix, so no existing asset is overwritten.
#         P4_meow_lookup.csv   eco_id <-> ECOREGION (P4_07 uses it)
# =============================================================================

# setwd("path/to/seagrass-agc-countrywide-indonesia/R")   # folder that holds these scripts
source("P4_00_config_functions.R")

NA_CODE <- -9999
up_dir  <- out_file("gee_upload")
dir.create(up_dir, showWarnings = FALSE, recursive = TRUE)

samples <- read_samples()
if (any(samples$tSPC[samples$is_cov] > CFG$cover_cap, na.rm = TRUE))
  stop(sprintf("P4_samples.csv still has tSPC > %g: rerun P4_01 (cover cap) first.", CFG$cover_cap))

# =============================================================================
# STEP 1: MEOW ecoregions with an integer id (same layer and crop as join_ecoregion)
# =============================================================================
sf::sf_use_s2(FALSE)
meow <- sf::st_read(CFG$meow_path, quiet = TRUE)
meow <- meow[meow$TYPE == "MEOW", "ECOREGION"]
meow <- sf::st_transform(meow, 4326)
meow <- suppressWarnings(suppressMessages(
  sf::st_crop(meow, xmin = 94, ymin = -12, xmax = 142, ymax = 7)))
meow <- suppressWarnings(sf::st_collection_extract(sf::st_make_valid(meow), "POLYGON"))
meow <- meow |>
  dplyr::group_by(ECOREGION) |>
  dplyr::summarise(.groups = "drop") |>
  sf::st_cast("MULTIPOLYGON") |>
  dplyr::arrange(ECOREGION) |>
  dplyr::mutate(eco_id = dplyr::row_number(), .before = 1)

lookup <- sf::st_drop_geometry(meow)
print(as.data.frame(lookup), row.names = FALSE)
readr::write_csv(lookup, out_file("P4_meow_lookup.csv"))

shp <- file.path(up_dir, "P4_meow_indo.shp")
sf::st_write(meow, shp, delete_dsn = file.exists(shp), quiet = TRUE)
cat(sprintf("Saved: %s (%d ecoregions)\n", shp, nrow(meow)))

missing_eco <- setdiff(unique(samples$ecoregion), lookup$ECOREGION)
if (length(missing_eco) > 0)
  stop("Sample ecoregions not in the MEOW lookup: ", paste(missing_eco, collapse = ", "))

# =============================================================================
# STEP 2: Training table for GEE (GSE values included, so GEE and R fit on the
# same inputs)
# =============================================================================
gee_tbl <- samples |>
  dplyr::mutate(rid = dplyr::row_number(),
                eco_id = lookup$eco_id[match(ecoregion, lookup$ECOREGION)]) |>
  dplyr::transmute(rid, gee_id, year, xcoord, ycoord,
                   is_agc = as.integer(is_agc), is_cov = as.integer(is_cov),
                   AGC_pred = ifelse(is_agc, AGC_pred, NA_CODE),
                   tSPC     = ifelse(is_cov, tSPC, NA_CODE),
                   eco_id, g_cluster,
                   dplyr::across(dplyr::all_of(CFG$gse_vars)))

# Checks: no missing values may reach GEE (a text "NA" turns a column into strings)
if (any(as.matrix(gee_tbl[, CFG$gse_vars]) <= NA_CODE))
  stop("GSE fill code (-9999) in the training table: rerun P4_01 first.")
na_cols <- names(gee_tbl)[colSums(is.na(gee_tbl)) > 0]
if (length(na_cols) > 0) stop("Missing values in: ", paste(na_cols, collapse = ", "))
if (anyDuplicated(gee_tbl$gee_id) > 0)
  cat("NOTE: gee_id is not unique; GEE exports are joined on rid instead.\n")

chk <- gee_tbl |>
  dplyr::summarise(rows = dplyr::n(), agc = sum(is_agc), cover = sum(is_cov),
                   agc_min = min(AGC_pred[is_agc == 1]), agc_max = max(AGC_pred[is_agc == 1]),
                   cover_min = min(tSPC[is_cov == 1]), cover_max = max(tSPC[is_cov == 1]))
print(as.data.frame(chk), row.names = FALSE)
print(gee_tbl |> dplyr::count(year, name = "rows") |> as.data.frame(), row.names = FALSE)

csv <- file.path(up_dir, "P4_gee_training.csv")
readr::write_csv(gee_tbl, csv, na = "")
cat(sprintf("Saved: %s (%d rows x %d columns, %.1f MB)\n",
            csv, nrow(gee_tbl), ncol(gee_tbl), file.size(csv) / 1e6))

# Fingerprint: P4_G01 prints the same numbers from the GEE asset. If they differ,
# the asset is an older upload: delete it and upload this file.
fp <- gee_tbl |>
  dplyr::summarise(rows = dplyr::n(), agc = sum(is_agc), cover = sum(is_cov),
                   max_tSPC = max(tSPC[is_cov == 1]), sum_AGC = round(sum(AGC_pred[is_agc == 1]), 1))
cat("\nFINGERPRINT (compare with the P4_G01 console):\n")
print(as.data.frame(fp), row.names = FALSE)

cat("\nUpload to GEE (create the asset folder chapter_3/output/P4v2 first):\n",
    " 1. P4_gee_training.csv -> chapter_3/output/P4v2/P4_gee_training (X = xcoord, Y = ycoord)\n",
    " 2. P4_meow_indo.shp/.shx/.dbf/.prj/.cpg -> chapter_3/output/P4v2/P4_meow_indo\n",
    " 3. Check RF_COV in P4_G01 against P4_best_params.csv, then run P4_G01_train_apply.js\n")

cat("\nP4_06_gee_prepare.R COMPLETE\n")
