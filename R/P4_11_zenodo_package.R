# =============================================================================
# P4_11_zenodo_package.R
# Paper 4, design v3 (20 Nov 2026): builds the open data deposit for Zenodo
# from processed results only. Field data, sample coordinates, sample-level
# predictions and 10 m distance-to-sample layers are never copied: field data
# are withheld under Indonesian data-sharing policy.
#
# Run from the project root, after P4_09 (tables) and after the GeoTIFFs of
# GEE/P4_G10_export_zenodo.js (Drive folder P4v3_ZENODO) are downloaded into
# ZEN$maps_dir.
#
# Reads : <maps_dir>/P4_seagrass_*.tif        P4_G10 (a product may be split into tiles)
#         <out_dir>/P4_tables_v3/P4_T*_*.csv, P4_S*_*.csv and P4_tables_all.md (P4_09)
#         <gee_dir>/P4_grid500m_dist.csv     P4_G08 (extent area per distance class,
#                                            ~500 m cells), summed here to 0.1 degree
# Writes: <out_dir>/P4_zenodo_v1/            the files to upload, all in one folder
#                                            (Zenodo keeps no folders):
#           P4_seagrass_*.zip   one zip per map product holding all its tiles (a
#                               product that is a single file is copied as .tif)
#           P4_T*.csv, P4_S*.csv the tables of the paper and the Supplement
#           P4_tables_all.md    the same tables with their titles and notes
#           P4_distance_class_0.1deg.csv  extent by distance class per 0.1 degree
#                               cell (the data of Fig. 6a)
#           README.md           description of every file
#           MANIFEST.csv        size and MD5 of every other file
#         <out_dir>/P4_11_checks.txt         privacy and completeness checks (not uploaded)
#
# The script stops before writing anything if a table looks like sample-level
# data (coordinate or sample ID columns, or more rows than ZEN$max_rows).
# =============================================================================

source("R/P4_00_config_functions.R")
cat("P4_11_zenodo_package.R | version 2026-11-21 (v3)\n")

# -----------------------------------------------------------------------------
# CONFIG -- edit here only
# -----------------------------------------------------------------------------
ZEN <- list(
  version   = "1.0",
  # deposit folder (about 5 GB of zips): set a folder on a large local disk if the
  # project folder is synced to the cloud, for example "D:/P4_zenodo_v1"
  dir       = out_file("P4_zenodo_v1"),
  checks    = out_file("P4_11_checks.txt"),
  # GeoTIFFs of P4_G10 (Drive folder P4v3_ZENODO, all download zips extracted); the
  # folder may be outside the project, for example on an external disk (use / in paths)
  maps_dir  = "data/zenodo_maps",
  maps_recursive = TRUE,                           # also search subfolders (one per extracted zip)
  check_tiles = TRUE,                              # every product must have the same tiles (see 1a);
                                                   # FALSE only after checking the tiles in Drive
  tab_dir   = out_file("P4_tables_v3"),            # P4_09
  grid_dist = file.path(if (!is.null(CFG$gee_dir)) CFG$gee_dir else out_file("gee"),
                        "P4_grid500m_dist.csv"),   # P4_G08
  cell_deg  = 0.1,                                 # cell of the distance table (as Fig. 6a)
  years     = 2017:2024,
  zip_maps  = TRUE,                                # FALSE = tiles in one folder per product
  max_files = 100,                                 # Zenodo limit of files per record
  max_rows  = 5000,                                # longer tables are treated as sample-level
  # columns that identify samples or their locations; a table holding one is refused
  forbidden_cols = "^(x|y|lon|lat|long|longitude|latitude|xcoord|ycoord|x_coord|y_coord|rid|gee_id|sample_id|block_id|composition_id|quadrat|quadrat_id|transect|transect_id|site_id|cluster_id|meadow_id|geometry)$",
  # file names that belong to inputs or sample-level outputs; never copied
  forbidden_files = "master_raw|samples|oof|training|measured|paper1|cover_agb|qgis|gse_extract|s2_extract|points|inventory_units|gee_upload",
  doi       = "10.5281/zenodo.XXXXXXX",            # reserved DOI from the Zenodo form
  code_url  = "https://github.com/Emha89/seagrass-agc-countrywide-indonesia",
  title     = "Seagrass above-ground carbon maps of Indonesia at 10 m, 2017-2024",
  paper     = paste("What pixel-based maps add to area-based estimates of seagrass",
                    "above-ground carbon: a national test across Indonesia"),
  authors   = "Hafizt, M., McMahon, K., Roelfsema, C., Hernawan, U.E., Adi, N.S., Lyons, M., Phinn, S."
)

# Distance classes of P4_G03 (0-10) and the validation design that represents
# each class, as in P4_08 (DIST_DESIGN)
DIST_GROUP <- c("le10m", "10_100m", "10_100m", "100m_1km", "100m_1km", "100m_1km",
                "1_20km", "1_20km", "1_20km", "gt20km", "gt20km")
GROUP_LEVELS <- c("le10m", "10_100m", "100m_1km", "1_20km", "gt20km")
GROUP_DESIGN <- c(le10m = "pixel", `10_100m` = "100 m block", `100m_1km` = "1 km block",
                  `1_20km` = "meadow cluster", gt20km = "ecoregion")

MAP_PRODUCTS <- c("P4_seagrass_extent", sprintf("P4_seagrass_AGC_%d", ZEN$years),
                  sprintf("P4_seagrass_AGC_mean_%d_%d", min(ZEN$years), max(ZEN$years)))

log_lines <- character(0)
note <- function(...) { s <- sprintf(...); cat(s, "\n"); log_lines <<- c(log_lines, s) }
problems <- character(0)
fail <- function(...) { s <- sprintf(...); cat("  PROBLEM:", s, "\n"); problems <<- c(problems, s) }

# -----------------------------------------------------------------------------
# 1. CHECK THE INPUTS (nothing is written before every check has passed)
# -----------------------------------------------------------------------------
note("== 1. Inputs")

# 1a. map products
tif <- list.files(ZEN$maps_dir, pattern = "^P4_seagrass_.*\\.tif$", full.names = TRUE,
                  recursive = ZEN$maps_recursive)
dup <- unique(basename(tif)[duplicated(basename(tif))])
if (length(dup) > 0) fail("the same tile is in more than one folder (keep one copy): %s",
                          paste(head(dup, 5), collapse = ", "))
product_of <- function(f) sub("(-\\d+-\\d+)?\\.tif$", "", basename(f))
tile_of    <- function(f) { m <- regmatches(basename(f), regexpr("-\\d+-\\d+(?=\\.tif$)", basename(f), perl = TRUE))
                            if (length(m) == 0) "(single file)" else m }
by_product <- split(tif, vapply(tif, product_of, ""))
for (p in MAP_PRODUCTS) {
  n <- length(by_product[[p]])
  if (n == 0) fail("map product missing in %s: %s (run P4_G10 and download P4v3_ZENODO)", ZEN$maps_dir, p)
  else note("map %-36s %d file(s), %.1f MB", p, n, sum(file.size(by_product[[p]])) / 1e6)
}
# all products share one mask (the final extent), so GEE writes the same tiles for
# each; a tile missing from one product usually means a download zip not extracted
present <- intersect(MAP_PRODUCTS, names(by_product))
if (ZEN$check_tiles && length(present) > 1) {
  tiles <- lapply(by_product[present], function(f) unique(vapply(f, tile_of, "")))
  all_tiles <- sort(unique(unlist(tiles)))
  note("tiles: %d different tile positions over %d products", length(all_tiles), length(present))
  for (p in present) {
    miss <- setdiff(all_tiles, tiles[[p]])
    if (length(miss) > 0) fail("%s lacks %d tile(s) that other products have: %s (extract every download zip; compare with the P4v3_ZENODO folder in Drive)",
                               p, length(miss), paste(head(miss, 6), collapse = ", "))
  }
}
extra <- setdiff(names(by_product), MAP_PRODUCTS)
if (length(extra) > 0) fail("unexpected map files (not copied, remove or rename): %s", paste(extra, collapse = ", "))
if (length(tif) > 0 && requireNamespace("terra", quietly = TRUE)) {
  for (f in tif) {
    r <- try(terra::rast(f), silent = TRUE)
    if (inherits(r, "try-error")) { fail("cannot read %s", basename(f)); next }
    if (terra::nlyr(r) != 1) fail("%s has %d bands (expected 1)", basename(f), terra::nlyr(r))
    if (!isTRUE(terra::is.lonlat(r))) fail("%s is not in geographic coordinates (EPSG:4326)", basename(f))
  }
}

# 1b. tables of the paper and the Supplement
tab_files <- list.files(ZEN$tab_dir, pattern = "^P4_(T\\d+|S\\d+)_.+\\.csv$", full.names = TRUE)
if (length(tab_files) == 0) fail("no tables in %s (run P4_09)", ZEN$tab_dir)
tab_md <- file.path(ZEN$tab_dir, "P4_tables_all.md")
for (f in tab_files) {
  if (grepl(ZEN$forbidden_files, basename(f), ignore.case = TRUE)) { fail("%s: name of an input file", basename(f)); next }
  d <- readr::read_csv(f, show_col_types = FALSE, guess_max = 100000)
  bad <- grep(ZEN$forbidden_cols, names(d), ignore.case = TRUE, value = TRUE)
  if (length(bad) > 0) fail("%s holds sample or location columns: %s", basename(f), paste(bad, collapse = ", "))
  if (nrow(d) > ZEN$max_rows) fail("%s has %d rows (more than %d, looks sample-level)", basename(f), nrow(d), ZEN$max_rows)
}
note("tables: %d CSV files from %s", length(tab_files), ZEN$tab_dir)

# 1c. distance classes per 0.1 degree cell (the data of Fig. 6a)
have_dist <- file.exists(ZEN$grid_dist)
if (!have_dist) note("distance table skipped, missing: %s (P4_G08)", ZEN$grid_dist)

if (length(problems) > 0) {
  writeLines(c(log_lines, "", "PROBLEMS (nothing written):", problems), ZEN$checks)
  stop(sprintf("%d problem(s); nothing was written. See %s", length(problems), ZEN$checks))
}

# -----------------------------------------------------------------------------
# 2. BUILD THE DEPOSIT FOLDER
# -----------------------------------------------------------------------------
note("\n== 2. Deposit folder %s", ZEN$dir)
if (dir.exists(ZEN$dir)) unlink(ZEN$dir, recursive = TRUE)
dir.create(ZEN$dir, recursive = TRUE)

# 2a. maps: one zip per product with all its tiles
zip_files <- function(zipfile, files) {
  if (requireNamespace("zip", quietly = TRUE)) {
    zip::zip(zipfile, files = files, mode = "cherry-pick")   # tiles without their folders
    return(TRUE)
  }
  if (nzchar(Sys.which("zip")) || nzchar(Sys.getenv("R_ZIPCMD"))) {
    return(utils::zip(zipfile, files, flags = "-j9Xq") == 0)
  }
  FALSE
}
manual_zip <- character(0)
for (p in MAP_PRODUCTS) {
  f <- sort(by_product[[p]])
  if (length(f) == 1) {
    invisible(file.copy(f, file.path(ZEN$dir, paste0(p, ".tif"))))
    next
  }
  zf <- file.path(normalizePath(ZEN$dir), paste0(p, ".zip"))
  if (!(ZEN$zip_maps && zip_files(zf, normalizePath(f)))) {
    dir.create(file.path(ZEN$dir, p))
    invisible(file.copy(f, file.path(ZEN$dir, p)))
    manual_zip <- c(manual_zip, p)
  }
}
if (length(manual_zip) > 0)
  note("No zip tool found (install the R package 'zip'). Compress these folders by hand into <name>.zip, then delete the folders: %s",
       paste(manual_zip, collapse = ", "))

# 2b. tables
invisible(file.copy(tab_files, ZEN$dir))
if (file.exists(tab_md)) invisible(file.copy(tab_md, ZEN$dir))

# 2c. extent by distance class per 0.1 degree cell
if (have_dist) {
  gd <- readr::read_csv(ZEN$grid_dist, show_col_types = FALSE)
  dcols <- sprintf("area_d%02d", 0:10)
  stopifnot(all(c("lon", "lat", dcols) %in% names(gd)))
  cell <- ZEN$cell_deg
  cells <- gd |>
    dplyr::mutate(lon_center = floor(lon / cell) * cell + cell / 2,
                  lat_center = floor(lat / cell) * cell + cell / 2) |>
    dplyr::group_by(lon_center, lat_center) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(dcols), sum), .groups = "drop")
  grp <- sapply(GROUP_LEVELS, function(g) rowSums(as.matrix(cells[, dcols[DIST_GROUP == g], drop = FALSE]))) / 1e6
  if (is.null(dim(grp))) grp <- matrix(grp, nrow = 1, dimnames = list(NULL, GROUP_LEVELS))
  out <- tibble::tibble(lon_center = round(cells$lon_center, 4), lat_center = round(cells$lat_center, 4),
                        extent_km2 = round(rowSums(grp), 4))
  for (g in GROUP_LEVELS) out[[paste0("km2_", g)]] <- round(grp[, g], 4)
  out$main_class <- GROUP_LEVELS[max.col(grp, ties.method = "first")]
  out <- dplyr::filter(out, extent_km2 > 0)
  readr::write_csv(out, file.path(ZEN$dir, "P4_distance_class_0.1deg.csv"))
  note("distance table: %d cells of %.1f degree, %.2f km2 of extent", nrow(out), cell, sum(out$extent_km2))
}

# -----------------------------------------------------------------------------
# 3. README
# -----------------------------------------------------------------------------
map_line <- function(p, desc) {
  f <- if (file.exists(file.path(ZEN$dir, paste0(p, ".tif")))) paste0(p, ".tif") else paste0(p, ".zip")
  sprintf("| `%s` | %s |", f, desc)
}
yr_desc <- "Seagrass AGC of the year, band `AGC_x10` (UInt16): AGC in gC m-2 x 10. NoData 65535."
readme <- c(
  sprintf("# %s", ZEN$title), "",
  sprintf("Version %s. Data deposit of the article \"%s\" (%s).", ZEN$version, ZEN$paper, ZEN$authors), "",
  "## Content", "",
  paste("Annual maps of seagrass above-ground carbon (AGC, gC m-2) for Indonesia from 2017 to 2024 at 10 m,",
        "their mean, and the fixed seagrass extent on which they were mapped. AGC was predicted with a random",
        "forest on the 64 bands of the Google Satellite Embeddings (AlphaEarth Foundations), trained on field",
        "samples of 2017-2023. No field data were collected in 2024. The deposit also holds the tables of the",
        "article and its Supplement as CSV, and the seagrass extent by distance to the nearest field sample",
        "per 0.1 degree cell."), "",
  paste("Map accuracy depends on the distance to field samples. The maps estimated unsampled sites within",
        "sampled meadows well, but far from field samples the predictions approach the national mean. Use the",
        "maps for national and large-area totals and for stratification, and see the article for the full",
        "assessment."), "",
  "## Files", "",
  "| File | Content |", "|---|---|",
  map_line("P4_seagrass_extent",
           "Seagrass extent, band `extent` (UInt8): 1 = core extent, 2 = add-on extent (shallow pixels without Allen Coral Atlas geomorphic data). NoData 0."),
  vapply(sprintf("P4_seagrass_AGC_%d", ZEN$years), map_line, "", desc = yr_desc),
  map_line(sprintf("P4_seagrass_AGC_mean_%d_%d", min(ZEN$years), max(ZEN$years)),
           "Mean of the eight annual AGC maps, band `AGC_x10` (UInt16, AGC x 10). NoData 65535."),
  "| `P4_T*.csv`, `P4_S*.csv` | Tables of the article (T) and the Supplement (S), one file per table, UTF-8 |",
  "| `P4_tables_all.md` | The same tables with their titles and notes |",
  if (have_dist) paste("| `P4_distance_class_0.1deg.csv` | Seagrass extent (km2) per 0.1 degree cell by distance to the",
                       "nearest AGC field sample: `km2_le10m`, `km2_10_100m`, `km2_100m_1km`, `km2_1_20km`, `km2_gt20km`;",
                       "`main_class` holds the largest share; `lon_center`, `lat_center` are cell centres (EPSG:4326) |"),
  "| `MANIFEST.csv` | Size and MD5 checksum of every file |", "",
  "## Map format", "",
  "- Cloud Optimized GeoTIFF, EPSG:4326 (WGS 84), native map grid of about 10 m.",
  "- Divide `AGC_x10` by 10 to obtain AGC in gC m-2. Pixels outside the extent are NoData.",
  "- Pixel area changes with latitude, so compute areas per pixel (for example `terra::cellSize()` in R).",
  "- A large product is split into tiles (suffix `-0000000000-0000000000`). Unzip all tiles of a product and read them as one layer, for example with `terra::vrt()` or `gdalbuildvrt`.",
  "", "```r",
  "library(terra)",
  "unzip(\"P4_seagrass_AGC_2023.zip\", exdir = \"agc2023\")",
  "r   <- vrt(list.files(\"agc2023\", \"\\\\.tif$\", full.names = TRUE))",
  "r   <- crop(r, ext(115.0, 116.0, -9.0, -8.0))      # crop first: the national 10 m grid is large",
  "agc <- r / 10                                       # gC m-2",
  "tot <- global(agc * cellSize(agc, unit = \"m\"), \"sum\", na.rm = TRUE)[1, 1] / 1e6   # t C",
  "```", "",
  "## Not included", "",
  paste("Field data, sample locations and sample-level predictions are not included, under Indonesian",
        "data-sharing policy. They are available from the corresponding author on reasonable request, with",
        "permission from the data custodians."), "",
  "## Code", "", sprintf("R and Google Earth Engine scripts: %s", ZEN$code_url), "",
  "## Licence", "", "Creative Commons Attribution 4.0 International (CC BY 4.0).", "",
  "## Citation", "",
  sprintf("%s (%s). %s (Version %s) [Data set]. Zenodo. https://doi.org/%s",
          ZEN$authors, format(Sys.Date(), "%Y"), ZEN$title, ZEN$version, ZEN$doi)
)
writeLines(readme, file.path(ZEN$dir, "README.md"), useBytes = TRUE)

# -----------------------------------------------------------------------------
# 4. MANIFEST AND FINAL CHECKS
# -----------------------------------------------------------------------------
all_f <- list.files(ZEN$dir, recursive = TRUE, full.names = TRUE)
all_f <- all_f[basename(all_f) != "MANIFEST.csv"]
man <- tibble::tibble(file = sub(paste0("^", ZEN$dir, "/"), "", all_f),
                      bytes = file.size(all_f), md5 = unname(tools::md5sum(all_f)))
readr::write_csv(man, file.path(ZEN$dir, "MANIFEST.csv"))

note("\n== 3. Deposit")
n_files <- length(all_f) + 1
note("files: %d (Zenodo allows %d per record); total %.1f MB", n_files, ZEN$max_files, sum(man$bytes) / 1e6)
if (n_files > ZEN$max_files) note("PROBLEM: too many files; set ZEN$zip_maps = TRUE")
leak <- grep(ZEN$forbidden_files, man$file, ignore.case = TRUE, value = TRUE)
note("files with input-like names: %s", if (length(leak) == 0) "none" else paste(leak, collapse = ", "))
if (grepl("X{7}", ZEN$doi)) note("ZEN$doi is still a placeholder: reserve the DOI on Zenodo, set it here and rerun")
print(as.data.frame(man), row.names = FALSE)
writeLines(log_lines, ZEN$checks)
cat(sprintf("\nUpload every file in %s to Zenodo.\n", ZEN$dir))
