# =============================================================================
# P4_08b_qgis_layers.R
# Paper 4, design v3 (17 Nov 2026): GIS layers for drawing Fig. 2 (study area)
# in QGIS. Nothing is recomputed: the layers are the samples, DBSCAN meadows,
# tested units and ecoregions already used by P4_01 to P4_09.
#
# Run from the project root, after P4_04 (inventory units) and P4_07 (ecoregion
# summary); the extent layers also need the P4_G08 grid export in gee_dir.
#
# Reads : P4_samples.csv (P4_01), P4_inventory_units.csv (P4_04),
#         P4_ecoregion_summary.csv (P4_07), gee_upload/P4_meow_indo.shp (P4_06),
#         optional <gee_dir>/P4_grid500m_agcmean.csv (P4_G08)
# Writes (new files only), in <out_dir>/P4_qgis_v3/:
#   P4_Fig2_qgis_layers_v3.gpkg   one GeoPackage, EPSG:4326, with the layers
#     ecoregions        MEOW ecoregion polygons, with extent area, AGC samples
#                       and Tier 2 stratum (n >= CFG$min_n_eco)
#     ecoregion_lines   the same ecoregions as boundary lines
#     meadows           one polygon per DBSCAN meadow: convex hull of its AGC
#                       samples plus a buffer (P8B$buffer_m), so small meadows
#                       stay visible at national scale; tested = held out as an
#                       inventory unit (n >= CFG$min_n_meadow)
#     sites_tested      100 m blocks tested as sites (n >= CFG$min_n_site in the
#                       main survey year), as squares
#     samples_agc       AGC samples (6,616), with year, ecoregion, AGC, cover,
#                       meadow (0 = unclustered) and the Sentinel-2 set flag
#     samples_cover     quadrats with seagrass cover but no AGC value
#     unclustered       AGC samples that DBSCAN left outside any meadow
#     extent_cells      seagrass extent per ~500 m cell (points, km2), when the
#                       P4_G08 grid file is present
#   P4_Fig2_extent_500m_v3.tif    the same extent cells as a raster (km2 per cell)
#   P4_Fig2_qgis_layer_list_v3.csv  every layer with its geometry, features,
#                                   content and a suggested style
# Base layers not written here: land (e.g. Natural Earth 10 m land or GADM
# Indonesia) and the EEZ (Marine Regions, World EEZ v12).
# The point layers hold field sample locations: keep the GeoPackage with the
# field data (not for redistribution).
# =============================================================================

source("R/P4_00_config_functions.R")
cat("P4_08b_qgis_layers.R | version 2026-11-17 (v3, Fig 2 layers for QGIS)\n")

for (p in c("sf", "terra")) {
  if (!requireNamespace(p, quietly = TRUE)) stop(sprintf("Package '%s' is needed: install.packages('%s')", p, p))
}

# -----------------------------------------------------------------------------
# CONFIG -- edit here only
# -----------------------------------------------------------------------------
P8B <- list(
  dir      = out_file("P4_qgis_v3"),
  gpkg     = "P4_Fig2_qgis_layers_v3.gpkg",
  raster   = "P4_Fig2_extent_500m_v3.tif",
  list_csv = "P4_Fig2_qgis_layer_list_v3.csv",
  buffer_m = 500,                                            # buffer around the meadow hulls (m)
  crs_m    = "+proj=merc +datum=WGS84 +units=m +no_defs"     # metric CRS for hulls and buffers
)

dir.create(P8B$dir, showWarnings = FALSE, recursive = TRUE)
gpkg <- file.path(P8B$dir, P8B$gpkg)
if (file.exists(gpkg)) file.remove(gpkg)          # rewrites the GeoPackage of this script only
sf::sf_use_s2(FALSE)
LAYERS <- list()
add_layer <- function(x, name, geometry, content, style) {
  sf::st_write(x, gpkg, layer = name, quiet = TRUE, append = file.exists(gpkg))
  LAYERS[[name]] <<- tibble::tibble(file = P8B$gpkg, layer = name, geometry = geometry,
                                    features = nrow(x), content = content, suggested_style = style)
  cat(sprintf("  %-16s %6d features\n", name, nrow(x)))
}

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
samples <- read_samples()
units_f <- out_file("P4_inventory_units.csv")
eco_f   <- out_file("P4_ecoregion_summary.csv")
meow_f  <- out_file("gee_upload/P4_meow_indo.shp")
units <- if (file.exists(units_f)) readr::read_csv(units_f, show_col_types = FALSE) else {
  cat("  -> missing: P4_inventory_units.csv (run P4_04); tested flags and sites skipped\n"); NULL }
eco_sum <- if (file.exists(eco_f)) readr::read_csv(eco_f, show_col_types = FALSE) else {
  cat("  -> missing: P4_ecoregion_summary.csv (run P4_07); ecoregion attributes skipped\n"); NULL }

pts <- sf::st_as_sf(samples, coords = c("xcoord", "ycoord"), crs = 4326, remove = FALSE)
pts$meadow <- ifelse(is.na(pts$dbscan_cluster), 0L, as.integer(pts$dbscan_cluster))
s2_flag <- if ("has_s2" %in% names(pts)) pts$has_s2 %in% TRUE else rep(NA, nrow(pts))
pts$s2_set <- s2_flag

agc <- pts[pts$is_agc, c("gee_id", "year", "ecoregion", "AGC_pred", "AGC_low", "AGC_up", "tSPC",
                         "meadow", "g_cluster", "g_b100", "s2_set", "xcoord", "ycoord")]
cov_only <- pts[pts$is_cov & !pts$is_agc, c("gee_id", "year", "ecoregion", "tSPC", "xcoord", "ycoord")]
noise <- agc[agc$meadow == 0, ]

# -----------------------------------------------------------------------------
# Ecoregions (polygons and lines)
# -----------------------------------------------------------------------------
cat("Writing", gpkg, "\n")
if (file.exists(meow_f)) {
  eco <- sf::st_read(meow_f, quiet = TRUE)
  names(eco)[names(eco) == "ECOREGION"] <- "ecoregion"
  if (!is.null(eco_sum)) {
    eco <- dplyr::left_join(eco, eco_sum |>
                              dplyr::transmute(eco_id, extent_km2 = round(area_km2, 2), n_agc,
                                               tier2_stratum = t2_stratum %in% TRUE,
                                               field_mean_gC_m2 = round(field_mean_gC_m2, 1)),
                            by = "eco_id")
    eco <- eco[!is.na(eco$extent_km2), ]           # the 12 ecoregions that hold seagrass extent
  }
  add_layer(eco, "ecoregions", "polygon",
            "MEOW ecoregions holding seagrass extent: extent area, AGC samples, Tier 2 stratum (n >= 30)",
            "no fill or a light fill by tier2_stratum; label with ecoregion")
  lines <- suppressWarnings(sf::st_cast(sf::st_boundary(sf::st_make_valid(eco)), "MULTILINESTRING"))
  add_layer(lines[, intersect(c("eco_id", "ecoregion"), names(lines))], "ecoregion_lines", "line",
            "boundaries of the same ecoregions", "thin dark grey line (0.3 mm)")
} else cat("  -> missing: gee_upload/P4_meow_indo.shp (run P4_06); ecoregion layers skipped\n")

# -----------------------------------------------------------------------------
# Meadows: convex hull of the AGC samples of each DBSCAN cluster, buffered
# -----------------------------------------------------------------------------
agc_m <- sf::st_transform(agc[agc$meadow > 0, ], P8B$crs_m)
cl_ids <- sort(unique(agc_m$meadow))
hulls <- lapply(cl_ids, function(k) {
  g <- sf::st_union(sf::st_geometry(agc_m[agc_m$meadow == k, ]))
  sf::st_buffer(sf::st_convex_hull(g), P8B$buffer_m)
})
g_of <- vapply(cl_ids, function(k) agc_m$g_cluster[agc_m$meadow == k][1], "")
tested_all <- if (!is.null(units)) unique(units$unit[units$level == "meadow" & units$sample_set == "all"]) else character(0)
tested_s2  <- if (!is.null(units)) unique(units$unit[units$level == "meadow" & units$sample_set == "s2"]) else character(0)
meadows <- sf::st_sf(
  meadow     = cl_ids,
  unit       = g_of,
  ecoregion  = vapply(cl_ids, function(k) names(which.max(table(agc_m$ecoregion[agc_m$meadow == k]))), ""),
  n_agc      = vapply(cl_ids, function(k) sum(agc_m$meadow == k), 0L),
  years      = vapply(cl_ids, function(k) paste(sort(unique(agc_m$year[agc_m$meadow == k])), collapse = ", "), ""),
  tested     = g_of %in% tested_all,
  tested_s2  = g_of %in% tested_s2,
  geometry   = do.call(c, hulls), crs = P8B$crs_m) |>
  sf::st_transform(4326)
add_layer(meadows, "meadows", "polygon",
          sprintf("DBSCAN meadows (5 km radius, >= 5 samples): hull of the AGC samples + %d m buffer; tested = held out as a meadow unit", P8B$buffer_m),
          "outline only, solid for tested = TRUE and dashed otherwise")

# -----------------------------------------------------------------------------
# Sites tested: 100 m blocks (equirectangular metric grid of P4_00 block_id)
# -----------------------------------------------------------------------------
if (!is.null(units)) {
  su <- units |> dplyr::filter(level == "site", sample_set == "all") |>
    dplyr::distinct(unit, ecoregion, year, n)
  blk <- sub("\\|.*$", "", su$unit)                       # b100_<ix>_<iy>
  ij  <- do.call(rbind, strsplit(sub("^b100_", "", blk), "_"))
  size <- CFG$block_sizes[["b100"]]
  to_deg <- function(m) m / (EARTH_R * pi / 180)
  squares <- lapply(seq_len(nrow(ij)), function(i) {
    x0 <- to_deg(as.numeric(ij[i, 1]) * size); x1 <- to_deg((as.numeric(ij[i, 1]) + 1) * size)
    y0 <- to_deg(as.numeric(ij[i, 2]) * size); y1 <- to_deg((as.numeric(ij[i, 2]) + 1) * size)
    sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))
  })
  sites <- sf::st_sf(site = su$unit, ecoregion = su$ecoregion, year = su$year, n_agc = su$n,
                     geometry = sf::st_sfc(squares, crs = 4326))
  add_layer(sites, "sites_tested", "polygon",
            "100 m blocks tested as sites (>= 10 AGC samples in their main survey year)",
            "small filled square or its centroid as a marker at national scale")
}

# -----------------------------------------------------------------------------
# Samples
# -----------------------------------------------------------------------------
add_layer(agc, "samples_agc", "point",
          "AGC samples: year, ecoregion, AGC and its 95% interval, cover, meadow (0 = unclustered), Sentinel-2 set",
          "small dot coloured by year; or by s2_set for the Sentinel-2 comparison")
add_layer(cov_only, "samples_cover", "point",
          "quadrats with seagrass cover but no AGC value (presences used for the extent model)",
          "small open circle, grey")
add_layer(noise, "unclustered", "point",
          "AGC samples outside any DBSCAN meadow (each its own group in validation)",
          "small cross")

# -----------------------------------------------------------------------------
# Extent per ~500 m cell (P4_G08), as points and as a raster
# -----------------------------------------------------------------------------
grid_f <- file.path(CFG$gee_dir, "P4_grid500m_agcmean.csv")
if (file.exists(grid_f)) {
  ga <- readr::read_csv(grid_f, show_col_types = FALSE) |> dplyr::mutate(area_km2 = area_m2 / 1e6)
  add_layer(sf::st_as_sf(ga[, c("lon", "lat", "area_km2")], coords = c("lon", "lat"), crs = 4326),
            "extent_cells", "point", "seagrass extent per ~500 m cell (km2)",
            "graduated marker size, or use the raster")
  res <- stats::median(diff(sort(unique(round(ga$lon, 7)))))
  tmpl <- terra::rast(xmin = min(ga$lon) - res / 2, xmax = max(ga$lon) + res / 2,
                      ymin = min(ga$lat) - res / 2, ymax = max(ga$lat) + res / 2,
                      resolution = res, crs = "EPSG:4326")
  r <- terra::rasterize(as.matrix(ga[, c("lon", "lat")]), tmpl, values = ga$area_km2, fun = "sum")
  names(r) <- "extent_km2"
  terra::writeRaster(r, file.path(P8B$dir, P8B$raster), overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  LAYERS[["extent_raster"]] <- tibble::tibble(file = P8B$raster, layer = "extent_km2", geometry = "raster",
                                              features = sum(!is.na(terra::values(r))),
                                              content = "seagrass extent per ~500 m cell (km2)",
                                              suggested_style = "single-band pseudocolour, transparent where empty")
  cat(sprintf("Extent cells: %d (%.2f km2), raster cell %.5f deg\n", nrow(ga), sum(ga$area_km2), res))
} else cat("  -> missing: P4_grid500m_agcmean.csv (P4_G08); extent layers skipped\n")

lst <- dplyr::bind_rows(LAYERS)
readr::write_csv(lst, file.path(P8B$dir, P8B$list_csv))
print_tbl(lst[, c("layer", "geometry", "features")])
cat(sprintf("AGC samples %d | cover-only %d | unclustered %d | meadows %d (tested %d)\n",
            nrow(agc), nrow(cov_only), nrow(noise), nrow(meadows), sum(meadows$tested)))
cat("Written:", P8B$dir, "\n")
cat("\nP4_08b_qgis_layers.R COMPLETE\n")
