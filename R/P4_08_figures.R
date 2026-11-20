# =============================================================================
# P4_08_figures.R
# Paper 4, design v3 (15 Nov 2026): all figures of the manuscript and the
# Supplement, built from the R results (P4_04, P4_05, P4_07) and the GEE
# figure exports (P4_G08 grid, P4_G08b zoom windows; the P4v2 exports are reused).
#
# Main figures (numbers in FIGNUM below; renumber there when the text is written)
#   Fig 1  Workflow diagram (drawn separately, not by this script)
#   Fig 2  Study area: ecoregions, extent, AGC samples; extent share vs sample share
#   Fig 3  Mean AGC 2017-2024: national grid + four 10 m zoom windows
#   Fig 4  Totals by method: national 2017-2024 (95 % CI), survey composition, ecoregions (aim 1)
#   Fig 5  Inventory-level validation: unit mean predicted vs field for sites, meadows and
#          ecoregions held out, and MEC of every method in both sample sets (aim 2)
#   Fig 6  Where new samples are needed: distance to AGC samples (map) and the extent of
#          each ecoregion by distance, with Tier 2 strata status (aim 3)
# Supplement (provisional order of first citation)
#   S1 presence probability threshold by year       S2 extent components and mask diagnostic
#   S3 extent under the geomorphic mask options      S4 annual AGC maps
#   S5 temporal behaviour (CV map, GSE similarity, annual AGC at the window samples)
#   S6 sample-level accuracy against distance, all methods, both sample sets
#   S7 mean GSE prediction by observed AGC (shrinkage toward the mean)
#   S8 cumulative extent by distance to AGC samples
# Dropped from P4v2: the 22 ROIs (Fig 7), the cover model (S6, S9), the
# leave-one-ecoregion-out bias panel and the group-mean shrinkage panel (Fig 5 shows
# the held-out meadows and ecoregions).
#
# Reads (out_dir): P4_samples.csv, P4_metrics_summary.csv, P4_inventory_units.csv,
#   P4_inventory_metrics.csv, P4_national_three_methods.csv, P4_ecoregion_summary.csv,
#   P4_tier2_national.csv, P4_distance_coverage_class.csv, P4_distance_coverage_eco.csv,
#   P4_error_by_agc_bin.csv, P4_pa_threshold_by_year.csv,
#   P4_gse_similarity_pairs.csv, P4_extent_sensitivity_ecoregion.csv, P4_pixel_national.csv,
#   gee_upload/P4_meow_indo.shp
# Reads (CFG$gee_dir): P4_grid500m_agcmean.csv, P4_grid500m_dist.csv, P4_grid500m_cv.csv,
#   P4_grid500m_agc_YYYY.csv, P4_zoom_<site>_layers.tif, _s2rgb.tif, _agc_annual.tif, _maskdiag.tif
# Writes (FIG$dir = <out_dir>/P4_figures_v3): every figure (PNG preview + RSE file:
#   maps as LZW TIFF 600 dpi, charts as vector PDF), and the CSVs behind numbers quoted
#   in the text: P4_Fig3_window_stats.csv, P4_FigS5_window_summary.csv,
#   P4_FigS5_window_years.csv, P4_Fig6_ecoregion_priority.csv.
# A missing input only skips the figures that need it. All export lines are active.
#
# Expected warning (harmless): "Raster pixels are placed at uneven horizontal
# intervals" in the zoom windows (columns without seagrass pixels leave gaps).
# =============================================================================

source("R/P4_00_config_functions.R")
cat("P4_08_figures.R | version 2026-11-16 (v3, labels of S1 and S3)\n")

for (p in c("terra", "patchwork", "maps")) {
  if (!requireNamespace(p, quietly = TRUE)) stop(sprintf("Package '%s' is needed: install.packages('%s')", p, p))
}
suppressPackageStartupMessages({ library(terra); library(patchwork) })

# -----------------------------------------------------------------------------
# FIGURE CONFIG -- edit here only
# -----------------------------------------------------------------------------
FIG <- list(
  dir        = out_file("P4_figures_v3"),
  gee_dir    = if (!is.null(CFG$gee_dir)) CFG$gee_dir else out_file("gee"),
  grid_m     = 500,                    # cell size of the P4_G08 grid exports
  nat_cell   = 0.1,                    # degrees, national map cells (Fig 3, 6, S2, S5)
  ext_cell   = 0.25,                   # degrees, extent cells in Fig 2 and the annual maps (S4)
  nat_xlim   = c(94.5, 141.5),
  nat_ylim   = c(-11.5, 6.5),
  years      = 2017:2024,
  agc_lim    = c(15, 40),              # colour limits of the AGC maps (gC m-2); values outside use
                                       # the end colours, legend ends read "<= 15" and ">= 40"
  s2_fact    = 2,                      # aggregate the Sentinel-2 background (2 = 20 m) for speed
  kkp_km2    = 6601.56,                # KKP 2025 national seagrass map (official area)
  turbid_eco = c("Sunda Shelf/Java Sea"),   # S3: Terrestrial Reef Flat removed only here
  trf_label  = "Terrestrial Reef Flat",     # ACA geomorphic class removed from the extent
  s2_set_label = "AGC samples 2019-2023 with a clear Sentinel-2 pixel",
  # zoom windows (P4_G08b): file stem, label shown in the figures
  windows = tibble::tribble(
    ~stem,      ~label,
    "Rote2",    "Rote",
    "Derawan2", "Derawan",
    "Seribu",   "Seribu Islands",
    "Selayar",  "Selayar"),
  base_size = 9
)
# Figure numbers by key: change here when the order of first citation is known
FIGNUM <- c(study = "2", agcmap = "3", totals = "4", inventory = "5", sampling = "6",
            pa = "S1", extent = "S2", kkp = "S3", annual = "S4", temporal = "S5",
            samplelevel = "S6", shrinkage = "S7", cumdist = "S8")

DESIGN_LEVELS <- c("pixel", "b100", "b1000", "cluster", "ecoregion")
DESIGN_LABEL  <- c(pixel = "Pixel (<= 10 m)", b100 = "100 m block (10-100 m)",
                   b1000 = "1 km block (100 m-1 km)", cluster = "Meadow cluster (1-20 km)",
                   ecoregion = "Ecoregion (> 20 km)")
DESIGN_COL    <- c(pixel = "#1a9850", b100 = "#91cf60", b1000 = "#d9ef8b",
                   cluster = "#fc8d59", ecoregion = "#d73027")
DIST_DESIGN   <- c("pixel", "b100", "b100", "b1000", "b1000", "b1000",
                   "cluster", "cluster", "cluster", "ecoregion", "ecoregion")   # classes 0-10
METHOD_COL    <- c("Tier 1" = "#762a83", "Tier 2" = "#e08214", "Local mean" = "#8c510a",
                   "Pixel-based GSE" = "#1b7837", "Pixel-based S2" = "#2166ac",
                   "Pixel-based S2 (median)" = "#67a9cf")
SET_LABEL     <- c(all = "All AGC samples", s2 = FIG$s2_set_label)
LEVEL_LABEL   <- c(site = "Site\n(100 m block)", meadow = "Meadow\n(DBSCAN cluster)", ecoregion = "Ecoregion")
wrap_cap <- function(x, width = 135) paste(strwrap(x, width = width), collapse = "\n")   # captions in mm-wide figures
# AGC maps: viridis (perceptually uniform, readable with colour-vision deficiency)
AGC_PAL       <- c("#440154", "#472D7B", "#3B528B", "#2C728E", "#21908C", "#27AD81", "#5DC863", "#AADC32", "#FDE725")
agc_scale <- function(...) {
  b <- pretty(FIG$agc_lim); b <- b[b >= FIG$agc_lim[1] & b <= FIG$agc_lim[2]]
  lab <- as.character(b)
  lab[b == FIG$agc_lim[1]] <- paste0("<= ", FIG$agc_lim[1])
  lab[b == FIG$agc_lim[2]] <- paste0(">= ", FIG$agc_lim[2])
  scale_fill_gradientn(colours = AGC_PAL, limits = FIG$agc_lim, oob = scales::squish,
                       breaks = b, labels = lab, ...)
}
YEAR_COL      <- setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666"),
                          2017:2024)

theme_p4 <- function() theme_bw(base_size = FIG$base_size) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"),
        legend.key.size = unit(3.5, "mm"))
theme_map <- function() theme_p4() + theme(axis.title = element_blank(), panel.grid = element_blank())

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
gfile <- function(n) file.path(FIG$gee_dir, n)
ofile <- function(n) out_file(n)
rd    <- function(p) readr::read_csv(p, show_col_types = FALSE, guess_max = 100000)
need  <- function(...) {
  f <- c(...); ok <- file.exists(f)
  if (!all(ok)) cat("  -> skipped, missing:", paste(basename(f[!ok]), collapse = ", "), "\n")
  all(ok)
}
fig <- function(key, txt) cat(sprintf("\n==== Figure %s: %s ====\n", FIGNUM[[key]], txt))
dir.create(FIG$dir, showWarnings = FALSE, recursive = TRUE)
fig_out <- function(n) file.path(FIG$dir, n)
# file name from the key: P4_Fig<number>_<key>.<ext>
fig_file <- function(key, ext) fig_out(sprintf("P4_Fig%s_%s.%s", FIGNUM[[key]], key, ext))
save_map <- function(key, p, w, h, dpi_png = 400) {
  ggsave(fig_file(key, "png"), p, width = w, height = h, units = "mm", dpi = dpi_png)
  ggsave(fig_file(key, "tif"), p, width = w, height = h, units = "mm", dpi = 600,
         device = "tiff", compression = "lzw")                                   # RSE submission (map)
}
save_chart <- function(key, p, w, h) {
  ggsave(fig_file(key, "png"), p, width = w, height = h, units = "mm", dpi = 400)
  ggsave(fig_file(key, "pdf"), p, width = w, height = h, units = "mm", device = cairo_pdf)  # RSE (vector)
}
grid_file <- function(what) gfile(sprintf("P4_grid%dm_%s.csv", FIG$grid_m, what))
zoom_file <- function(stem, what) gfile(sprintf("P4_zoom_%s_%s.tif", stem, what))

# Sum grid cells into coarser cells (degrees)
agg_cells <- function(d, cell, cols) {
  d |>
    dplyr::mutate(gx = floor(lon / cell) * cell + cell / 2, gy = floor(lat / cell) * cell + cell / 2) |>
    dplyr::group_by(gx, gy) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(cols), sum), .groups = "drop")
}

world <- ggplot2::map_data("world")
land  <- function() geom_polygon(data = world, aes(long, lat, group = group), inherit.aes = FALSE,
                                 fill = "grey85", colour = "grey55", linewidth = 0.1)
lon_lab <- function(x) ifelse(x == 0, "0°", paste0(abs(x), "°", ifelse(x < 0, "W", "E")))
lat_lab <- function(y) ifelse(y == 0, "0°", paste0(abs(y), "°", ifelse(y < 0, "S", "N")))
nat_axes  <- function() list(scale_x_continuous(labels = lon_lab), scale_y_continuous(labels = lat_lab))
nat_coord <- function() c(list(coord_quickmap(xlim = FIG$nat_xlim, ylim = FIG$nat_ylim, expand = FALSE)), nat_axes())
pct_lab   <- function(p) ifelse(p < 0.05, "< 0.1%", sprintf("%.1f%%", p))
KM2 <- "km²"

scalebar <- function(x0, y0, km, h, col = "black", size = 2.4) {
  dx <- km / (111.32 * cos(y0 * pi / 180))
  list(annotate("rect", xmin = x0, xmax = x0 + dx, ymin = y0, ymax = y0 + h, fill = col),
       annotate("text", x = x0 + dx / 2, y = y0 + 1.8 * h, label = paste(km, "km"), size = size, colour = col,
                vjust = 0))
}

s2_layer <- function(stem) {
  r <- terra::rast(zoom_file(stem, "s2rgb"))
  if (FIG$s2_fact > 1) r <- terra::aggregate(r, FIG$s2_fact, mean, na.rm = TRUE)
  v <- terra::values(r)
  st <- function(x) { q <- stats::quantile(x, c(0.02, 0.98), na.rm = TRUE); pmin(pmax((x - q[1]) / (q[2] - q[1]), 0), 1) }
  col <- grDevices::rgb(st(v[, 1]), st(v[, 2]), st(v[, 3]))
  col[is.na(v[, 1])] <- "#000000"
  e <- as.vector(terra::ext(r))
  annotation_raster(matrix(col, nrow = terra::nrow(r), ncol = terra::ncol(r), byrow = TRUE),
                    xmin = e[1], xmax = e[2], ymin = e[3], ymax = e[4], interpolate = FALSE)
}
rast_df <- function(file, band) {
  r <- terra::rast(file)[[band]]
  d <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE); names(d)[3] <- "value"; d
}
zoom_coord <- function(stem) {
  e <- as.vector(terra::ext(terra::rast(zoom_file(stem, "s2rgb"))))
  coord_fixed(ratio = 1 / cos(mean(e[3:4]) * pi / 180), xlim = e[1:2], ylim = e[3:4], expand = FALSE)
}
zoom_ext <- function(stem) as.vector(terra::ext(terra::rast(zoom_file(stem, "s2rgb"))))
zoom_theme <- function() theme_map() + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                                             plot.title = element_text(size = FIG$base_size))
zoom_bar <- function(stem) {
  e <- zoom_ext(stem); w <- e[2] - e[1]; h <- e[4] - e[3]
  scalebar(e[1] + 0.05 * w, e[3] + 0.05 * h, 2, 0.012 * h, col = "white")
}

# ---- Inputs used by several figures -------------------------------------------
samples <- read_samples()
agc_s   <- samples |> dplyr::filter(is_agc)
nat_mean <- mean(agc_s$AGC_pred)
win <- FIG$windows |>
  dplyr::mutate(have = file.exists(zoom_file(stem, "s2rgb")) & file.exists(zoom_file(stem, "layers")),
                tag = letters[seq_len(dplyr::n()) + 1])
if (any(!win$have)) cat("Zoom windows without files:", paste(win$stem[!win$have], collapse = ", "), "\n")
win <- win |> dplyr::filter(have)
win_box <- purrr::map_dfr(seq_len(nrow(win)), function(i) {
  e <- zoom_ext(win$stem[i])
  tibble::tibble(stem = win$stem[i], label = win$label[i], tag = win$tag[i],
                 xmin = e[1], xmax = e[2], ymin = e[3], ymax = e[4])
})
ms <- if (file.exists(ofile("P4_metrics_summary.csv"))) rd(ofile("P4_metrics_summary.csv")) else NULL
if (!is.null(ms) && !"sample_set" %in% names(ms)) ms$sample_set <- "all"
mec_px <- if (!is.null(ms)) {
  m <- ms |> dplyr::filter(sample_set == "all", method == "Pixel-based GSE"); setNames(m$MEC_mean, m$design)
} else NULL
design_lab_mec <- function(d) if (is.null(mec_px) || length(mec_px) == 0) DESIGN_LABEL[d] else
  sprintf("%s: MEC %.2f", DESIGN_LABEL[d], mec_px[d])

# =============================================================================
# FIGURE 2: study area
# =============================================================================
fig("study", "study area")
meow_shp <- ofile("gee_upload/P4_meow_indo.shp")
if (need(grid_file("agcmean"), meow_shp, ofile("P4_ecoregion_summary.csv"))) {
  sf::sf_use_s2(FALSE)
  meow <- sf::st_read(meow_shp, quiet = TRUE)
  ga   <- rd(grid_file("agcmean"))
  ext_c <- agg_cells(ga, FIG$ext_cell, "area_m2") |> dplyr::mutate(area_km2 = area_m2 / 1e6)
  pts <- agc_s |> dplyr::distinct(xcoord, ycoord, year)
  p2a <- ggplot() +
    geom_sf(data = meow, fill = NA, colour = "grey40", linewidth = 0.2, linetype = 2) +
    land() +
    geom_tile(data = ext_c, aes(gx, gy, fill = area_km2), width = FIG$ext_cell, height = FIG$ext_cell) +
    geom_point(data = pts, aes(xcoord, ycoord, colour = factor(year)), size = 0.7, shape = 16) +
    scale_fill_gradientn(colours = c("#c6dbef", "#6baed6", "#2171b5", "#08306b"), trans = "log10",
                         breaks = c(0.01, 0.1, 1, 10), labels = c("0.01", "0.1", "1", "10"),
                         name = expression("Extent (km"^2*" per cell)")) +
    scale_colour_manual(values = YEAR_COL, name = "AGC survey year") +
    guides(colour = guide_legend(override.aes = list(size = 2.2))) +
    coord_sf(xlim = FIG$nat_xlim, ylim = FIG$nat_ylim, expand = FALSE) + nat_axes() +
    scalebar(96, -10.8, 500, 0.15) +
    labs(title = sprintf("a  Seagrass extent (%.0f %s) and AGC field samples (n = %d)",
                         sum(ga$area_m2) / 1e6, KM2, nrow(agc_s))) +
    theme_map() + theme(legend.position = "right")
  es <- rd(ofile("P4_ecoregion_summary.csv")) |>
    dplyr::filter(eco_id > 0) |>
    dplyr::select(ecoregion, `Extent area` = share_area_pct, `AGC samples` = share_samples_pct) |>
    tidyr::pivot_longer(-ecoregion, names_to = "what", values_to = "pct")
  ord <- es |> dplyr::filter(what == "Extent area") |> dplyr::arrange(pct) |> dplyr::pull(ecoregion)
  p2b <- ggplot(es, aes(x = factor(ecoregion, levels = ord), y = pct, fill = what)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = c("Extent area" = "#2171b5", "AGC samples" = "#d95f02"), name = NULL) +
    labs(title = "b  Share of extent area and of AGC samples", x = NULL, y = "Share (%)") +
    theme_p4() + theme(legend.position = "bottom")
  fig2 <- p2a / p2b + plot_layout(heights = c(1, 1))
  print(fig2)
  save_map("study", fig2, 190, 160)
}

# =============================================================================
# FIGURE 3: mean AGC 2017-2024, national grid + zoom windows
# =============================================================================
fig("agcmap", "mean AGC map")
if (need(grid_file("agcmean")) && nrow(win) > 0) {
  ga <- rd(grid_file("agcmean"))
  nat <- agg_cells(ga, FIG$nat_cell, c("area_m2", "agc_mean_g")) |> dplyr::mutate(agc = agc_mean_g / area_m2)
  p3a <- ggplot() + land() +
    geom_tile(data = nat, aes(gx, gy, fill = agc), width = FIG$nat_cell, height = FIG$nat_cell) +
    geom_rect(data = win_box, aes(xmin = xmin - 0.25, xmax = xmax + 0.25, ymin = ymin - 0.25, ymax = ymax + 0.25),
              fill = NA, colour = "black", linewidth = 0.4) +
    geom_text(data = win_box, aes(x = xmax + 0.45, y = ymax + 0.35, label = tag), size = 3, fontface = "bold") +
    agc_scale(name = expression("AGC (gC m"^-2*")")) +
    nat_coord() + scalebar(96, -10.8, 500, 0.15) +
    labs(title = sprintf("a  Mean AGC 2017-2024 from GSE, area-weighted per %.1f° cell", FIG$nat_cell)) +
    theme_map()
  wst <- purrr::map_dfr(seq_len(nrow(win)), function(i) {
    st <- win$stem[i]; e <- zoom_ext(st)
    d <- rast_df(zoom_file(st, "layers"), "AGC_mean")
    fm <- agc_s |> dplyr::filter(dplyr::between(xcoord, e[1], e[2]), dplyr::between(ycoord, e[3], e[4]))
    tibble::tibble(panel = win$tag[i], window = win$label[i], stem = st, n_field = nrow(fm),
                   field_mean_gC_m2 = if (nrow(fm) > 0) mean(fm$AGC_pred) else NA_real_,
                   n_map_px = nrow(d), map_mean_gC_m2 = mean(d$value))
  })
  print(as.data.frame(wst), row.names = FALSE, digits = 4)
  readr::write_csv(wst, fig_out("P4_Fig3_window_stats.csv"))
  zp <- purrr::map(seq_len(nrow(win)), function(i) {
    st <- win$stem[i]; e <- zoom_ext(st)
    d <- rast_df(zoom_file(st, "layers"), "AGC_mean")
    sp <- agc_s |> dplyr::filter(dplyr::between(xcoord, e[1], e[2]), dplyr::between(ycoord, e[3], e[4])) |>
      dplyr::distinct(xcoord, ycoord)
    fm <- agc_s |> dplyr::filter(dplyr::between(xcoord, e[1], e[2]), dplyr::between(ycoord, e[3], e[4]))
    ttl <- if (nrow(fm) > 0) sprintf("%s  %s (n = %d)", win$tag[i], win$label[i], nrow(fm)) else
      sprintf("%s  %s", win$tag[i], win$label[i])
    sub <- if (nrow(fm) > 0) sprintf("field %.1f | map %.1f", mean(fm$AGC_pred), mean(d$value)) else
      sprintf("no field data | map %.1f", mean(d$value))
    ggplot() + s2_layer(st) +
      geom_raster(data = d, aes(x, y, fill = value)) +
      geom_point(data = sp, aes(xcoord, ycoord), shape = 16, size = 0.25, colour = "red", alpha = 0.6) +
      agc_scale(guide = "none") +
      zoom_coord(st) + zoom_bar(st) +
      labs(title = ttl, subtitle = sub) +
      zoom_theme() + theme(plot.subtitle = element_text(size = FIG$base_size - 2))
  })
  fig3 <- p3a / wrap_plots(zp, nrow = 1) + plot_layout(heights = c(1.3, 1)) +
    plot_annotation(caption = paste0("b-e: 10 m pixels over Sentinel-2 (median 2020-2024); red dots: AGC field samples; n: field samples in the window.\n",
                                     "field | map: mean AGC (gC m-2) of the field samples | of the mapped pixels in the window."),
                    theme = theme(plot.caption = element_text(size = FIG$base_size - 2)))
  print(fig3)
  save_map("agcmap", fig3, 190, 135)
}

# =============================================================================
# FIGURE 4: totals by method (aim 1)
#   a national totals 2017-2024 with 95 % CI; b survey composition; c ecoregion totals
# =============================================================================
fig("totals", "national and ecoregion totals by method")
if (need(ofile("P4_national_three_methods.csv"), ofile("P4_ecoregion_summary.csv"))) {
  n3 <- rd(ofile("P4_national_three_methods.csv"))
  col4 <- c("Pixel-based GSE" = "#1b7837", "Tier 1 (pooled mean)" = "#762a83",
            "Tier 1 (survey-year mean)" = "#c2a5cf", "Tier 2 (ecoregion strata)" = "#e08214")
  n3 <- n3 |> dplyr::mutate(method = factor(method, levels = names(col4)))
  line_m <- n3 |> dplyr::filter(method != "Tier 1 (survey-year mean)")
  pt_m   <- n3 |> dplyr::filter(method == "Tier 1 (survey-year mean)")
  p4a <- ggplot(n3, aes(year, total_kt, colour = method, fill = method)) +
    geom_ribbon(data = line_m, aes(ymin = ci_low_kt, ymax = ci_up_kt), alpha = 0.12, colour = NA) +
    geom_line(data = line_m) + geom_point(data = line_m, size = 1) +
    geom_linerange(data = dplyr::filter(pt_m, !is.na(ci_low_kt)), aes(ymin = ci_low_kt, ymax = ci_up_kt),
                   position = position_nudge(x = 0.12)) +
    geom_point(data = pt_m, size = 1.6, position = position_nudge(x = 0.12)) +
    annotate("text", x = 2024.3, y = max(n3$ci_up_kt, na.rm = TRUE), label = "2024: no field data",
             size = 2.4, hjust = 1, vjust = 1, colour = "grey30") +
    scale_colour_manual(values = col4, name = NULL, drop = FALSE) +
    scale_fill_manual(values = col4, name = NULL, drop = FALSE) +
    scale_x_continuous(breaks = FIG$years) +
    guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2)) +
    labs(title = "a  National seagrass AGC by method (95 % CI)", x = NULL, y = "AGC (kt C)") +
    theme_p4() + theme(legend.position = "bottom")

  es <- rd(ofile("P4_ecoregion_summary.csv")) |> dplyr::filter(eco_id > 0)
  es <- es |> dplyr::mutate(eco_lab = sprintf("%s (%.0f %s, n = %d)", ecoregion, area_km2, KM2,
                                              as.integer(dplyr::coalesce(n_agc, 0))))
  ord <- es |> dplyr::arrange(area_km2) |> dplyr::pull(ecoregion)
  lab_of <- setNames(es$eco_lab, es$ecoregion)
  e4 <- dplyr::bind_rows(
    es |> dplyr::transmute(ecoregion, method = "Pixel-based GSE", kt = px_total_tC / 1e3,
                           lo = px_total_min_tC / 1e3, hi = px_total_max_tC / 1e3, shape = "estimate"),
    es |> dplyr::transmute(ecoregion, method = "Tier 1, national mean", kt = t1_national_factor_tC / 1e3,
                           lo = NA_real_, hi = NA_real_, shape = "estimate"),
    es |> dplyr::transmute(ecoregion, method = "Tier 2 (ecoregion strata)", kt = t2_total_tC / 1e3,
                           lo = (t2_total_tC - t2_ci_half_tC) / 1e3, hi = (t2_total_tC + t2_ci_half_tC) / 1e3,
                           shape = ifelse(t2_stratum %in% TRUE, "estimate", "national mean (n < 30)"))) |>
    dplyr::mutate(ecoregion = factor(ecoregion, levels = ord))
  col4b <- c("Pixel-based GSE" = "#1b7837", "Tier 1, national mean" = "#c2a5cf",
             "Tier 2 (ecoregion strata)" = "#e08214")
  p4b <- ggplot(e4, aes(kt, ecoregion, colour = method)) +
    geom_linerange(aes(xmin = lo, xmax = hi), na.rm = TRUE, position = position_dodge(width = 0.6)) +
    geom_point(aes(shape = shape), size = 1.6, na.rm = TRUE, position = position_dodge(width = 0.6)) +
    scale_y_discrete(labels = lab_of) +
    scale_colour_manual(values = col4b, name = NULL) +
    scale_shape_manual(values = c("estimate" = 16, "national mean (n < 30)" = 1), name = NULL) +
    labs(title = "c  Mean annual AGC per ecoregion", x = "AGC (kt C)", y = NULL) +
    guides(colour = guide_legend(nrow = 1, order = 1), shape = guide_legend(nrow = 1, order = 2)) +
    theme_p4() + theme(legend.position = "bottom", legend.box = "vertical", plot.title.position = "plot")

  comp <- agc_s |> dplyr::count(year, ecoregion)
  top  <- comp |> dplyr::group_by(ecoregion) |> dplyr::summarise(n = sum(n)) |> dplyr::arrange(dplyr::desc(n)) |>
    dplyr::slice_head(n = 6) |> dplyr::pull(ecoregion)
  comp <- comp |> dplyr::mutate(eco = ifelse(ecoregion %in% top, ecoregion, "Other")) |>
    dplyr::group_by(year, eco) |> dplyr::summarise(n = sum(n), .groups = "drop") |>
    dplyr::mutate(eco = factor(eco, levels = c(top, "Other")))
  p4c <- ggplot(comp, aes(factor(year), n, fill = eco)) + geom_col(width = 0.75) +
    scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "grey70"), name = NULL) +
    guides(fill = guide_legend(ncol = 2)) +
    labs(title = "b  AGC samples by survey year", x = NULL, y = "Samples") +
    theme_p4() + theme(legend.position = "bottom", legend.text = element_text(size = FIG$base_size - 2))
  # the ecoregion panel is wrapped so its long labels do not shift the panels above
  fig4 <- (p4a + p4c + plot_layout(widths = c(1.25, 1))) / wrap_elements(full = p4b) +
    plot_layout(heights = c(1, 1.25))
  print(fig4)
  save_chart("totals", fig4, 190, 215)
}

# =============================================================================
# FIGURE 5: inventory-level validation (aim 2)
#   a unit means, out-of-fold (unit held out) vs field, all AGC samples, area-based and GSE
#   b MEC of every method by unit level, all AGC samples and the Sentinel-2 set
# =============================================================================
fig("inventory", "inventory-level validation")
if (need(ofile("P4_inventory_units.csv"), ofile("P4_inventory_metrics.csv"))) {
  iu <- rd(ofile("P4_inventory_units.csv"))
  im <- rd(ofile("P4_inventory_metrics.csv"))
  meth_all <- intersect(c("Tier 1", "Tier 2", "Local mean", "Pixel-based GSE"), names(iu))
  pa <- iu |> dplyr::filter(sample_set == "all") |>
    tidyr::pivot_longer(dplyr::all_of(meth_all), names_to = "method", values_to = "pred") |>
    dplyr::filter(!is.na(pred)) |>
    dplyr::mutate(level = factor(LEVEL_LABEL[level], levels = LEVEL_LABEL),
                  method = factor(method, levels = meth_all))
  st <- im |> dplyr::filter(sample_set == "all", method %in% meth_all) |>
    dplyr::mutate(level = factor(LEVEL_LABEL[level], levels = LEVEL_LABEL),
                  method = factor(method, levels = meth_all),
                  txt = sprintf("MEC %.2f\nslope %.2f", MEC_mean, slope_mean))
  # common limits per row: field and predicted means of that level
  lims <- pa |> dplyr::group_by(level) |>
    dplyr::summarise(lo = min(c(obs, pred)), hi = max(c(obs, pred)), .groups = "drop")
  st <- st |> dplyr::left_join(lims, by = "level")
  # panels equal by definition (local mean = Tier 2 with the meadow held out; Tier 2 =
  # Tier 1 with the ecoregion held out) are left empty with a note
  dup <- pa |> dplyr::arrange(level, method, unit) |> dplyr::group_by(level, method) |>
    dplyr::summarise(sig = paste(sprintf("%.9g", pred), collapse = "|"), .groups = "drop") |>
    dplyr::arrange(level, method) |> dplyr::group_by(level, sig) |>
    dplyr::mutate(same_as = dplyr::first(as.character(method)), is_dup = dplyr::row_number() > 1) |>
    dplyr::ungroup() |> dplyr::filter(is_dup) |> dplyr::select(level, method, same_as) |>
    dplyr::left_join(lims, by = "level")
  pa <- pa |> dplyr::anti_join(dup, by = c("level", "method"))
  st <- st |> dplyr::anti_join(dup, by = c("level", "method"))
  blank <- lims |> tidyr::crossing(method = factor(meth_all, levels = meth_all)) |>
    tidyr::pivot_longer(c(lo, hi), values_to = "v")
  p5a <- ggplot(pa, aes(obs, pred)) +
    geom_blank(data = blank, aes(v, v)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
    geom_hline(yintercept = nat_mean, linetype = 3, colour = "firebrick") +
    geom_point(aes(size = n, colour = method), alpha = 0.55, shape = 16) +
    geom_text(data = st, aes(x = lo, y = hi, label = txt), hjust = 0, vjust = 1, size = 2.3,
              lineheight = 0.9, inherit.aes = FALSE) +
    geom_text(data = dup, aes(x = hi, y = lo, label = sprintf("same as %s", same_as)),
              hjust = 1, vjust = 0, size = 2.4, colour = "grey30", inherit.aes = FALSE) +
    facet_grid(level ~ method, scales = "free") +
    scale_colour_manual(values = METHOD_COL, guide = "none") +
    scale_size(range = c(0.5, 2.8), trans = "log10", breaks = c(10, 100, 1000),
               name = "Samples in unit\n(log scale)") +
    labs(title = "a  Unit means held out of training, all AGC samples",
         x = expression("Field mean AGC (gC m"^-2*")"), y = expression("Out-of-fold mean (gC m"^-2*")")) +
    theme_p4() + theme(legend.position = "right", strip.text = element_text(size = FIG$base_size - 1.5))

  mb <- im |>
    dplyr::mutate(level = factor(LEVEL_LABEL[level], levels = LEVEL_LABEL),
                  method = factor(method, levels = rev(names(METHOD_COL))),
                  set = factor(SET_LABEL[sample_set], levels = SET_LABEL))
  p5b <- ggplot(mb, aes(MEC_mean, method, colour = method, shape = set)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_linerange(aes(xmin = MEC_q025, xmax = MEC_q975), position = position_dodge(width = 0.6)) +
    geom_point(size = 1.8, position = position_dodge(width = 0.6)) +
    facet_wrap(~ level, nrow = 1) +
    scale_colour_manual(values = METHOD_COL, guide = "none") +
    scale_shape_manual(values = c(16, 2), name = NULL) +
    labs(title = "b  MEC of the unit means (mean and 2.5-97.5 percentiles over CV repeats)",
         x = "MEC", y = NULL) +
    theme_p4() + theme(legend.position = "bottom")
  fig5 <- p5a / p5b + plot_layout(heights = c(1.9, 1)) +
    plot_annotation(caption = wrap_cap(paste0(
      "Units: site = 100 m block, samples of its main survey year (n >= ", CFG$min_n_site, "), held out in the 100 m block design; ",
      "meadow = DBSCAN cluster (n >= ", CFG$min_n_meadow, "), cluster design; ecoregion (n >= ", CFG$min_n_eco, "), ecoregion design.\n",
      "Dashed: 1:1; dotted: national field mean. With the meadow held out the local mean equals Tier 2, ",
      "and with the ecoregion held out Tier 2 equals Tier 1, by definition.")),
      theme = theme(plot.caption = element_text(size = FIG$base_size - 2, hjust = 0)))
  print(fig5)
  save_chart("inventory", fig5, 190, 220)
}

# =============================================================================
# FIGURE 6: where new samples are needed (aim 3)
#   a distance to the nearest AGC sample (validation design of the largest share per cell)
#   b extent of each ecoregion by distance to AGC samples, ordered by area beyond 5 km
# =============================================================================
fig("sampling", "distance to samples and sampling priority")
p6 <- list()
if (need(grid_file("dist"))) {
  gd <- rd(grid_file("dist"))
  dcols <- sprintf("area_d%02d", 0:10)
  cells <- agg_cells(gd, FIG$nat_cell, dcols)
  by_des <- sapply(DESIGN_LEVELS, function(d) rowSums(as.matrix(cells[, dcols[DIST_DESIGN == d], drop = FALSE])))
  cells$design <- factor(DESIGN_LEVELS[max.col(by_des, ties.method = "first")], levels = DESIGN_LEVELS)
  share <- 100 * colSums(by_des) / sum(by_des)
  labs_des <- setNames(sprintf("%s; %s of extent", design_lab_mec(DESIGN_LEVELS), pct_lab(share)), DESIGN_LEVELS)
  p6$a <- ggplot() + land() +
    geom_tile(data = cells, aes(gx, gy, fill = design), width = FIG$nat_cell, height = FIG$nat_cell) +
    scale_fill_manual(values = DESIGN_COL, labels = labs_des, drop = FALSE,
                      name = "Distance to the nearest AGC sample (validation design):\nsample-level MEC of GSE; share of extent") +
    guides(fill = guide_legend(ncol = 2, title.position = "top")) +
    nat_coord() + scalebar(96, -10.8, 500, 0.15) +
    labs(title = sprintf("a  Distance class of the largest share of extent per %.1f° cell", FIG$nat_cell)) +
    theme_map() + theme(legend.position = "bottom")
}
if (need(ofile("P4_distance_coverage_eco.csv"), ofile("P4_ecoregion_summary.csv"))) {
  dce <- rd(ofile("P4_distance_coverage_eco.csv")) |> dplyr::filter(eco_id > 0)
  es6 <- rd(ofile("P4_ecoregion_summary.csv")) |>
    dplyr::select(ecoregion, n_agc, t2_stratum, share_area_pct, share_samples_pct)
  pr <- dce |> dplyr::left_join(es6, by = "ecoregion") |>
    dplyr::mutate(n_agc = as.integer(dplyr::coalesce(n_agc, 0)),
                  `Within 1 km` = within_1km_pct,
                  `1-5 km` = 100 - within_1km_pct - beyond_5km_pct,
                  `5-50 km` = beyond_5km_pct - beyond_50km_pct,
                  `> 50 km` = beyond_50km_pct,
                  beyond5_km2 = area_km2 * beyond_5km_pct / 100) |>
    dplyr::arrange(dplyr::desc(beyond5_km2))
  readr::write_csv(pr |> dplyr::select(ecoregion, area_km2, n_agc, t2_stratum, share_area_pct, share_samples_pct,
                                       within_1km_pct, beyond_5km_pct, beyond_50km_pct, beyond5_km2),
                   fig_out("P4_Fig6_ecoregion_priority.csv"))
  lab6 <- setNames(sprintf("%s (%.0f %s; n = %d%s)", pr$ecoregion, pr$area_km2, KM2, pr$n_agc,
                           ifelse(pr$t2_stratum %in% TRUE, "", ", no stratum")), pr$ecoregion)
  dl <- pr |> tidyr::pivot_longer(c(`Within 1 km`, `1-5 km`, `5-50 km`, `> 50 km`), names_to = "dist", values_to = "pct") |>
    dplyr::mutate(km2 = area_km2 * pct / 100,
                  dist = factor(dist, levels = c("> 50 km", "5-50 km", "1-5 km", "Within 1 km")),
                  ecoregion = factor(ecoregion, levels = rev(pr$ecoregion)))
  p6$b <- ggplot(dl, aes(km2, ecoregion, fill = dist)) +
    geom_col(width = 0.7) +
    scale_y_discrete(labels = lab6) +
    scale_fill_manual(values = c("Within 1 km" = "#1a9850", "1-5 km" = "#d9ef8b", "5-50 km" = "#fc8d59",
                                 "> 50 km" = "#d73027"), name = "Distance to the nearest AGC sample",
                      guide = guide_legend(reverse = TRUE)) +
    labs(title = "b  Extent of each ecoregion by distance to AGC samples (ordered by area beyond 5 km)",
         x = expression("Extent area (km"^2*")"), y = NULL) +
    theme_p4() + theme(legend.position = "bottom", plot.title.position = "plot")
}
if (length(p6) > 0) {
  # panel b is wrapped so its long ecoregion labels do not shrink the map
  fig6 <- if (length(p6) == 2) p6$a / wrap_elements(full = p6$b) + plot_layout(heights = c(1.05, 1)) else p6[[1]]
  fig6 <- fig6 + plot_annotation(caption = paste0(
    "n: AGC field samples; no stratum: fewer than ", CFG$min_n_eco, " samples, so Tier 2 uses the national mean."),
    theme = theme(plot.caption = element_text(size = FIG$base_size - 2, hjust = 0)))
  print(fig6)
  if (!is.null(p6$a)) save_map("sampling", fig6, 190, 185) else save_chart("sampling", fig6, 190, 100)
}

# =============================================================================
# SUPPLEMENT: presence probability threshold per year
# =============================================================================
fig("pa", "presence threshold by year")
if (need(ofile("P4_pa_threshold_by_year.csv"))) {
  pa_t <- rd(ofile("P4_pa_threshold_by_year.csv")) |> dplyr::mutate(year = as.character(year)) |>
    dplyr::filter(!is.na(best_F1)) |>
    tidyr::pivot_longer(c(best_F1, F1_050, F1_070), names_to = "what", values_to = "F1") |>
    dplyr::mutate(what = dplyr::recode(what, best_F1 = "best threshold", F1_050 = "threshold 0.5", F1_070 = "threshold 0.7"))
  figS1 <- ggplot(pa_t, aes(year, F1, colour = what, group = what)) + geom_point(size = 1.8) + geom_line() +
    scale_colour_manual(values = c("best threshold" = "grey30", "threshold 0.5" = "#2171b5", "threshold 0.7" = "#d73027"),
                        name = NULL) +
    labs(x = NULL, y = "F1 (presence and absence points)", caption = "Points of 2017 were not used to train the presence model.") +
    theme_p4() + theme(legend.position = "bottom")
  print(figS1)
  save_chart("pa", figS1, 90, 70)
}

# =============================================================================
# SUPPLEMENT: extent components and mask diagnostic
# =============================================================================
fig("extent", "extent components and mask diagnostic")
s2p <- list()
if (need(grid_file("agcmean"))) {
  ga <- rd(grid_file("agcmean"))
  ec <- agg_cells(ga, FIG$nat_cell, c("area_m2", "area_core_m2")) |>
    dplyr::mutate(addon_pct = 100 * (1 - area_core_m2 / area_m2))
  s2p$a <- ggplot() + land() +
    geom_tile(data = ec, aes(gx, gy, fill = addon_pct), width = FIG$nat_cell, height = FIG$nat_cell) +
    scale_fill_gradientn(colours = c("#2171b5", "#bdbdbd", "#e6550d"), limits = c(0, 100),
                         name = "Add-on share (%)") +
    nat_coord() +
    labs(title = sprintf("a  Final extent: core %.0f %s + add-on %.0f %s (shallow pixels without ACA geomorphic data)",
                         sum(ga$area_core_m2) / 1e6, KM2, (sum(ga$area_m2) - sum(ga$area_core_m2)) / 1e6, KM2)) +
    theme_map()
}
md_stem <- if (file.exists(zoom_file("Rote", "maskdiag")) && file.exists(zoom_file("Rote", "s2rgb"))) "Rote" else NA
if (!is.na(md_stem)) {
  M <- terra::rast(zoom_file(md_stem, "maskdiag"))
  d <- terra::as.data.frame(M, xy = TRUE, na.rm = TRUE) |>
    dplyr::filter(depth_raw > 0, depth_raw <= 500) |>
    dplyr::mutate(reason = dplyr::case_when(mask_code == 0 ~ "in extent", mask_code == 1 ~ "MaxProb < 0.7",
                                            mask_code == 4 ~ FIG$trf_label, TRUE ~ "two or more"),
                  reason = factor(reason, levels = c("in extent", "MaxProb < 0.7", FIG$trf_label, "two or more")))
  e <- zoom_ext(md_stem)
  sw <- samples |> dplyr::filter(dplyr::between(xcoord, e[1], e[2]), dplyr::between(ycoord, e[3], e[4])) |>
    dplyr::distinct(xcoord, ycoord)
  sh <- d |> dplyr::count(reason, .drop = FALSE) |> dplyr::mutate(p = 100 * n / sum(n))
  lab_r <- setNames(sprintf("%s (%.0f%%)", sh$reason, sh$p), as.character(sh$reason))
  bx <- c(range(d$x), range(d$y)); pad <- 0.08 * max(bx[2] - bx[1], bx[4] - bx[3])
  bx <- c(max(e[1], bx[1] - pad), min(e[2], bx[2] + pad), max(e[3], bx[3] - pad), min(e[4], bx[4] + pad))
  s2p$b <- ggplot() + s2_layer(md_stem) +
    geom_raster(data = d, aes(x, y, fill = reason), alpha = 0.75) +
    geom_point(data = sw, aes(xcoord, ycoord), shape = 3, size = 0.3) +
    scale_fill_manual(values = setNames(c("#1a9850", "#fdae61", "#d73027", "#762a83"),
                                        c("in extent", "MaxProb < 0.7", FIG$trf_label, "two or more")),
                      labels = lab_r, drop = FALSE, name = "Share of shallow area") +
    coord_fixed(ratio = 1 / cos(mean(bx[3:4]) * pi / 180), xlim = bx[1:2], ylim = bx[3:4], expand = FALSE) +
    scalebar(bx[1] + 0.06 * (bx[2] - bx[1]), bx[3] + 0.04 * (bx[4] - bx[3]), 1, 0.008 * (bx[4] - bx[3]), col = "white") +
    labs(title = "b  Rote: shallow pixels (ACA depth 0-5 m) by extent mask; + = field samples") +
    zoom_theme() + theme(legend.position = "right")
}
if (length(s2p) > 0) {
  figS2 <- if (length(s2p) == 2) s2p$a / (s2p$b + plot_spacer() + plot_layout(widths = c(1.3, 1))) +
    plot_layout(heights = c(1, 1.3)) else wrap_plots(s2p, ncol = 1)
  print(figS2)
  save_map("extent", figS2, 190, 170)
}

# =============================================================================
# SUPPLEMENT: extent under the geomorphic mask options vs KKP 2025
# =============================================================================
fig("kkp", "extent and the geomorphic mask")
if (need(ofile("P4_extent_sensitivity_ecoregion.csv"), ofile("P4_pixel_national.csv"))) {
  ee <- rd(ofile("P4_extent_sensitivity_ecoregion.csv"))
  final_km2 <- mean(rd(ofile("P4_pixel_national.csv"))$area_km2)
  core_km2  <- if (file.exists(grid_file("agcmean"))) sum(rd(grid_file("agcmean"))$area_core_m2) / 1e6 else NA_real_
  f30 <- core_km2 / sum(ee$current_km2)
  reg_km2 <- f30 * (sum(ee$only16_km2) + sum((ee$nogeo_km2 - ee$only16_km2)[!ee$ecoregion %in% FIG$turbid_eco]))
  s3 <- tibble::tibble(
    option = c("Core extent only", sprintf("Final: %s removed nationally", FIG$trf_label),
               sprintf("%s removed only in %s", FIG$trf_label, paste(FIG$turbid_eco, collapse = ", ")),
               "No geomorphic mask"),
    km2 = c(core_km2, final_km2, reg_km2, f30 * sum(ee$nogeo_km2)),
    kind = c("mapped", "mapped", "estimate", "estimate")) |>
    dplyr::mutate(option = factor(option, levels = rev(option)), pct_kkp = 100 * km2 / FIG$kkp_km2)
  print(as.data.frame(s3), row.names = FALSE, digits = 4)
  figS3 <- ggplot(s3, aes(km2, option, fill = kind)) + geom_col(width = 0.6) +
    geom_vline(xintercept = FIG$kkp_km2, linetype = 2, colour = "firebrick") +
    annotate("text", x = FIG$kkp_km2, y = 4.75, label = sprintf(" KKP 2025: %.0f %s", FIG$kkp_km2, KM2),
             colour = "firebrick", size = 2.5, hjust = 0, vjust = 0.5) +
    geom_text(aes(label = sprintf("%.0f %s (%.0f%% of KKP)", km2, KM2, pct_kkp)), hjust = -0.05, size = 2.4) +
    scale_fill_manual(values = c(mapped = "#2171b5", estimate = "#9ecae1"),
                      labels = c(mapped = "mapped (10 m)", estimate = "estimated from 30 m, scaled to 10 m"), name = NULL) +
    scale_x_continuous(limits = c(0, max(s3$km2, na.rm = TRUE) * 1.45), expand = expansion(mult = c(0, 0.02))) +
    scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
    labs(x = expression("Extent area (km"^2*")"), y = NULL) + theme_p4() +
    theme(legend.position = "bottom")
  print(figS3)
  save_chart("kkp", figS3, 170, 65)
}

# =============================================================================
# SUPPLEMENT: annual AGC maps
# =============================================================================
fig("annual", "annual AGC maps")
yf <- gfile(sprintf("P4_grid%dm_agc_%d.csv", FIG$grid_m, FIG$years))
if (need(yf)) {
  s1 <- purrr::map_dfr(FIG$years, function(y) {
    agg_cells(rd(gfile(sprintf("P4_grid%dm_agc_%d.csv", FIG$grid_m, y))), FIG$ext_cell, c("area_m2", "agc_g")) |>
      dplyr::mutate(year = y, agc = agc_g / area_m2)
  })
  tot <- s1 |> dplyr::group_by(year) |> dplyr::summarise(kt = sum(agc_g) / 1e9, .groups = "drop")
  figS4 <- ggplot() + land() +
    geom_tile(data = s1, aes(gx, gy, fill = agc), width = FIG$ext_cell, height = FIG$ext_cell) +
    geom_text(data = tot, aes(x = 140.5, y = 5.5, label = sprintf("%.0f kt C", kt)), hjust = 1, size = 2.3) +
    facet_wrap(~ year, ncol = 2) +
    agc_scale(name = expression("AGC (gC m"^-2*")")) +
    guides(fill = guide_colourbar(barwidth = unit(60, "mm"), title.position = "top", title.hjust = 0.5)) +
    nat_coord() + theme_map() + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                                      legend.position = "bottom")
  print(figS4)
  save_map("annual", figS4, 190, 175, dpi_png = 300)
}

# =============================================================================
# SUPPLEMENT: temporal behaviour (main Fig 8 in P4v2)
# =============================================================================
fig("temporal", "temporal behaviour")
p9 <- list()
if (need(grid_file("cv"))) {
  gc <- rd(grid_file("cv"))
  cvc <- agg_cells(gc, FIG$nat_cell, c("area_cv_m2", "cv_x_area")) |> dplyr::mutate(cv = cv_x_area / area_cv_m2)
  p9$a <- ggplot() + land() +
    geom_tile(data = cvc, aes(gx, gy, fill = cv), width = FIG$nat_cell, height = FIG$nat_cell) +
    scale_fill_gradientn(colours = c("#f7fbff", "#9ecae1", "#4292c6", "#08519c", "#08306b"),
                         limits = c(0, 30), oob = scales::squish, name = "CV (%)") +
    nat_coord() + scalebar(96, -10.8, 500, 0.15) +
    labs(title = sprintf("a  AGC CV between years (national mean %.1f %%)",
                         sum(gc$cv_x_area) / sum(gc$area_cv_m2))) +
    theme_map() + theme(legend.position = "right")
}
if (need(ofile("P4_gse_similarity_pairs.csv"))) {
  sp <- rd(ofile("P4_gse_similarity_pairs.csv")) |> dplyr::filter(grepl("^sim_\\d{4}_\\d{4}$", pair)) |>
    dplyr::mutate(pair = sub("sim_(\\d{4})_(\\d{4})", "\\1-\\2", pair))
  p9$b <- ggplot(sp, aes(x = pair)) +
    geom_linerange(aes(ymin = p05, ymax = p95), colour = "grey60") +
    geom_crossbar(aes(y = p50, ymin = p25, ymax = p75), width = 0.5, fill = "grey90", linewidth = 0.3) +
    labs(title = "b  GSE similarity", x = "Consecutive years", y = "Dot product (1 = identical)") +
    theme_p4() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
ann <- win |> dplyr::filter(file.exists(zoom_file(stem, "agc_annual")))
if (nrow(ann) > 0) {
  # one line per window: the samples of its best-sampled survey year (main cohort)
  yr9 <- purrr::map_dfr(seq_len(nrow(ann)), function(i) {
    st <- ann$stem[i]; A <- terra::rast(zoom_file(st, "agc_annual")); e <- as.vector(terra::ext(A))
    sw <- agc_s |> dplyr::filter(dplyr::between(xcoord, e[1], e[2]), dplyr::between(ycoord, e[3], e[4]))
    if (nrow(sw) == 0) return(NULL)
    v <- terra::extract(A, as.matrix(sw[, c("xcoord", "ycoord")]))
    ok <- !is.na(v[[1]]); sw <- sw[ok, ]; v <- v[ok, , drop = FALSE]
    purrr::map_dfr(sort(unique(sw$year)), function(y) {
      k <- sw$year == y
      if (sum(k) < 20) return(NULL)
      tibble::tibble(label = ann$label[i], cohort = factor(y), year = FIG$years,
                     map = colMeans(v[k, , drop = FALSE]), field = mean(sw$AGC_pred[k]),
                     survey = FIG$years == y, n = sum(k))
    })
  })
  if (nrow(yr9) > 0) yr9 <- yr9 |> dplyr::group_by(label) |> dplyr::filter(n == max(n)) |>
    dplyr::filter(cohort == dplyr::first(cohort)) |> dplyr::ungroup()
  wsum <- purrr::map_dfr(seq_len(nrow(ann)), function(i) {
    st <- ann$stem[i]; A <- terra::rast(zoom_file(st, "agc_annual")); e <- as.vector(terra::ext(A))
    sw <- agc_s |> dplyr::filter(dplyr::between(xcoord, e[1], e[2]), dplyr::between(ycoord, e[3], e[4]))
    if (nrow(sw) == 0) return(NULL)
    v <- as.matrix(terra::extract(A, as.matrix(sw[, c("xcoord", "ycoord")])))
    ok <- !is.na(v[, 1]); sw <- sw[ok, ]; v <- v[ok, , drop = FALSE]
    main_y <- as.integer(names(which.max(table(sw$year))))
    k <- sw$year == main_y; sw <- sw[k, ]; v <- v[k, , drop = FALSE]
    jy <- match(sw$year, FIG$years)
    own <- v[cbind(seq_len(nrow(v)), jy)]
    oth <- vapply(seq_len(nrow(v)), function(k) mean(v[k, -jy[k]]), numeric(1))
    tibble::tibble(window = ann$label[i], stem = st, n = nrow(sw), survey_years = paste(sort(unique(sw$year)), collapse = ", "),
                   field_mean_gC_m2 = mean(sw$AGC_pred), map_survey_year_gC_m2 = mean(own),
                   map_other_years_gC_m2 = mean(oth), map_all_years_gC_m2 = mean(v))
  })
  if (nrow(wsum) > 0) {
    print(as.data.frame(wsum), row.names = FALSE, digits = 4)
    readr::write_csv(wsum, fig_out(sprintf("P4_Fig%s_window_summary.csv", FIGNUM[["temporal"]])))
  }
  if (nrow(yr9) > 0) {
    readr::write_csv(yr9 |> dplyr::mutate(cohort = as.character(cohort)),
                     fig_out(sprintf("P4_Fig%s_window_years.csv", FIGNUM[["temporal"]])))
    yr9 <- yr9 |> dplyr::mutate(cohort = as.character(cohort)) |> dplyr::arrange(label, cohort, year) |>
      dplyr::mutate(panel = sprintf("%s, %s samples\n(n = %d)", label, cohort, n))
    yr9$panel <- factor(yr9$panel, levels = unique(yr9$panel))
    p9$c <- ggplot(yr9, aes(year, map)) +
      geom_hline(aes(yintercept = field), linetype = 2, linewidth = 0.4, colour = "grey35") +
      geom_hline(yintercept = nat_mean, linetype = 3, colour = "firebrick") +
      geom_line(colour = "grey20") +
      geom_point(data = dplyr::filter(yr9, !survey), shape = 1, size = 1.4, colour = "grey20") +
      geom_point(data = dplyr::filter(yr9, survey), shape = 16, size = 2.4, colour = "#2171b5") +
      facet_wrap(~ panel, nrow = 1) +
      scale_x_continuous(breaks = FIG$years, labels = substr(FIG$years, 3, 4)) +
      labs(title = "c  Mapped AGC at the sample pixels of the main survey year, by map year",
           caption = "Blue point: survey year of the samples; dashed: their field mean; dotted: national field mean.",
           x = "Map year", y = expression("AGC (gC m"^-2*")")) +
      theme_p4()
  }
}
if (length(p9) > 0) {
  figS5 <- if (all(c("a", "b", "c") %in% names(p9))) (p9$a + p9$b + plot_layout(widths = c(2.6, 1))) / p9$c +
    plot_layout(heights = c(0.85, 1)) else wrap_plots(p9, ncol = 1)
  print(figS5)
  save_map("temporal", figS5, 190, 140)
}

# =============================================================================
# SUPPLEMENT: sample-level accuracy against distance (main Fig 4 in P4v2), both sample sets
# =============================================================================
fig("samplelevel", "sample-level accuracy against distance")
if (!is.null(ms)) {
  ml <- dplyr::bind_rows(lapply(c("RMSE", "MAE", "R2", "MEC", "bias_pct"), function(k) {
    ms |> dplyr::transmute(sample_set, design, method, median_dist_m, metric = k,
                           mean = .data[[paste0(k, "_mean")]], lo = .data[[paste0(k, "_q025")]],
                           hi = .data[[paste0(k, "_q975")]])
  })) |>
    dplyr::mutate(metric = factor(metric, levels = c("RMSE", "MAE", "R2", "MEC", "bias_pct"),
                                  labels = c("RMSE~(gC~m^-2)", "MAE~(gC~m^-2)", "R^2", "MEC", "Bias~('%')")),
                  method = factor(method, levels = names(METHOD_COL)),
                  set = factor(sample_set, levels = names(SET_LABEL),
                               labels = c("All AGC samples", "Sentinel-2 set (2019-2023)")))
  brk <- ms |> dplyr::filter(sample_set == "all", method == "Pixel-based GSE") |>
    dplyr::mutate(design = factor(design, levels = DESIGN_LEVELS)) |> dplyr::arrange(design)
  brk_lab <- sprintf("%s\n%s", c("Pixel", "100 m", "1 km", "Meadow", "Eco-\nregion")[as.integer(brk$design)],
                     ifelse(brk$median_dist_m < 1000, sprintf("%.0f m", brk$median_dist_m),
                            sprintf("%.1f km", brk$median_dist_m / 1000)))
  figS6 <- ggplot(ml, aes(median_dist_m, mean, colour = method, group = method)) +
    geom_hline(data = data.frame(metric = factor(c("MEC", "Bias~('%')"), levels = levels(ml$metric)), y = 0),
               aes(yintercept = y), linetype = 2, colour = "grey40") +
    geom_line() + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.05) + geom_point(size = 1.3) +
    facet_grid(metric ~ set, scales = "free_y", labeller = labeller(metric = label_parsed)) +
    scale_x_log10(breaks = brk$median_dist_m, labels = brk_lab, expand = expansion(mult = 0.08)) +
    scale_colour_manual(values = METHOD_COL, name = NULL) +
    labs(x = "Validation design and median distance from test to nearest training sample (log scale)",
         y = "Mean over repeats (2.5-97.5 percentiles)",
         caption = paste0("The local mean equals Tier 2 in the meadow cluster design, and Tier 2 equals Tier 1 ",
                          "in the ecoregion design, so their lines overlap there.")) +
    theme_p4() + theme(legend.position = "bottom", axis.text.x = element_text(size = FIG$base_size - 2.5))
  print(figS6)
  save_chart("samplelevel", figS6, 190, 230)
}

# =============================================================================
# SUPPLEMENT: error by observed AGC (shrinkage of the GSE predictions toward the mean).
# The group means of meadows and ecoregions are in Fig 5 (inventory level).
# =============================================================================
fig("shrinkage", "prediction by observed AGC bin")
if (need(ofile("P4_error_by_agc_bin.csv"))) {
  eb <- rd(ofile("P4_error_by_agc_bin.csv")) |>
    dplyr::mutate(lo = as.numeric(sub("^\\[(\\d+),.*", "\\1", bin)), design = factor(design, levels = DESIGN_LEVELS))
  figS7 <- ggplot(eb, aes(lo + 5, mean_pred, colour = design)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
    geom_hline(yintercept = nat_mean, linetype = 3, colour = "firebrick") +
    geom_line() + geom_point(aes(size = n_per_repeat)) +
    scale_colour_manual(values = DESIGN_COL[c("pixel", "cluster")], labels = DESIGN_LABEL, name = NULL) +
    scale_size_continuous(range = c(0.6, 3), name = "Samples per bin") +
    labs(x = expression("Observed AGC, bin centre (gC m"^-2*")"),
         y = expression("Mean out-of-fold GSE prediction (gC m"^-2*")"),
         caption = "Bin width 10 gC m-2. Dashed: 1:1; dotted: national field mean.") + theme_p4()
  print(figS7)
  save_chart("shrinkage", figS7, 170, 85)
}

# =============================================================================
# SUPPLEMENT: cumulative extent by distance to AGC samples
# =============================================================================
fig("cumdist", "cumulative extent by distance")
if (need(ofile("P4_distance_coverage_class.csv"))) {
  dc <- rd(ofile("P4_distance_coverage_class.csv")) |>
    dplyr::mutate(dist_label = factor(dist_label, levels = dist_label), design = factor(design, levels = DESIGN_LEVELS))
  figS8 <- ggplot(dc, aes(dist_label, area_cum_pct, group = 1)) +
    geom_col(aes(y = area_pct, fill = design), width = 0.7) + geom_line() + geom_point(size = 1) +
    scale_fill_manual(values = DESIGN_COL, labels = DESIGN_LABEL, name = "Validation design") +
    labs(x = "Distance to the nearest AGC sample", y = "Share of extent (%)",
         caption = "Bars: share per distance class; line: cumulative share") +
    theme_p4() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  print(figS8)
  save_chart("cumdist", figS8, 160, 80)
}

cat("\nP4_08_figures.R COMPLETE\n")
