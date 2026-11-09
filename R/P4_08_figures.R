# =============================================================================
# P4_08_figures.R
# Paper 4 (revised design), added 7 Nov 2026: all figures of the manuscript and
# the Supplement, built from the R results (P4_04, P4_05, P4_07) and the GEE
# figure exports (P4_G08 grid, P4_G08b zoom windows).
#
# Update 9 Nov 2026: figures renumbered in the order of
# their first citation; output goes to a new folder <out_dir>/P4_figures_v2, so
# nothing written by the 7 Nov version is overwritten. Other changes: "Expected
# skill" -> "Expected performance" (Fig 5); the ACA geomorphic class removed
# from the extent is labelled "Terrestrial Reef Flat" (Fig S2b, S4); new panel
# S3c (cover error by 10 % bin, P4_cover_error_by_bin.csv); Fig 3 and Fig 9c
# write their window numbers to CSV, so the text can be traced (P4_10).
#
# Update 9 Nov 2026b: RSE submission ggsave lines added under each PNG line
# (TIFF 600 dpi for maps, vector PDF for charts). Figures are unchanged.
#
# Main figures
#   Fig 1  Study area: ecoregions, extent, AGC samples; extent share vs sample share
#   Fig 2  Workflow diagram (drawn in the manuscript, not by this script)
#   Fig 3  Mean AGC 2017-2024: national grid + four 10 m zoom windows       (was Fig 2)
#   Fig 4  RMSE, R2 and MEC against distance to the nearest training sample (was Fig 3)
#   Fig 5  Expected performance: extent by distance to AGC samples          (was Fig 4)
#   Fig 6  National totals 2017-2024 (three methods) + survey composition   (was Fig 5)
#   Fig 7  Ecoregions: totals by method + leave-one-ecoregion-out bias      (was Fig 6)
#   Fig 8  22 ROIs: field mean vs map mean, final map and out-of-fold        (was Fig 7)
#   Fig 9  Temporal: AGC CV map, GSE similarity between years, AGC by year at
#          the sample pixels of the zoom windows                             (was Fig 8)
# Supplement (was -> now)
#   S1 PA threshold by year (S6); S2 extent components and mask diagnostic (S2);
#   S3 error by AGC bin, shrinkage of group means, cover error by bin (S4 + new c);
#   S4 extent under the geomorphic mask options vs the KKP 2025 area (S7);
#   S5 annual AGC maps (S1); S6 cover class confusion (S5);
#   S7 cumulative extent by distance to AGC samples (S8);
#   S8 cover classes (Tier 2 strata) in the zoom windows (S3)
#
# Reads (out_dir): P4_samples.csv, P4_metrics_summary.csv,
#   P4_national_three_methods.csv, P4_ecoregion_summary.csv,
#   P4_lolo_region_means.csv, P4_roi_comparison.csv, P4_shrinkage_roi.csv,
#   P4_distance_coverage_class.csv, P4_distance_coverage_design.csv,
#   P4_error_by_agc_bin.csv, P4_shrinkage_groups.csv, P4_cover_error_by_bin.csv,
#   P4_cover_confusion.csv,
#   P4_pa_threshold_by_year.csv, P4_gse_similarity_pairs.csv,
#   P4_extent_sensitivity_ecoregion.csv, P4_extent_sensitivity_national.csv,
#   P4_pixel_national.csv, gee_upload/P4_meow_indo.shp
# Reads (gee_dir): P4_grid500m_agcmean.csv, P4_grid500m_dist.csv,
#   P4_grid500m_cv.csv, P4_grid500m_agc_YYYY.csv,
#   P4_zoom_<site>_layers.tif, _s2rgb.tif, _agc_annual.tif, _maskdiag.tif
# Writes (FIG$dir = <out_dir>/P4_figures_v2): P4_Fig3_window_stats.csv,
#   P4_Fig9c_window_years.csv, P4_Fig9c_window_summary.csv; figures only when
#   the ggsave lines are uncommented.
# A missing input only skips the figures that need it.
# All figures are printed; ggsave lines are commented (sizes in mm for RSE:
# single column 90 mm, double column 190 mm). Each figure has a PNG line for
# preview and an RSE submission line (Guide for Authors, checked 9 Nov 2026):
# maps as LZW TIFF at 600 dpi (Fig 1, 3, 5, 9, S2, S5, S8), charts as vector
# PDF through cairo_pdf (Fig 4, 6, 7, 8, S1, S3, S4, S6, S7).
#
# Expected warnings (harmless; the figures are correct):
#   "Raster pixels are placed at uneven horizontal intervals" - zoom windows:
#     columns without seagrass pixels leave gaps between the 10 m columns;
#     every pixel is still drawn at its own position.
#   "Removed 1 rows containing missing values (geom_segment)" - one ROI label
#     line in Fig 8 that ggrepel does not need.
# =============================================================================

# setwd("path/to/seagrass-agc-countrywide-indonesia/R")   # folder that holds these scripts
source("P4_00_config_functions.R")
cat("P4_08_figures.R | version 2026-11-09b (RSE file formats)\n")

for (p in c("terra", "patchwork", "ggrepel", "maps")) {
  if (!requireNamespace(p, quietly = TRUE)) stop(sprintf("Package '%s' is needed: install.packages('%s')", p, p))
}
suppressPackageStartupMessages({ library(terra); library(patchwork) })

# -----------------------------------------------------------------------------
# FIGURE CONFIG -- edit here only
# -----------------------------------------------------------------------------
FIG <- list(
  dir        = out_file("P4_figures_v2"),    # new folder (9 Nov 2026): figures and window CSVs
  gee_dir    = if (!is.null(CFG$gee_dir)) CFG$gee_dir else out_file("gee"),
  grid_m     = 500,                    # cell size of the P4_G08 grid exports
  nat_cell   = 0.1,                    # degrees, national map cells (Fig 3, 5, 9, S2)
  ext_cell   = 0.25,                   # degrees, extent cells in Fig 1 and annual maps in S5
  nat_xlim   = c(94.5, 141.5),
  nat_ylim   = c(-11.5, 6.5),
  years      = 2017:2024,
  agc_lim    = c(15, 40),              # colour limits of AGC maps (gC m-2), Fig 3 and S5. Fitted to
                                       # the mapped range (0.1 deg cells 16-37; 2-98 % of the annual
                                       # 0.25 deg cells 22-36); values outside are drawn at the end
                                       # colour and the legend ends read "<= 15" and ">= 40" (ASCII:
                                       # the PDF device cannot draw the Unicode signs)
  s2_fact    = 2,                      # aggregate Sentinel-2 background (2 = 20 m) for speed
  kkp_km2    = 6601.56,                # KKP 2025 national seagrass map (official area)
  turbid_eco = c("Sunda Shelf/Java Sea"),   # S4: Terrestrial Reef Flat removed only here (future work option)
  trf_label  = "Terrestrial Reef Flat",     # ACA geomorphic class removed from the extent (code 16)
  # zoom windows (P4_G08b): file stem, label shown in the figures
  windows = tibble::tribble(
    ~stem,      ~label,
    "Rote2",    "Rote",
    "Derawan2", "Derawan",
    "Seribu",   "Seribu Islands",
    "Selayar",  "Selayar"),
  # 22 ROIs (Supplementary Table S4)
  roi22 = sprintf("ROI-%03d", c(95, 70, 272, 115, 83, 210, 192, 317, 153, 214, 320,
                                41, 149, 3, 156, 57, 96, 159, 148, 103, 53, 102)),
  base_size = 9
)

DESIGN_LEVELS <- c("pixel", "b100", "b1000", "cluster", "ecoregion")
DESIGN_LABEL  <- c(pixel = "Pixel (<= 10 m)", b100 = "100 m block (10-100 m)",
                   b1000 = "1 km block (100 m-1 km)", cluster = "Meadow cluster (1-20 km)",
                   ecoregion = "Ecoregion (> 20 km)")
DESIGN_COL    <- c(pixel = "#1a9850", b100 = "#91cf60", b1000 = "#d9ef8b",
                   cluster = "#fc8d59", ecoregion = "#d73027")
DIST_DESIGN   <- c("pixel", "b100", "b100", "b1000", "b1000", "b1000",
                   "cluster", "cluster", "cluster", "ecoregion", "ecoregion")   # classes 0-10
METHOD_COL    <- c("Pixel-based" = "#1b7837", "Tier 1" = "#762a83", "Tier 2" = "#e08214")
# AGC maps: viridis (perceptually uniform, readable with colour-vision deficiency, lightness
# rises with AGC). Changed 7 Nov 2026 from a yellow-green ramp on 0-60, which showed the
# narrow mapped range (mostly 22-34 gC m-2) as nearly one colour.
AGC_PAL       <- c("#440154", "#472D7B", "#3B528B", "#2C728E", "#21908C", "#27AD81", "#5DC863", "#AADC32", "#FDE725")
# Legend of the AGC maps: end labels show that values beyond the limits use the end colours
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
fig <- function(txt) cat(sprintf("\n==== %s ====\n", txt))
dir.create(FIG$dir, showWarnings = FALSE, recursive = TRUE)
fig_out <- function(n) file.path(FIG$dir, n)      # all figure files and CSVs of this script
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
# Axis labels in degrees (plain Unicode degree sign, so no plotmath symbol font is needed)
lon_lab <- function(x) ifelse(x == 0, "0\u00b0", paste0(abs(x), "\u00b0", ifelse(x < 0, "W", "E")))
lat_lab <- function(y) ifelse(y == 0, "0\u00b0", paste0(abs(y), "\u00b0", ifelse(y < 0, "S", "N")))
nat_axes  <- function() list(scale_x_continuous(labels = lon_lab), scale_y_continuous(labels = lat_lab))
nat_coord <- function() c(list(coord_quickmap(xlim = FIG$nat_xlim, ylim = FIG$nat_ylim, expand = FALSE)), nat_axes())
pct_lab   <- function(p) ifelse(p < 0.05, "< 0.1%", sprintf("%.1f%%", p))
KM2 <- "km\u00b2"

# Scale bar for lon/lat panels (km at the latitude of the bar)
scalebar <- function(x0, y0, km, h, col = "black", size = 2.4) {
  dx <- km / (111.32 * cos(y0 * pi / 180))
  list(annotate("rect", xmin = x0, xmax = x0 + dx, ymin = y0, ymax = y0 + h, fill = col),
       annotate("text", x = x0 + dx / 2, y = y0 + 1.8 * h, label = paste(km, "km"), size = size, colour = col,
                vjust = 0))
}

# Sentinel-2 background of a zoom window as a colour matrix (2-98 % linear stretch)
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
mec_px <- if (!is.null(ms)) {
  m <- ms |> dplyr::filter(method == "Pixel-based"); setNames(m$MEC_mean, m$design)
} else NULL
design_lab_mec <- function(d) if (is.null(mec_px)) DESIGN_LABEL[d] else
  sprintf("%s: MEC %.2f", DESIGN_LABEL[d], mec_px[d])

# =============================================================================
# FIGURE 1: study area
# =============================================================================
fig("Figure 1: study area")
meow_shp <- ofile("gee_upload/P4_meow_indo.shp")
if (need(grid_file("agcmean"), meow_shp, ofile("P4_ecoregion_summary.csv"))) {
  sf::sf_use_s2(FALSE)
  meow <- sf::st_read(meow_shp, quiet = TRUE)
  ga   <- rd(grid_file("agcmean"))
  ext_c <- agg_cells(ga, FIG$ext_cell, "area_m2") |> dplyr::mutate(area_km2 = area_m2 / 1e6)
  pts <- agc_s |> dplyr::distinct(xcoord, ycoord, year)
  p1a <- ggplot() +
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
  p1b <- ggplot(es, aes(x = factor(ecoregion, levels = ord), y = pct, fill = what)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = c("Extent area" = "#2171b5", "AGC samples" = "#d95f02"), name = NULL) +
    labs(title = "b  Share of extent area and of AGC samples", x = NULL, y = "Share (%)") +
    theme_p4() + theme(legend.position = "bottom")
  fig1 <- p1a / p1b + plot_layout(heights = c(1, 1))
  print(fig1)
  # ggsave(fig_out("P4_Fig1_study_area.png"), fig1, width = 190, height = 160, units = "mm", dpi = 400)
  # ggsave(fig_out("P4_Fig1_study_area.tif"), fig1, width = 190, height = 160, units = "mm", dpi = 600, device = "tiff", compression = "lzw")   # RSE submission (map, 600 dpi)
}

# =============================================================================
# FIGURE 3: mean AGC 2017-2024, national grid + zoom windows
# =============================================================================
fig("Figure 3: mean AGC map")
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
    labs(title = sprintf("a  Mean AGC 2017-2024, area-weighted per %.1f\u00b0 cell", FIG$nat_cell)) +
    theme_map()
  # Numbers behind the panel subtitles, written for tracing the text (P4_10)
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
  # ggsave(fig_out("P4_Fig3_agc_map.png"), fig3, width = 190, height = 135, units = "mm", dpi = 400)
  # ggsave(fig_out("P4_Fig3_agc_map.tif"), fig3, width = 190, height = 135, units = "mm", dpi = 600, device = "tiff", compression = "lzw")   # RSE submission (map, 600 dpi)
}

# =============================================================================
# FIGURE 4: RMSE, R2 and MEC against distance
# =============================================================================
fig("Figure 4: accuracy against distance")
if (!is.null(ms)) {
  ml <- dplyr::bind_rows(lapply(c("RMSE", "R2", "MEC"), function(k) {
    ms |> dplyr::transmute(design, method, median_dist_m, metric = k,
                           mean = .data[[paste0(k, "_mean")]], lo = .data[[paste0(k, "_q025")]],
                           hi = .data[[paste0(k, "_q975")]])
  })) |>
    dplyr::mutate(metric = factor(metric, levels = c("RMSE", "R2", "MEC"),
                                  labels = c("RMSE~(gC~m^-2)", "R^2", "MEC")),
                  method = factor(method, levels = names(METHOD_COL)))
  brk <- ms |> dplyr::filter(method == "Pixel-based") |>
    dplyr::mutate(design = factor(design, levels = DESIGN_LEVELS)) |> dplyr::arrange(design)
  brk_lab <- sprintf("%s\n%s", c("Pixel", "100 m\nblock", "1 km\nblock", "Meadow\ncluster", "Eco-\nregion")[as.integer(brk$design)],
                     ifelse(brk$median_dist_m < 1000, sprintf("%.0f m", brk$median_dist_m),
                            sprintf("%.1f km", brk$median_dist_m / 1000)))
  fig4 <- ggplot(ml, aes(median_dist_m, mean, colour = method, group = method)) +
    geom_hline(data = data.frame(metric = factor("MEC", levels = levels(ml$metric)), y = 0),
               aes(yintercept = y), linetype = 2, colour = "grey40") +
    geom_line() + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.05) + geom_point(size = 1.6) +
    facet_wrap(~ metric, nrow = 1, scales = "free_y", labeller = label_parsed) +
    scale_x_log10(breaks = brk$median_dist_m, labels = brk_lab, expand = expansion(mult = 0.08)) +
    scale_colour_manual(values = METHOD_COL, name = NULL) +
    labs(x = "Validation design and median distance from test to nearest training sample (log scale)",
         y = "Mean over repeats (2.5-97.5 percentiles)") +
    theme_p4() + theme(legend.position = "bottom", axis.text.x = element_text(size = FIG$base_size - 2.5))
  print(fig4)
  # ggsave(fig_out("P4_Fig4_accuracy_distance.png"), fig4, width = 190, height = 85, units = "mm", dpi = 400)
  # ggsave(fig_out("P4_Fig4_accuracy_distance.pdf"), fig4, width = 190, height = 85, units = "mm", device = cairo_pdf)   # RSE submission (vector)
}

# =============================================================================
# FIGURE 5: expected performance (distance to AGC samples)
# =============================================================================
fig("Figure 5: expected performance")
if (need(grid_file("dist"))) {
  gd <- rd(grid_file("dist"))
  dcols <- sprintf("area_d%02d", 0:10)
  cells <- agg_cells(gd, FIG$nat_cell, dcols)
  by_des <- sapply(DESIGN_LEVELS, function(d) rowSums(as.matrix(cells[, dcols[DIST_DESIGN == d]])))
  cells$design <- factor(DESIGN_LEVELS[max.col(by_des, ties.method = "first")], levels = DESIGN_LEVELS)
  share <- 100 * colSums(by_des) / sum(by_des)
  labs_des <- setNames(sprintf("%s; %s of extent", design_lab_mec(DESIGN_LEVELS), pct_lab(share)), DESIGN_LEVELS)
  p5a <- ggplot() + land() +
    geom_tile(data = cells, aes(gx, gy, fill = design), width = FIG$nat_cell, height = FIG$nat_cell) +
    geom_rect(data = win_box, aes(xmin = xmin - 0.25, xmax = xmax + 0.25, ymin = ymin - 0.25, ymax = ymax + 0.25),
              fill = NA, colour = "black", linewidth = 0.4) +
    geom_text(data = win_box, aes(x = xmax + 0.45, y = ymax + 0.35, label = tag), size = 3, fontface = "bold") +
    scale_fill_manual(values = DESIGN_COL, labels = labs_des, drop = FALSE,
                      name = "Distance to nearest AGC sample (validation design): pixel-based MEC; share of extent") +
    guides(fill = guide_legend(ncol = 2, title.position = "top")) +
    nat_coord() + scalebar(96, -10.8, 500, 0.15) +
    labs(title = sprintf("a  Expected performance: validation design of the largest share of extent per %.1f\u00b0 cell", FIG$nat_cell)) +
    theme_map() + theme(legend.position = "bottom")
  zp5 <- if (nrow(win) > 0) purrr::map(seq_len(nrow(win)), function(i) {
    st <- win$stem[i]
    d <- rast_df(zoom_file(st, "layers"), "dist_class") |>
      dplyr::mutate(design = factor(DIST_DESIGN[value + 1], levels = DESIGN_LEVELS))
    ggplot() + s2_layer(st) +
      geom_raster(data = d, aes(x, y, fill = design)) +
      scale_fill_manual(values = DESIGN_COL, guide = "none", drop = FALSE) +
      zoom_coord(st) + zoom_bar(st) +
      labs(title = sprintf("%s  %s", win$tag[i], win$label[i])) + zoom_theme()
  }) else list()
  fig5 <- p5a / wrap_plots(zp5, nrow = 1) + plot_layout(heights = c(1.6, 1))
  print(fig5)
  # ggsave(fig_out("P4_Fig5_expected_performance.png"), fig5, width = 190, height = 165, units = "mm", dpi = 400)
  # ggsave(fig_out("P4_Fig5_expected_performance.tif"), fig5, width = 190, height = 165, units = "mm", dpi = 600, device = "tiff", compression = "lzw")   # RSE submission (map, 600 dpi)
}

# =============================================================================
# FIGURE 6: national totals by method + survey composition
# =============================================================================
fig("Figure 6: national totals")
if (need(ofile("P4_national_three_methods.csv"))) {
  n3 <- rd(ofile("P4_national_three_methods.csv"))
  line_m <- n3 |> dplyr::filter(method != "Tier 1 (survey-year mean)")
  pt_m   <- n3 |> dplyr::filter(method == "Tier 1 (survey-year mean)")
  col5 <- c("Pixel-based" = "#1b7837", "Tier 1 (pooled mean)" = "#762a83",
            "Tier 1 (survey-year mean)" = "#c2a5cf", "Tier 2 (pooled class means)" = "#e08214")
  p6a <- ggplot(n3, aes(year, total_kt, colour = method, fill = method)) +
    geom_ribbon(data = line_m, aes(ymin = ci_low_kt, ymax = ci_up_kt), alpha = 0.12, colour = NA) +
    geom_line(data = line_m) + geom_point(data = line_m, size = 1) +
    geom_pointrange(data = pt_m, aes(ymin = ci_low_kt, ymax = ci_up_kt), size = 0.25,
                    position = position_nudge(x = 0.12), na.rm = TRUE) +
    annotate("text", x = 2024.3, y = max(n3$ci_up_kt, na.rm = TRUE), label = "2024: no field data",
             size = 2.4, hjust = 1, vjust = 1, colour = "grey30") +
    scale_colour_manual(values = col5, name = NULL) + scale_fill_manual(values = col5, name = NULL) +
    scale_x_continuous(breaks = FIG$years) +
    labs(title = "a  National seagrass AGC by method (95 % CI)", x = NULL, y = "AGC (kt C)") +
    theme_p4() + theme(legend.position = "right")
  comp <- agc_s |> dplyr::count(year, ecoregion)
  top  <- comp |> dplyr::group_by(ecoregion) |> dplyr::summarise(n = sum(n)) |> dplyr::arrange(dplyr::desc(n)) |>
    dplyr::slice_head(n = 6) |> dplyr::pull(ecoregion)
  comp <- comp |> dplyr::mutate(eco = ifelse(ecoregion %in% top, ecoregion, "Other")) |>
    dplyr::group_by(year, eco) |> dplyr::summarise(n = sum(n), .groups = "drop") |>
    dplyr::mutate(eco = factor(eco, levels = c(top, "Other")))
  p6b <- ggplot(comp, aes(factor(year), n, fill = eco)) + geom_col(width = 0.75) +
    scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "grey70"), name = NULL) +
    labs(title = "b  AGC field samples by survey year and ecoregion", x = NULL, y = "Samples") +
    theme_p4() + theme(legend.position = "right")
  fig6 <- p6a / p6b + plot_layout(heights = c(1.4, 1))
  print(fig6)
  # ggsave(fig_out("P4_Fig6_national_totals.png"), fig6, width = 190, height = 150, units = "mm", dpi = 400)
  # ggsave(fig_out("P4_Fig6_national_totals.pdf"), fig6, width = 190, height = 150, units = "mm", device = cairo_pdf)   # RSE submission (vector)
}

# =============================================================================
# FIGURE 7: ecoregions
# =============================================================================
fig("Figure 7: ecoregions")
if (need(ofile("P4_ecoregion_summary.csv"), ofile("P4_lolo_region_means.csv"))) {
  es <- rd(ofile("P4_ecoregion_summary.csv")) |> dplyr::filter(eco_id > 0)
  es <- es |> dplyr::mutate(eco_lab = sprintf("%s (%.0f %s, n = %d)", ecoregion, area_km2, KM2, n_agc))
  ord <- es |> dplyr::arrange(area_km2) |> dplyr::pull(ecoregion)
  lab_of <- setNames(es$eco_lab, es$ecoregion)
  e6 <- dplyr::bind_rows(
    es |> dplyr::transmute(ecoregion, method = "Pixel-based", kt = px_total_tC / 1e3, lo = NA_real_, hi = NA_real_),
    es |> dplyr::transmute(ecoregion, method = "Tier 1, ecoregion field mean", kt = t1_eco_total_tC / 1e3,
                           lo = (t1_eco_total_tC - t1_eco_ci_half_tC) / 1e3, hi = (t1_eco_total_tC + t1_eco_ci_half_tC) / 1e3),
    es |> dplyr::transmute(ecoregion, method = "Tier 1, national mean", kt = t1_national_factor_tC / 1e3, lo = NA_real_, hi = NA_real_),
    es |> dplyr::transmute(ecoregion, method = "Tier 2", kt = t2_total_tC / 1e3, lo = NA_real_, hi = NA_real_)) |>
    dplyr::mutate(ecoregion = factor(ecoregion, levels = ord))
  col6 <- c("Pixel-based" = "#1b7837", "Tier 1, ecoregion field mean" = "#762a83",
            "Tier 1, national mean" = "#c2a5cf", "Tier 2" = "#e08214")
  p7a <- ggplot(e6, aes(kt, ecoregion, colour = method)) +
    geom_linerange(aes(xmin = lo, xmax = hi), na.rm = TRUE, position = position_dodge(width = 0.6)) +
    geom_point(size = 1.5, na.rm = TRUE, position = position_dodge(width = 0.6)) +
    scale_y_discrete(labels = lab_of) +
    scale_colour_manual(values = col6, name = NULL) +
    labs(title = "a  Mean annual AGC per ecoregion", x = "AGC (kt C)", y = NULL) +
    guides(colour = guide_legend(nrow = 2)) +
    theme_p4() + theme(legend.position = "bottom", plot.title.position = "plot")
  lolo <- rd(ofile("P4_lolo_region_means.csv")) |> dplyr::filter(reported) |>
    dplyr::select(ecoregion, n, field, `Pixel-based` = bias_pct_pixel, `Tier 1` = bias_pct_tier1, `Tier 2` = bias_pct_tier2) |>
    tidyr::pivot_longer(c(`Pixel-based`, `Tier 1`, `Tier 2`), names_to = "method", values_to = "bias")
  ord_l <- lolo |> dplyr::distinct(ecoregion, field) |> dplyr::arrange(field) |> dplyr::pull(ecoregion)
  p7b <- ggplot(lolo, aes(bias, factor(ecoregion, levels = ord_l), colour = method)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(size = 1.6, position = position_dodge(width = 0.5)) +
    scale_colour_manual(values = METHOD_COL, name = NULL) +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 10, base = 10),
                       breaks = c(-50, -20, 0, 20, 50, 100, 300)) +
    labs(title = "b  Leave-one-ecoregion-out bias of the regional mean",
         x = "Predicted minus field mean (% of field mean)", y = NULL) +
    theme_p4() + theme(legend.position = "bottom", plot.title.position = "plot")
  fig7 <- p7a / p7b + plot_layout(heights = c(1.2, 1)) +
    plot_annotation(caption = paste("a: ecoregion field mean (Tier 1) only where n >= 30, with 95 % CI; n = AGC field samples.",
                                    "b: ecoregions ordered by field mean AGC; pseudo-log x axis."),
                    theme = theme(plot.caption = element_text(size = FIG$base_size - 2)))
  print(fig7)
  # ggsave(fig_out("P4_Fig7_ecoregions.png"), fig7, width = 190, height = 180, units = "mm", dpi = 400)
  # ggsave(fig_out("P4_Fig7_ecoregions.pdf"), fig7, width = 190, height = 180, units = "mm", device = cairo_pdf)   # RSE submission (vector)
}

# =============================================================================
# FIGURE 8: 22 ROIs, field vs map mean (in-sample and out-of-fold)
# =============================================================================
fig("Figure 8: ROIs")
if (need(ofile("P4_roi_comparison.csv"), ofile("P4_shrinkage_roi.csv"))) {
  r_in  <- rd(ofile("P4_roi_comparison.csv")) |> dplyr::filter(loc_id %in% FIG$roi22) |>
    dplyr::transmute(loc_id, n = n_year, field = field_mean_gC_m2, map = px_mean_gC_m2,
                     panel = "a  Final map (ROI samples in training)")
  r_oof <- rd(ofile("P4_shrinkage_roi.csv")) |> dplyr::filter(loc_id %in% FIG$roi22) |>
    dplyr::transmute(loc_id, n, field = mean_obs, map = mean_pred_oof,
                     panel = "b  Out-of-fold (meadow cluster held out)")
  cat(sprintf("ROIs found: %d of 22 (final map), %d of 22 (out-of-fold)\n", nrow(r_in), nrow(r_oof)))
  r7 <- dplyr::bind_rows(r_in, r_oof) |> dplyr::mutate(n_class = ifelse(n >= 10, "n >= 10", "n < 10"))
  st7 <- r7 |> dplyr::filter(n >= 10) |> dplyr::group_by(panel) |>
    dplyr::summarise(slope = coef(lm(map ~ field))[2], r2 = cor(field, map)^2,
                     mec = 1 - sum((field - map)^2) / sum((field - mean(field))^2),
                     rmse = sqrt(mean((field - map)^2)), k = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(txt = sprintf("n >= 10 (%d ROIs): slope %.2f, R\u00b2 %.2f, MEC %.2f, RMSE %.1f", k, slope, r2, mec, rmse))
  lim <- c(0, max(c(r7$field, r7$map), na.rm = TRUE) * 1.05)
  fig8 <- ggplot(r7, aes(field, map)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
    geom_hline(yintercept = nat_mean, linetype = 3, colour = "firebrick") +
    geom_smooth(data = dplyr::filter(r7, n >= 10), method = "lm", formula = y ~ x, se = FALSE,
                colour = "#1b7837", linewidth = 0.5) +
    geom_point(aes(size = n, shape = n_class), colour = "#1b7837", alpha = 0.8) +
    ggrepel::geom_text_repel(aes(label = sub("ROI-", "", loc_id)), size = 2, colour = "grey30",
                             max.overlaps = 30, min.segment.length = 0.2) +
    geom_text(data = st7, aes(x = lim[1], y = lim[2], label = txt), hjust = 0, vjust = 1, size = 2.4,
              inherit.aes = FALSE) +
    facet_wrap(~ panel) +
    scale_shape_manual(values = c("n >= 10" = 16, "n < 10" = 1), name = NULL) +
    scale_size_continuous(range = c(1, 4), name = "Samples") +
    coord_equal(xlim = lim, ylim = lim) +
    labs(x = expression("Field mean AGC (gC m"^-2*")"), y = expression("Map mean AGC (gC m"^-2*")"),
         caption = "Dashed: 1:1. Dotted: national field mean (Tier 1 factor).") +
    theme_p4() + theme(legend.position = "right")
  print(fig8)
  # ggsave(fig_out("P4_Fig8_roi.png"), fig8, width = 190, height = 95, units = "mm", dpi = 400)
  # ggsave(fig_out("P4_Fig8_roi.pdf"), fig8, width = 190, height = 95, units = "mm", device = cairo_pdf)   # RSE submission (vector)
}

# =============================================================================
# FIGURE 9: temporal behaviour
# =============================================================================
fig("Figure 9: temporal")
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
    labs(title = "b  GSE similarity", x = "Consecutive years",
         y = "Dot product (1 = identical)") +
    theme_p4() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
ann <- win |> dplyr::filter(file.exists(zoom_file(stem, "agc_annual")))
if (nrow(ann) > 0) {
  # One line per survey cohort (samples of one survey year in one window, n >= 20):
  # mean mapped AGC at the cohort's pixels in every map year, and the cohort's field mean
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
  # Numbers for the text (9 Nov 2026): per window, all samples with a mapped value. Survey year:
  # mean map value in each sample's own survey year; other years: mean over the other map years
  wsum <- purrr::map_dfr(seq_len(nrow(ann)), function(i) {
    st <- ann$stem[i]; A <- terra::rast(zoom_file(st, "agc_annual")); e <- as.vector(terra::ext(A))
    sw <- agc_s |> dplyr::filter(dplyr::between(xcoord, e[1], e[2]), dplyr::between(ycoord, e[3], e[4]))
    if (nrow(sw) == 0) return(NULL)
    v <- as.matrix(terra::extract(A, as.matrix(sw[, c("xcoord", "ycoord")])))
    ok <- !is.na(v[, 1]); sw <- sw[ok, ]; v <- v[ok, , drop = FALSE]
    jy <- match(sw$year, FIG$years)
    own <- v[cbind(seq_len(nrow(v)), jy)]
    oth <- vapply(seq_len(nrow(v)), function(k) mean(v[k, -jy[k]]), numeric(1))
    tibble::tibble(window = ann$label[i], stem = st, n = nrow(sw), survey_years = paste(sort(unique(sw$year)), collapse = ", "),
                   field_mean_gC_m2 = mean(sw$AGC_pred), map_survey_year_gC_m2 = mean(own),
                   map_other_years_gC_m2 = mean(oth), map_all_years_gC_m2 = mean(v))
  })
  if (nrow(wsum) > 0) {
    print(as.data.frame(wsum), row.names = FALSE, digits = 4)
    readr::write_csv(wsum, fig_out("P4_Fig9c_window_summary.csv"))
  }
  if (nrow(yr9) > 0) {
    readr::write_csv(yr9 |> dplyr::mutate(cohort = as.character(cohort)), fig_out("P4_Fig9c_window_years.csv"))
    # cohorts in year order; the sample count of each cohort goes into the panel strip
    yr9 <- yr9 |> dplyr::mutate(cohort = factor(as.character(cohort), levels = sort(unique(as.character(cohort)))))
    lab9 <- yr9 |> dplyr::distinct(label, cohort, n) |> dplyr::arrange(label, cohort) |> dplyr::group_by(label) |>
      dplyr::summarise(txt = sprintf("%s (%s)", dplyr::first(label),
                                     paste(sprintf("%s: n = %d", cohort, n), collapse = "; ")), .groups = "drop")
    yr9 <- yr9 |> dplyr::mutate(label = factor(label, levels = lab9$label, labels = lab9$txt))
    p9$c <- ggplot(yr9, aes(year, map, colour = cohort)) +
      geom_hline(aes(yintercept = field, colour = cohort), linetype = 2, linewidth = 0.4) +
      geom_hline(yintercept = nat_mean, linetype = 3, colour = "firebrick") +
      geom_line() +
      geom_point(data = dplyr::filter(yr9, !survey), shape = 1, size = 1.5) +
      geom_point(data = dplyr::filter(yr9, survey), shape = 16, size = 2.4) +
      facet_wrap(~ label, nrow = 1) +
      scale_colour_manual(values = YEAR_COL, name = "Survey year of the samples") +
      scale_x_continuous(breaks = FIG$years, labels = substr(FIG$years, 3, 4)) +
      labs(title = "c  Mapped AGC at the sample pixels, by map year",
           caption = "Filled point: survey year of the samples; dashed: their field mean; dotted: national field mean.",
           x = "Map year", y = expression("AGC (gC m"^-2*")")) +
      theme_p4() + theme(legend.position = "bottom")
  }
}
if (length(p9) > 0) {
  fig9 <- if (all(c("a", "b", "c") %in% names(p9))) (p9$a + p9$b + plot_layout(widths = c(2.6, 1))) / p9$c +
    plot_layout(heights = c(0.85, 1)) else wrap_plots(p9, ncol = 1)
  print(fig9)
  # ggsave(fig_out("P4_Fig9_temporal.png"), fig9, width = 190, height = 140, units = "mm", dpi = 400)
  # ggsave(fig_out("P4_Fig9_temporal.tif"), fig9, width = 190, height = 140, units = "mm", dpi = 600, device = "tiff", compression = "lzw")   # RSE submission (map, 600 dpi)
}

# =============================================================================
# SUPPLEMENT S1: presence probability threshold per year
# =============================================================================
fig("Figure S1: PA threshold by year")
if (need(ofile("P4_pa_threshold_by_year.csv"))) {
  pa <- rd(ofile("P4_pa_threshold_by_year.csv")) |> dplyr::mutate(year = as.character(year)) |>
    dplyr::filter(!is.na(best_F1)) |>
    tidyr::pivot_longer(c(best_F1, F1_050, F1_070), names_to = "what", values_to = "F1") |>
    dplyr::mutate(what = dplyr::recode(what, best_F1 = "best threshold", F1_050 = "threshold 0.5", F1_070 = "threshold 0.7"))
  figS1 <- ggplot(pa, aes(year, F1, colour = what, group = what)) + geom_point(size = 1.8) + geom_line() +
    scale_colour_manual(values = c("best threshold" = "grey30", "threshold 0.5" = "#2171b5", "threshold 0.7" = "#d73027"),
                        name = NULL) +
    labs(x = NULL, y = "F1 (PA points)", caption = "Only 2017 is out-of-sample for the PA model; other years trained it.") +
    theme_p4() + theme(legend.position = "bottom")
  print(figS1)
  # ggsave(fig_out("P4_FigS1_pa_threshold.png"), figS1, width = 90, height = 70, units = "mm", dpi = 400)
  # ggsave(fig_out("P4_FigS1_pa_threshold.pdf"), figS1, width = 90, height = 70, units = "mm", device = cairo_pdf)   # RSE submission (vector)
}

# =============================================================================
# SUPPLEMENT S2: extent components and mask diagnostic
# =============================================================================
fig("Figure S2: extent components and mask diagnostic")
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
  # share of the shallow area per mask reason (legend), and a crop to the shallow pixels
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
  # the Rote panel is a narrow strip: keep it and its legend in the left part of the lower row
  figS2 <- if (length(s2p) == 2) s2p$a / (s2p$b + plot_spacer() + plot_layout(widths = c(1.3, 1))) +
    plot_layout(heights = c(1, 1.3)) else wrap_plots(s2p, ncol = 1)
  print(figS2)
  # ggsave(fig_out("P4_FigS2_extent_masks.png"), figS2, width = 190, height = 170, units = "mm", dpi = 400)
  # ggsave(fig_out("P4_FigS2_extent_masks.tif"), figS2, width = 190, height = 170, units = "mm", dpi = 600, device = "tiff", compression = "lzw")   # RSE submission (map, 600 dpi)
}

# =============================================================================
# SUPPLEMENT S3: error by observed AGC and shrinkage of group means
# =============================================================================
fig("Figure S3: error by AGC bin and shrinkage")
if (need(ofile("P4_error_by_agc_bin.csv"), ofile("P4_shrinkage_groups.csv"))) {
  eb <- rd(ofile("P4_error_by_agc_bin.csv")) |>
    dplyr::mutate(lo = as.numeric(sub("^\\[(\\d+),.*", "\\1", bin)), design = factor(design, levels = DESIGN_LEVELS))
  pS3a <- ggplot(eb, aes(lo + 5, mean_pred, colour = design)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
    geom_hline(yintercept = nat_mean, linetype = 3, colour = "firebrick") +
    geom_line() + geom_point(aes(size = n_per_repeat)) +
    scale_colour_manual(values = DESIGN_COL[c("pixel", "cluster")], labels = DESIGN_LABEL, name = NULL) +
    scale_size_continuous(range = c(0.6, 3), name = "Samples per bin") +
    labs(title = "a  Mean out-of-fold prediction per bin of observed AGC (bin width 10)",
         x = expression("Observed AGC, bin centre (gC m"^-2*")"),
         y = expression("Mean prediction (gC m"^-2*")")) + theme_p4()
  sg <- rd(ofile("P4_shrinkage_groups.csv"))
  pS3b <- ggplot(sg, aes(mean_obs, mean_pred)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
    geom_hline(yintercept = nat_mean, linetype = 3, colour = "firebrick") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "steelblue", linewidth = 0.5) +
    geom_point(aes(size = n), alpha = 0.6) + facet_wrap(~ level) +
    scale_size_continuous(range = c(0.6, 3.5), name = "Samples") +
    labs(title = "b  Group means: field vs out-of-fold pixel-based", x = expression("Field mean (gC m"^-2*")"),
         y = expression("Predicted mean (gC m"^-2*")")) + theme_p4()
  figS3 <- pS3a / pS3b
  if (file.exists(ofile("P4_cover_error_by_bin.csv"))) {
    # c (9 Nov 2026): mean error of the cover model per 10 % bin of observed cover (P4_04 STEP 6)
    cb <- rd(ofile("P4_cover_error_by_bin.csv")) |>
      dplyr::mutate(lo = as.numeric(sub("^\\[(\\d+),.*", "\\1", bin)), design = factor(design, levels = DESIGN_LEVELS))
    pS3c <- ggplot(cb, aes(lo + 5, mean_err, colour = design)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
      geom_line() + geom_point(aes(size = n_per_repeat)) +
      scale_colour_manual(values = DESIGN_COL[c("pixel", "cluster")], labels = DESIGN_LABEL, name = NULL) +
      scale_size_continuous(range = c(0.6, 3), name = "Samples per bin") +
      scale_x_continuous(breaks = seq(0, 100, 20)) +
      labs(title = "c  Mean out-of-fold error of the cover model per bin of observed cover (bin width 10)",
           x = "Observed cover, bin centre (%)", y = "Mean error, prediction - observation (% cover)") + theme_p4()
    figS3 <- figS3 / pS3c
  } else cat("  -> panel c skipped, missing: P4_cover_error_by_bin.csv\n")
  print(figS3)
  # ggsave(fig_out("P4_FigS3_error_shrinkage.png"), figS3, width = 190, height = 240, units = "mm", dpi = 400)
  # ggsave(fig_out("P4_FigS3_error_shrinkage.pdf"), figS3, width = 190, height = 240, units = "mm", device = cairo_pdf)   # RSE submission (vector)
}

# =============================================================================
# SUPPLEMENT S4: extent under the geomorphic mask options vs KKP 2025
# =============================================================================
fig("Figure S4: extent and the geomorphic mask")
if (need(ofile("P4_extent_sensitivity_ecoregion.csv"), ofile("P4_pixel_national.csv"))) {
  ee <- rd(ofile("P4_extent_sensitivity_ecoregion.csv"))
  final_km2 <- mean(rd(ofile("P4_pixel_national.csv"))$area_km2)
  core_km2  <- if (file.exists(grid_file("agcmean"))) sum(rd(grid_file("agcmean"))$area_core_m2) / 1e6 else NA_real_
  f30 <- core_km2 / sum(ee$current_km2)                       # 30 m -> 10 m scale (P4_07 STEP 12)
  reg_km2 <- f30 * (sum(ee$only16_km2) + sum((ee$nogeo_km2 - ee$only16_km2)[!ee$ecoregion %in% FIG$turbid_eco]))
  s7 <- tibble::tibble(
    option = c("Core extent only", sprintf("Final: %s removed nationally", FIG$trf_label),
               sprintf("%s removed only in %s", FIG$trf_label, paste(FIG$turbid_eco, collapse = ", ")),
               "No geomorphic mask"),
    km2 = c(core_km2, final_km2, reg_km2, f30 * sum(ee$nogeo_km2)),
    kind = c("mapped", "mapped", "estimate", "estimate")) |>
    dplyr::mutate(option = factor(option, levels = rev(option)), pct_kkp = 100 * km2 / FIG$kkp_km2)
  print(as.data.frame(s7), row.names = FALSE, digits = 4)
  figS4 <- ggplot(s7, aes(km2, option, fill = kind)) + geom_col(width = 0.6) +
    geom_vline(xintercept = FIG$kkp_km2, linetype = 2, colour = "firebrick") +
    annotate("text", x = FIG$kkp_km2, y = 4.75, label = sprintf(" KKP 2025: %.0f %s", FIG$kkp_km2, KM2),
             colour = "firebrick", size = 2.5, hjust = 0, vjust = 0.5) +
    geom_text(aes(label = sprintf("%.0f %s (%.0f%% of KKP)", km2, KM2, pct_kkp)), hjust = -0.05, size = 2.4) +
    scale_fill_manual(values = c(mapped = "#2171b5", estimate = "#9ecae1"),
                      labels = c(mapped = "mapped (10 m)", estimate = "estimate (30 m sensitivity run, P4_G05, scaled to 10 m)"), name = NULL) +
    scale_x_continuous(limits = c(0, max(s7$km2) * 1.45), expand = expansion(mult = c(0, 0.02))) +
    scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
    labs(x = expression("Extent area (km"^2*")"), y = NULL) + theme_p4() +
    theme(legend.position = "bottom")
  print(figS4)
  # ggsave(fig_out("P4_FigS4_extent_masks_kkp.png"), figS4, width = 170, height = 65, units = "mm", dpi = 400)
  # ggsave(fig_out("P4_FigS4_extent_masks_kkp.pdf"), figS4, width = 170, height = 65, units = "mm", device = cairo_pdf)   # RSE submission (vector)
}

# =============================================================================
# SUPPLEMENT S5: annual AGC maps
# =============================================================================
fig("Figure S5: annual AGC maps")
yf <- gfile(sprintf("P4_grid%dm_agc_%d.csv", FIG$grid_m, FIG$years))
if (need(yf)) {
  s1 <- purrr::map_dfr(FIG$years, function(y) {
    agg_cells(rd(gfile(sprintf("P4_grid%dm_agc_%d.csv", FIG$grid_m, y))), FIG$ext_cell, c("area_m2", "agc_g")) |>
      dplyr::mutate(year = y, agc = agc_g / area_m2)
  })
  tot <- s1 |> dplyr::group_by(year) |> dplyr::summarise(kt = sum(agc_g) / 1e9, .groups = "drop")
  figS5 <- ggplot() + land() +
    geom_tile(data = s1, aes(gx, gy, fill = agc), width = FIG$ext_cell, height = FIG$ext_cell) +
    geom_text(data = tot, aes(x = 140.5, y = 5.5, label = sprintf("%.0f kt C", kt)), hjust = 1, size = 2.3) +
    facet_wrap(~ year, ncol = 2) +
    agc_scale(name = expression("AGC (gC m"^-2*")")) +
    guides(fill = guide_colourbar(barwidth = unit(60, "mm"), title.position = "top", title.hjust = 0.5)) +
    nat_coord() + theme_map() + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                                      legend.position = "bottom")
  print(figS5)
  # ggsave(fig_out("P4_FigS5_annual_maps.png"), figS5, width = 190, height = 175, units = "mm", dpi = 300)
  # ggsave(fig_out("P4_FigS5_annual_maps.tif"), figS5, width = 190, height = 175, units = "mm", dpi = 600, device = "tiff", compression = "lzw")   # RSE submission (map, 600 dpi)
}

# =============================================================================
# SUPPLEMENT S6: cover class confusion (out-of-fold)
# =============================================================================
fig("Figure S6: cover class confusion")
if (need(ofile("P4_cover_confusion.csv"))) {
  cf <- rd(ofile("P4_cover_confusion.csv")) |>
    dplyr::mutate(design = factor(design, levels = DESIGN_LEVELS, labels = DESIGN_LABEL[DESIGN_LEVELS]),
                  obs_class = factor(obs_class, levels = CFG$cover_labels),
                  pred_class = factor(pred_class, levels = CFG$cover_labels)) |>
    tidyr::complete(design, obs_class, pred_class, fill = list(n = 0)) |>
    dplyr::group_by(design, obs_class) |> dplyr::mutate(row_pct = 100 * n / sum(n)) |> dplyr::ungroup()
  figS6 <- ggplot(cf, aes(pred_class, obs_class, fill = row_pct)) +
    geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.0f", row_pct)), size = 2.4) +
    facet_wrap(~ design, nrow = 1) +
    scale_fill_gradient(low = "white", high = "#2171b5", limits = c(0, 100), name = "% of observed class",
                        guide = guide_colourbar(barwidth = unit(40, "mm"))) +
    labs(x = "Predicted cover class (%)", y = "Observed cover class (%)") +
    theme_p4() + theme(legend.position = "bottom")
  print(figS6)
  # ggsave(fig_out("P4_FigS6_cover_confusion.png"), figS6, width = 190, height = 70, units = "mm", dpi = 400)
  # ggsave(fig_out("P4_FigS6_cover_confusion.pdf"), figS6, width = 190, height = 70, units = "mm", device = cairo_pdf)   # RSE submission (vector)
}

# =============================================================================
# SUPPLEMENT S7: cumulative extent by distance to AGC samples
# =============================================================================
fig("Figure S7: cumulative extent by distance")
if (need(ofile("P4_distance_coverage_class.csv"))) {
  dc <- rd(ofile("P4_distance_coverage_class.csv")) |>
    dplyr::mutate(dist_label = factor(dist_label, levels = dist_label), design = factor(design, levels = DESIGN_LEVELS))
  figS7 <- ggplot(dc, aes(dist_label, area_cum_pct, group = 1)) +
    geom_col(aes(y = area_pct, fill = design), width = 0.7) + geom_line() + geom_point(size = 1) +
    scale_fill_manual(values = DESIGN_COL, labels = DESIGN_LABEL, name = "Validation design") +
    labs(x = "Distance to the nearest AGC sample", y = "Share of extent (%)",
         caption = "Bars: share per distance class; line: cumulative share") +
    theme_p4() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  print(figS7)
  # ggsave(fig_out("P4_FigS7_distance_cumulative.png"), figS7, width = 160, height = 80, units = "mm", dpi = 400)
  # ggsave(fig_out("P4_FigS7_distance_cumulative.pdf"), figS7, width = 160, height = 80, units = "mm", device = cairo_pdf)   # RSE submission (vector)
}

# =============================================================================
# SUPPLEMENT S8: cover classes (Tier 2 strata) in the zoom windows
# =============================================================================
fig("Figure S8: cover classes in the zoom windows")
if (nrow(win) > 0) {
  zp8 <- purrr::map(seq_len(nrow(win)), function(i) {
    st <- win$stem[i]
    d <- rast_df(zoom_file(st, "layers"), "cover_class") |>
      dplyr::mutate(cls = factor(CFG$cover_labels[value], levels = CFG$cover_labels))
    ggplot() + s2_layer(st) + geom_raster(data = d, aes(x, y, fill = cls)) +
      scale_fill_manual(values = c("<30" = "#fee391", "30-59.9" = "#41ab5d", ">=60" = "#00441b"),
                        drop = FALSE, name = "Mapped cover (%)") +
      zoom_coord(st) + zoom_bar(st) + labs(title = sprintf("%s  %s", letters[i], win$label[i])) + zoom_theme()
  })
  figS8 <- wrap_plots(zp8, nrow = 1, guides = "collect") & theme(legend.position = "bottom")
  print(figS8)
  # ggsave(fig_out("P4_FigS8_cover_classes.png"), figS8, width = 190, height = 70, units = "mm", dpi = 400)
  # ggsave(fig_out("P4_FigS8_cover_classes.tif"), figS8, width = 190, height = 70, units = "mm", dpi = 600, device = "tiff", compression = "lzw")   # RSE submission (map, 600 dpi)
}

cat("\nP4_08_figures.R COMPLETE\n")
