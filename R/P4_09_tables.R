# =============================================================================
# P4_09_tables.R
# Paper 4, design v3 (15 Nov 2026): all tables of the manuscript and the
# Supplement, built from the R results (P4_01, P4_02, P4_04, P4_05, P4_07) and
# the GEE exports. Nothing is recomputed from raw data: every number comes from a
# result CSV, so a table and its source can be checked line by line.
#
# Main tables (numbers in TABNUM below; renumber there when the text is written)
#   Table 1  Field data by ecoregion: extent, samples, sites and meadows tested, Tier 2 stratum
#   Table 2  AGC totals by method, by ecoregion and for Indonesia (aim 1)
#   Table 3  Inventory-level validation: sites, meadows and ecoregions held out (aim 2)
# Supplement (provisional order of first citation)
#   S1  AGC samples by ecoregion and survey year
#   S2  Presence model cross-validation (P4_G09)
#   S3  AGC samples inside the extent and the masks that removed the rest
#   S4  Extent under the geomorphic mask options and the KKP 2025 area
#   S5  Sentinel-2 L2A availability at the AGC samples, by survey year and ecoregion
#   S6  Random forest settings of the four AGC models
#   S7  National AGC by year and method (95 % CI)
#   S8  Tier 1 by survey year and the survey composition
#   S9  Tier 2 strata: ecoregion means, effective n, totals and leave-one-year-out range
#   S10 Clustering, uncertainty parameters and checks of the totals
#   S11 Site-level accuracy by ecoregion
#   S12 Accuracy for each ecoregion when it is held out
#   S13 Sample-level accuracy by validation design, method and sample set
#   S14 Distance to AGC samples and between-year CV by ecoregion
# Dropped from P4v2: the 22 ROIs (S4), the cover model (S8 classes, S9) and the
# cover columns of Table 1.
#
# Output (new folder, P4_tables_v2 is not touched): <out_dir>/P4_tables_v3/
#   P4_T1_<key>.csv ... P4_S14_<key>.csv  one file per table (UTF-8 with BOM, opens in Excel)
#   P4_tables_all.md                all tables with their notes, in manuscript order
#   P4_09_checks.txt                consistency checks between the source files
#   P4_09_derived_numbers.csv       numbers computed here and quoted in the text
#   inputs/                         copies of the source CSVs (zip P4_tables_v3 and send)
# A missing input only skips the part that needs it. All export lines are active.
#
# Number format: values with a 95 % range are written "mean [low, high]". For the
# validation tables the range is the 2.5-97.5 percentile over the CV repeats
# (10 repeats; 3 seeds for the ecoregion design); for totals it is the 95 % CI.
# =============================================================================

source("R/P4_00_config_functions.R")
cat("P4_09_tables.R | version 2026-11-16 (v3)\n")

# -----------------------------------------------------------------------------
# TABLE CONFIG -- edit here only
# -----------------------------------------------------------------------------
TAB <- list(
  dir        = out_file("P4_tables_v3"),
  gee_dir    = if (!is.null(CFG$gee_dir)) CFG$gee_dir else out_file("gee"),
  trf_label  = "Terrestrial Reef Flat",       # ACA geomorphic class removed from the extent
  kkp_km2    = 6601.56,                     # KKP 2025 national seagrass map (official area)
  turbid_eco = c("Sunda Shelf/Java Sea"),   # S4: option with Terrestrial Reef Flat removed only here
  years      = 2017:2024,
  tol_pct    = 0.1,                         # consistency checks: allowed difference (%)
  # Sentinel-2 scenes over the AGC samples per survey year, printed by P4_G11 in the GEE
  # console (run of 14 Nov 2026). Earth Engine does not hold the full L2A archive for
  # 2017-2018, so the Sentinel-2 comparison uses the 2019-2023 surveys (CFG$s2_years).
  s2_scenes = tibble::tribble(
    ~year, ~l2a_scenes, ~l1c_scenes,
    2017,          9,         143,
    2018,         89,         839,
    2019,       1234,        1238,
    2020,        322,         324,
    2021,       3556,        3561,
    2022,        805,         805,
    2023,       1015,        1008)
)
# Table numbers by key: change here when the order of first citation is known
TABNUM <- c(data = "1", totals = "2", inventory = "3",
            samples_year = "S1", pa_cv = "S2", in_extent = "S3", extent_kkp = "S4", s2_avail = "S5",
            rf_settings = "S6", national_year = "S7", tier1_year = "S8", tier2_strata = "S9",
            uncertainty = "S10", site_eco = "S11", lolo = "S12", sample_level = "S13", dist_cv = "S14")
# 16 Nov 2026: order of first citation in the manuscript (site_eco, lolo and
# sample_level moved: Results 3.3 cites the site and ecoregion tables before the sample level)

DESIGN_LEVELS <- c("pixel", "b100", "b1000", "cluster", "ecoregion")
DESIGN_NAME <- c(pixel = "Pixel", b100 = "100 m block", b1000 = "1 km block",
                 cluster = "Meadow cluster", ecoregion = "Ecoregion")
DESIGN_RANGE <- c(pixel = "≤ 10 m", b100 = "10–100 m", b1000 = "100 m–1 km",
                  cluster = "1–20 km", ecoregion = "> 20 km")
METHODS <- c("Tier 1", "Tier 2", "Local mean", "Pixel-based GSE", "Pixel-based S2", "Pixel-based S2 (median)")
SET_NAME <- c(all = "All AGC samples", s2 = "Sentinel-2 set")
LEVEL_NAME <- c(site = "Site (100 m block)", meadow = "Meadow (DBSCAN cluster)", ecoregion = "Ecoregion")
MODEL_NAME <- c(AGC = "GSE, all AGC samples", AGC_GSE_SUB = "GSE, Sentinel-2 set",
                AGC_S2 = "Sentinel-2 percentiles, Sentinel-2 set", AGC_S2MED = "Sentinel-2 median, Sentinel-2 set")
KM2 <- "km²"; GCM2 <- "gC m⁻²"; R2 <- "R²"

dir.create(TAB$dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(TAB$dir, "inputs"), showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
used <- character(0)
src <- function(name) {                       # read a source CSV and remember it
  p <- out_file(name)
  if (!file.exists(p)) { cat("  -> missing:", name, "\n"); return(NULL) }
  used <<- union(used, name)
  readr::read_csv(p, show_col_types = FALSE, guess_max = 100000)
}
have <- function(...) all(file.exists(out_file(c(...))))
num <- function(x, d = 1, big = ",") {          # "+ 0" turns a rounded -0 into 0
  ifelse(is.na(x), "–", formatC(round(x, d) + 0, format = "f", digits = d, big.mark = big))
}
share_lab <- function(p) ifelse(!is.na(p) & p > 0 & p < 0.05, "< 0.1", num(p, 1))
rng <- function(m, lo, hi, d = 1) ifelse(is.na(m), "–",
                                         ifelse(is.na(lo) | is.na(hi), num(m, d),
                                                sprintf("%s [%s, %s]", num(m, d), num(lo, d), num(hi, d))))
pm  <- function(m, s, d = 1) ifelse(is.na(m), "–", sprintf("%s ± %s", num(m, d), num(s, d)))
span <- function(lo, hi, d = 1) ifelse(is.na(lo), "–", sprintf("%s–%s", num(lo, d), num(hi, d)))
dist_lab <- function(m) ifelse(m < 1000, sprintf("%.0f m", m),
                               ifelse(m < 1e5, sprintf("%.1f km", m / 1000), sprintf("%.0f km", m / 1000)))
in_text <- function(m) ifelse(m == "Local mean", "local mean", m)
# Methods that are identical by definition in a design (local mean = Tier 2 when the
# meadow is held out; Tier 2 = Tier 1 when the ecoregion is held out) are merged into
# one row "Tier 1 (= Tier 2 = local mean)". d must be ordered by method.
collapse_same <- function(d, group_cols, value_cols) {
  sig <- do.call(paste, c(lapply(value_cols, function(v) sprintf("%.9g", d[[v]])), sep = "|"))
  d$.sig <- sig
  d |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, ".sig")))) |>
    dplyr::mutate(method_label = if (dplyr::n() > 1)
      sprintf("%s (= %s)", dplyr::first(method), paste(in_text(method[-1]), collapse = " = ")) else method,
      .first = dplyr::row_number() == 1) |>
    dplyr::ungroup() |> dplyr::filter(.first) |> dplyr::select(-.sig, -.first)
}

TABLES <- list()                               # id -> list(title, data, notes, file)
CHECKS <- character(0)
check <- function(label, ok, detail = "") {
  line <- sprintf("[%s] %s%s", if (isTRUE(ok)) "OK  " else "CEK ", label,
                  if (nzchar(detail)) paste0(": ", detail) else "")
  cat(line, "\n"); CHECKS <<- c(CHECKS, line)
}
close_pct <- function(a, b) abs(a - b) <= TAB$tol_pct / 100 * max(abs(b), 1e-9)
# Numbers computed in this script and quoted in the text, written to
# P4_09_derived_numbers.csv
DERIVED <- tibble::tibble(name = character(0), value = numeric(0), how = character(0))
derive <- function(name, value, how) {
  DERIVED <<- dplyr::bind_rows(DERIVED, tibble::tibble(name = name, value = as.numeric(value), how = how))
}

save_table <- function(key, title, d, notes = character(0)) {
  id <- TABNUM[[key]]
  file <- sprintf("P4_%s_%s.csv", if (substr(id, 1, 1) == "S") id else paste0("T", id), key)
  title <- sprintf("Table %s. %s", id, title)
  readr::write_excel_csv(d, file.path(TAB$dir, file), na = "")
  TABLES[[id]] <<- list(title = title, data = d, notes = notes[nzchar(notes)], file = file)
  cat(sprintf("\n==== %s ====\n%s\n", file, title))
  print(as.data.frame(d), row.names = FALSE, right = FALSE)
  if (length(notes) > 0) cat(paste0("  ", notes, collapse = "\n"), "\n")
}
step <- function(txt) cat(sprintf("\n######## %s ########\n", txt))

# -----------------------------------------------------------------------------
# Inputs used by several tables
# -----------------------------------------------------------------------------
samples <- read_samples(); used <- union(used, "P4_samples.csv")
agc <- samples |> dplyr::filter(is_agc)
px_nat <- src("P4_pixel_national.csv")
nat3   <- src("P4_national_three_methods.csv")
eco    <- src("P4_ecoregion_summary.csv")
iu     <- src("P4_inventory_units.csv")
im     <- src("P4_inventory_metrics.csv")
ms     <- src("P4_metrics_summary.csv")
if (!is.null(ms) && !"sample_set" %in% names(ms)) ms$sample_set <- "all"
bp     <- src("P4_best_params.csv")
t2n    <- src("P4_tier2_national.csv")
t2e    <- src("P4_tier2_by_ecoregion.csv")
t1b    <- src("P4_tier1_agc_bounds_gee.csv")
t1g    <- src("P4_tier1_national_gee.csv")
awb    <- src("P4_lolo_area_weighted_bias.csv")
dcd    <- src("P4_distance_coverage_design.csv")
AREA <- if (!is.null(px_nat)) mean(px_nat$area_km2) else NA_real_
NAT_MEAN <- mean(agc$AGC_pred)
N_S2 <- if ("has_s2" %in% names(agc)) sum(agc$has_s2 %in% TRUE) else NA_integer_
cat(sprintf("AGC samples %d (Sentinel-2 set %s), extent %.2f %s, national field mean %.3f %s\n",
            nrow(agc), num(N_S2, 0), AREA, KM2, NAT_MEAN, GCM2))
derive("n_agc", nrow(agc), "AGC samples (P4_samples.csv, is_agc)")
derive("n_s2_set", N_S2, "AGC samples with has_s2 (survey years CFG$s2_years, clear Sentinel-2 pixel)")
derive("pct_s2_set", 100 * N_S2 / nrow(agc), "Sentinel-2 set as % of the AGC samples")

# =============================================================================
# CONSISTENCY CHECKS between the source files
# =============================================================================
step("CHECKS")
if (!is.null(eco) && !is.null(px_nat)) {
  check("Ecoregion areas add up to the national extent",
        close_pct(sum(eco$area_km2), AREA), sprintf("%.2f vs %.2f %s", sum(eco$area_km2), AREA, KM2))
  check("Ecoregion pixel-based totals add up to the national mean total",
        close_pct(sum(eco$px_total_tC), mean(px_nat$total_tC)),
        sprintf("%.0f vs %.0f t C", sum(eco$px_total_tC), mean(px_nat$total_tC)))
  check("AGC samples in the ecoregion table = AGC samples", sum(eco$n_agc, na.rm = TRUE) == nrow(agc),
        sprintf("%d vs %d", sum(eco$n_agc, na.rm = TRUE), nrow(agc)))
}
eby <- src("P4_ecoregion_by_year.csv")
if (!is.null(eby) && !is.null(px_nat)) {
  s <- eby |> dplyr::group_by(year) |> dplyr::summarise(t = sum(px_total_tC), a = sum(area_km2), .groups = "drop") |>
    dplyr::inner_join(px_nat |> dplyr::select(year, total_tC, area_km2), by = "year")
  check("Ecoregion-by-year totals add up to the national totals every year",
        all(close_pct(s$t, s$total_tC)) && nrow(s) == nrow(px_nat),
        sprintf("largest difference %.3f %%", max(abs(s$t - s$total_tC) / s$total_tC * 100)))
}
if (!is.null(nat3)) {
  t1p <- nat3 |> dplyr::filter(method == "Tier 1 (pooled mean)")
  check("Tier 1 pooled total = national field mean x extent",
        close_pct(t1p$total_tC[1], NAT_MEAN * AREA), sprintf("%.0f vs %.0f t C", t1p$total_tC[1], NAT_MEAN * AREA))
  pxn <- nat3 |> dplyr::filter(method == "Pixel-based GSE")
  check("Pixel-based totals match P4_pixel_national", all(close_pct(sort(pxn$total_tC), sort(px_nat$total_tC))))
  if (!is.null(t1b)) {
    b <- t1b |> dplyr::filter(response == "AGC_pred")
    check("Tier 1 from P4_tier1_agc_bounds_gee (AGC_pred) = Tier 1 pooled total",
          close_pct(b$total_tC, t1p$total_tC[1]), sprintf("%.0f vs %.0f t C", b$total_tC, t1p$total_tC[1]))
  }
}
if (!is.null(t2n) && !is.null(nat3)) {
  t2a <- nat3 |> dplyr::filter(method == "Tier 2 (ecoregion strata)") |>
    dplyr::inner_join(t2n |> dplyr::select(year, t2 = total_tC), by = "year")
  check("Tier 2 totals in P4_tier2_national match P4_national_three_methods",
        nrow(t2a) > 0 && all(close_pct(t2a$t2, t2a$total_tC)), sprintf("%d years", nrow(t2a)))
}
if (!is.null(t2n) && !is.null(t2e))
  check("Tier 2 national total = sum of the ecoregion terms (P4_tier2_by_ecoregion)",
        close_pct(sum(t2e$t2_total_tC), t2n$total_tC[1]), sprintf("%.0f vs %.0f t C", sum(t2e$t2_total_tC), t2n$total_tC[1]))
if (!is.null(dcd)) check("Distance shares by design add up to 100 %", close_pct(sum(dcd$area_pct), 100),
                         sprintf("%.3f %%", sum(dcd$area_pct)))
if (!is.null(t1g) && "pooled" %in% t1g$year) {
  pooled <- t1g |> dplyr::filter(year == "pooled")
  check("Tier 1 pooled mean in P4_tier1_national_gee = national field mean",
        close_pct(pooled$mean_gC_m2, NAT_MEAN), sprintf("%.3f vs %.3f", pooled$mean_gC_m2, NAT_MEAN))
}
if (!is.null(ms)) {
  n_all <- unique(ms$n[ms$sample_set == "all"]); n_s2 <- unique(ms$n[ms$sample_set == "s2"])
  check("Sample-level metrics: n of the full set = AGC samples", identical(as.integer(n_all), nrow(agc)),
        paste(n_all, collapse = ", "))
  if (length(n_s2) > 0) check("Sample-level metrics: n of the Sentinel-2 set = samples with has_s2",
                              identical(as.integer(n_s2), as.integer(N_S2)), sprintf("%s vs %s", paste(n_s2, collapse = ", "), N_S2))
}
if (!is.null(bp) && "AGC_S2" %in% bp$model)
  check("Tuning n of the Sentinel-2 models = samples with has_s2",
        all(bp$n[bp$model %in% c("AGC_GSE_SUB", "AGC_S2", "AGC_S2MED")] == N_S2), sprintf("%s", N_S2))
if (!is.null(iu) && !is.null(im)) {
  nu <- iu |> dplyr::count(sample_set, level, name = "n_iu") |>
    dplyr::inner_join(im |> dplyr::distinct(sample_set, level, n_units), by = c("sample_set", "level"))
  check("Inventory units per set and level = n_units of P4_inventory_metrics", all(nu$n_iu == nu$n_units),
        paste(sprintf("%s/%s %d", nu$sample_set, nu$level, nu$n_iu), collapse = ", "))
}

# =============================================================================
# Table 1: field data by ecoregion
# =============================================================================
step("Table 1")
cnt <- src("P4_counts_ecoregion_year.csv")
if (!is.null(eco)) {
  yrs <- if (!is.null(cnt)) {
    yc <- intersect(as.character(TAB$years), names(cnt))
    tibble::tibble(ecoregion = cnt$ecoregion,
                   years = apply(cnt[, yc], 1, function(v) paste(yc[v > 0], collapse = ", ")))
  } else tibble::tibble(ecoregion = character(0), years = character(0))
  clu <- agc |> dplyr::group_by(ecoregion) |>
    dplyr::summarise(clusters = dplyr::n_distinct(g_cluster[grepl("^c", g_cluster)]),
                     noise = sum(!grepl("^c", g_cluster)), .groups = "drop")
  units <- if (!is.null(iu)) iu |> dplyr::filter(sample_set == "all", level %in% c("site", "meadow")) |>
    dplyr::count(ecoregion, level) |> tidyr::pivot_wider(names_from = level, values_from = n, values_fill = 0) else
      tibble::tibble(ecoregion = character(0), site = integer(0), meadow = integer(0))
  for (k in c("site", "meadow")) if (!k %in% names(units)) units[[k]] <- 0L
  t1 <- eco |>
    dplyr::left_join(yrs, by = "ecoregion") |> dplyr::left_join(clu, by = "ecoregion") |>
    dplyr::left_join(units, by = "ecoregion") |>
    dplyr::arrange(eco_id == 0, dplyr::desc(area_km2))
  n_str <- sum(t1$t2_stratum %in% TRUE)
  tot_row <- tibble::tibble(
    Ecoregion = "Indonesia", a = sum(t1$area_km2), sa = 100, n = nrow(agc), sn = 100,
    years = paste(sort(unique(agc$year)), collapse = ", "),
    clusters = dplyr::n_distinct(agc$g_cluster[grepl("^c", agc$g_cluster)]),
    noise = sum(!grepl("^c", agc$g_cluster)), site = sum(t1$site, na.rm = TRUE), meadow = sum(t1$meadow, na.rm = TRUE),
    fm = NAT_MEAN, fs = sd(agc$AGC_pred),
    str = sprintf("%d strata, %s %% of extent", n_str, num(100 * sum(t1$area_km2[t1$t2_stratum %in% TRUE]) / sum(t1$area_km2), 1)))
  T1 <- dplyr::bind_rows(
    t1 |> dplyr::transmute(Ecoregion = ecoregion, a = area_km2, sa = share_area_pct, n = n_agc,
                           sn = share_samples_pct, years, clusters, noise, site, meadow,
                           fm = field_mean_gC_m2, fs = field_sd_gC_m2,
                           str = ifelse(t2_stratum %in% TRUE, "ecoregion mean",
                                        sprintf("national mean (n < %d)", CFG$min_n_eco))),
    tot_row) |>
    dplyr::transmute(Ecoregion,
                     !!sprintf("Extent (%s)", KM2) := num(a, 1),
                     `Extent share (%)` = share_lab(sa),
                     `AGC samples (n)` = num(dplyr::coalesce(n, 0), 0),
                     `Sample share (%)` = num(dplyr::coalesce(sn, 0), 1),
                     `Survey years` = ifelse(is.na(years) | years == "", "–", years),
                     `Meadow clusters (n)` = num(dplyr::coalesce(clusters, 0L), 0),
                     `Unclustered samples (n)` = num(dplyr::coalesce(noise, 0L), 0),
                     `Sites tested (n)` = num(dplyr::coalesce(site, 0L), 0),
                     `Meadows tested (n)` = num(dplyr::coalesce(meadow, 0L), 0),
                     !!sprintf("Field AGC, mean ± SD (%s)", GCM2) := pm(fm, fs, 1),
                     `Tier 2 stratum` = str)
  derive("t2_n_strata", n_str, sprintf("ecoregions with >= %d AGC samples (Tier 2 strata)", CFG$min_n_eco))
  derive("n_sites_tested", tot_row$site, "sites (100 m blocks, n >= CFG$min_n_site in the main survey year), all AGC samples")
  derive("n_meadows_tested", tot_row$meadow, "meadows (DBSCAN clusters, n >= CFG$min_n_meadow), all AGC samples")
  derive("n_meadow_clusters", tot_row$clusters, "DBSCAN clusters of AGC samples")
  save_table("data",
             "Seagrass extent, field samples, the sites and meadows tested, and the Tier 2 stratum of each marine ecoregion (MEOW).",
             T1, c(sprintf("Extent: final map extent (core + add-on), %s %s.", num(AREA, 2), KM2),
                   sprintf("Meadow clusters: DBSCAN clusters of AGC samples (eps %s m, minPts %s). Unclustered samples: DBSCAN noise points, each its own group in the cluster design.",
                           CFG$dbscan_eps, CFG$dbscan_minPts),
                   sprintf("Sites tested: 100 m blocks with at least %d AGC samples in their main survey year. Meadows tested: clusters with at least %d AGC samples.",
                           CFG$min_n_site, CFG$min_n_meadow),
                   sprintf("Tier 2 stratum: the ecoregion field mean where n >= %d, otherwise the national field mean.", CFG$min_n_eco)))
}

# =============================================================================
# Table 2: AGC totals by method, by ecoregion and for Indonesia
# =============================================================================
step("Table 2")
lyr <- src("P4_tier2_eco_loyo.csv")
if (!is.null(eco) && !is.null(px_nat)) {
  t2 <- eco |> dplyr::arrange(eco_id == 0, dplyr::desc(area_km2)) |>
    dplyr::mutate(diff_pct = 100 * (px_total_tC - t2_total_tC) / t2_total_tC)
  nat_px <- c(mean(px_nat$total_tC), min(px_nat$total_tC), max(px_nat$total_tC)) / 1e3
  t1n <- if (!is.null(t1b)) t1b |> dplyr::filter(response == "AGC_pred") else NULL
  t1_txt <- if (!is.null(t1n)) rng(t1n$total_tC / 1e3, t1n$ci_low_tC / 1e3, t1n$ci_up_tC / 1e3, 1) else
    num(sum(eco$t1_national_factor_tC) / 1e3, 1)
  t2_nat <- if (!is.null(t2n)) t2n[1, ] else NULL
  t2_txt <- if (!is.null(t2_nat)) rng(t2_nat$total_tC / 1e3, t2_nat$ci_low_tC / 1e3, t2_nat$ci_up_tC / 1e3, 1) else
    num(sum(eco$t2_total_tC) / 1e3, 1)
  T2 <- dplyr::bind_rows(
    t2 |> dplyr::transmute(
      Ecoregion = ecoregion,
      !!sprintf("Extent (%s)", KM2) := num(area_km2, 1),
      `AGC samples (n)` = num(dplyr::coalesce(n_agc, 0), 0),
      !!sprintf("Field mean (%s)", GCM2) := num(field_mean_gC_m2, 1),
      !!sprintf("Map mean (%s)", GCM2) := num(px_mean_gC_m2, 1),
      `Pixel-based GSE, mean [min, max] (kt C)` = rng(px_total_tC / 1e3, px_total_min_tC / 1e3, px_total_max_tC / 1e3, 2),
      `Tier 1, national mean (kt C)` = num(t1_national_factor_tC / 1e3, 2),
      `Tier 2 [95 % CI] (kt C)` = paste0(rng(t2_total_tC / 1e3, (t2_total_tC - t2_ci_half_tC) / 1e3,
                                             (t2_total_tC + t2_ci_half_tC) / 1e3, 2),
                                         ifelse(t2_stratum %in% TRUE, "", " *")),
      `Pixel-based vs Tier 2 (%)` = num(diff_pct, 1)),
    tibble::tibble(Ecoregion = "Indonesia",
                   !!sprintf("Extent (%s)", KM2) := num(sum(eco$area_km2), 1),
                   `AGC samples (n)` = num(sum(eco$n_agc, na.rm = TRUE), 0),
                   !!sprintf("Field mean (%s)", GCM2) := num(NAT_MEAN, 1),
                   !!sprintf("Map mean (%s)", GCM2) := num(sum(eco$px_total_tC) / sum(eco$area_km2), 1),
                   `Pixel-based GSE, mean [min, max] (kt C)` = rng(nat_px[1], nat_px[2], nat_px[3], 1),
                   `Tier 1, national mean (kt C)` = t1_txt,
                   `Tier 2 [95 % CI] (kt C)` = t2_txt,
                   `Pixel-based vs Tier 2 (%)` = num(100 * (nat_px[1] * 1e3 - sum(eco$t2_total_tC)) / sum(eco$t2_total_tC), 1)))
  str_rows <- t2 |> dplyr::filter(t2_stratum %in% TRUE)
  derive("px_total_mean_tC", nat_px[1] * 1e3, "mean of the annual pixel-based national totals 2017-2024")
  derive("px_total_min_tC", nat_px[2] * 1e3, "lowest annual pixel-based national total")
  derive("px_total_max_tC", nat_px[3] * 1e3, "highest annual pixel-based national total")
  derive("px_vs_t2_absdiff_max_pct", max(abs(str_rows$diff_pct)), "largest |pixel - Tier 2| / Tier 2 x 100 among the Tier 2 strata")
  derive("px_vs_t2_absdiff_median_pct", stats::median(abs(str_rows$diff_pct)), "median |pixel - Tier 2| / Tier 2 x 100 among the Tier 2 strata")
  derive("px_vs_t2_national_pct", 100 * (nat_px[1] * 1e3 - sum(eco$t2_total_tC)) / sum(eco$t2_total_tC),
         "(mean pixel-based total - Tier 2 total) / Tier 2 total x 100, Indonesia")
  pci <- px_nat$ci_half_tC[px_nat$field_data %in% TRUE][1]
  notes <- c("Pixel-based GSE: mean of the annual totals 2017-2024, with the lowest and highest annual total; map mean = mean total / extent.",
             if (!is.na(pci)) sprintf("Pixel-based 95 %% CI of each national annual total: ± %s kt C (1.96 x cluster-robust SE of the mean out-of-fold residual x extent; Table S10).", num(pci / 1e3, 1)) else "",
             "Tier 1: national field mean x extent; national 95 % CI with n / DEFF as the effective sample size.",
             sprintf("Tier 2: ecoregion field mean x ecoregion extent where n >= %d (stratum); * national field mean (fewer than %d samples). 95 %% CI: stratified, n / DEFF per stratum, with the national-mean areas as one pooled term (Table S9).",
                     CFG$min_n_eco, CFG$min_n_eco))
  if (!is.null(t2_nat)) notes <- c(notes, sprintf("Tier 2 when one survey year is left out of the strata: %s kt C.",
                                                  span(t2_nat$loyo_min_tC / 1e3, t2_nat$loyo_max_tC / 1e3, 1)))
  if (!is.null(t1b)) {
    lo <- t1b |> dplyr::filter(response == "AGC_low"); hi <- t1b |> dplyr::filter(response == "AGC_up")
    if (nrow(lo) == 1 && nrow(hi) == 1) {
      notes <- c(notes, sprintf("Tier 1 from the lower and upper AGC bounds of the cover-to-AGC models: %s and %s kt C.",
                                num(lo$total_tC / 1e3, 1), num(hi$total_tC / 1e3, 1)))
      derive("t1_agc_low_tC", lo$total_tC, "Tier 1 pooled with AGC_low (P4_tier1_agc_bounds_gee)")
      derive("t1_agc_up_tC", hi$total_tC, "Tier 1 pooled with AGC_up (P4_tier1_agc_bounds_gee)")
    }
  }
  if (!is.null(awb)) {
    a <- awb |> dplyr::filter(grepl("30", ecoregions)) |> dplyr::slice_head(n = 1)
    if (nrow(a) == 1) notes <- c(notes, sprintf(
      "Bias of the pixel-based map when each ecoregion is held out, weighted by ecoregion extent (%d ecoregions, %s %% of the extent): %s %s, i.e. %s kt C nationally, %s the 95 %% CI of ± %s kt C.",
      a$n_ecoregions, num(a$area_covered_pct, 1), num(a$bias_pixel_gC_m2, 1), GCM2, num(a$bias_pixel_national_kt, 1),
      if (isTRUE(a$pixel_bias_inside_ci)) "inside" else "outside", num(a$pixel_ci_half_kt, 1)))
  }
  save_table("totals",
             "Seagrass AGC by marine ecoregion and for Indonesia from the pixel-based map (GSE), Tier 1 and Tier 2.",
             T2, notes)
}

# =============================================================================
# Table 3: inventory-level validation
# =============================================================================
step("Table 3")
if (!is.null(im)) {
  T3d <- im |>
    dplyr::mutate(sample_set = factor(sample_set, levels = names(SET_NAME)),
                  level = factor(level, levels = names(LEVEL_NAME)),
                  method = factor(method, levels = METHODS)) |>
    dplyr::arrange(sample_set, level, method) |>
    dplyr::mutate(method = as.character(method)) |>
    collapse_same(c("sample_set", "level"), c("RMSE_mean", "MEC_mean", "bias_pct_mean", "slope_mean"))
  first_in <- function(x) c(TRUE, x[-1] != x[-length(x)])
  T3 <- T3d |>
    dplyr::mutate(new_set = first_in(as.character(sample_set)),
                  new_lev = first_in(paste(sample_set, level))) |>
    dplyr::transmute(`Sample set` = ifelse(new_set, SET_NAME[as.character(sample_set)], ""),
                     Unit = ifelse(new_lev, LEVEL_NAME[as.character(level)], ""),
                     `Units (n)` = ifelse(new_lev, num(n_units, 0), ""),
                     `Samples (n)` = ifelse(new_lev, num(n_samples, 0), ""),
                     Method = method_label,
                     !!sprintf("RMSE (%s)", GCM2) := num(RMSE_mean, 1),
                     !!sprintf("MAE (%s)", GCM2) := num(MAE_mean, 1),
                     !!R2 := num(R2_mean, 2),
                     MEC = rng(MEC_mean, MEC_q025, MEC_q975, 2),
                     `Bias (%)` = num(bias_pct_mean, 1),
                     Slope = num(slope_mean, 2))
  for (i in seq_len(nrow(im))) {
    r <- im[i, ]
    key <- sprintf("inv_%s_%s_%s", r$sample_set, r$level, gsub("[^a-z0-9]+", "_", tolower(r$method)))
    derive(paste0(sub("_+$", "", key), "_MEC"), r$MEC_mean, "P4_inventory_metrics MEC_mean")
  }
  nrep <- im |> dplyr::distinct(level, n_repeats)
  save_table("inventory",
             "Accuracy of the mean AGC of sites, meadows and ecoregions held out of training, by estimation method.",
             T3, c(sprintf("Units: site = 100 m block, samples of its main survey year (n >= %d), held out in the 100 m block design; meadow = DBSCAN cluster (n >= %d), held out in the cluster design; ecoregion (n >= %d), held out in the ecoregion design.",
                           CFG$min_n_site, CFG$min_n_meadow, CFG$min_n_eco),
                   sprintf("Each unit mean is the mean of the out-of-fold predictions of its samples; metrics compare the unit means with the field unit means. Mean over %s repeats (ecoregion: %s seeds); MEC with the 2.5-97.5 percentile range.",
                           paste(unique(nrep$n_repeats[nrep$level != "ecoregion"]), collapse = ", "),
                           paste(unique(nrep$n_repeats[nrep$level == "ecoregion"]), collapse = ", ")),
                   sprintf("Sentinel-2 set: AGC samples of the %d-%d surveys with a clear Sentinel-2 L2A pixel (n = %s); every method is trained and tested on these samples only (Table S5).",
                           min(CFG$s2_years), max(CFG$s2_years), num(N_S2, 0)),
                   "Tier 1: training national mean. Tier 2: training ecoregion mean (n >= 30), else Tier 1. Local mean: training mean of the same meadow (n >= 10), else Tier 2. Methods equal by definition in a design share one row.",
                   "MEC: 1 - SSE / sum((obs - mean(obs))^2). Bias: (mean pred - mean obs) / mean obs x 100. Slope: cov(obs, pred) / var(obs) of the unit means (1 = no shrinkage toward the mean)."))
}

# =============================================================================
# S1: AGC samples by ecoregion and survey year
# =============================================================================
step("S1")
if (!is.null(cnt)) {
  yc <- intersect(as.character(TAB$years), names(cnt))
  S1 <- dplyr::bind_rows(cnt |> dplyr::arrange(dplyr::desc(total)),
                         cnt |> dplyr::summarise(dplyr::across(c(dplyr::all_of(yc), total), sum)) |>
                           dplyr::mutate(ecoregion = "Total")) |>
    dplyr::rename(Ecoregion = ecoregion, Total = total) |>
    dplyr::mutate(dplyr::across(-Ecoregion, ~ num(.x, 0)))
  save_table("samples_year", "AGC field samples by marine ecoregion and survey year.", S1,
             "No AGC samples were collected in 2024.")
}

# =============================================================================
# S2: presence model cross-validation (GEE script P4_G09)
# =============================================================================
step("S2")
pa_cv_f <- file.path(TAB$gee_dir, "P4_gee_pa_cv.csv"); pa_n_f <- file.path(TAB$gee_dir, "P4_gee_pa_cv_counts.csv")
if (file.exists(pa_cv_f)) {
  pcv <- readr::read_csv(pa_cv_f, show_col_types = FALSE)
  file.copy(pa_cv_f, file.path(TAB$dir, "inputs"), overwrite = TRUE)
  check("Presence CV: TP + FP + FN + TN = n at every threshold", all(pcv$TP + pcv$FP + pcv$FN + pcv$TN == pcv$n),
        sprintf("n = %s", paste(unique(pcv$n), collapse = ", ")))
  S2t <- pcv |> dplyr::arrange(threshold) |> dplyr::transmute(
    Threshold = num(threshold, 1), `Points (n)` = num(n, 0), Presences = num(n_pres, 0), Absences = num(n_abs, 0),
    TP = num(TP, 0), FP = num(FP, 0), FN = num(FN, 0), TN = num(TN, 0),
    Accuracy = num(accuracy, 3), Kappa = num(kappa, 3), Precision = num(precision, 3), Recall = num(recall, 3),
    F1 = num(f1, 3))
  cnt_note <- character(0)
  if (file.exists(pa_n_f)) {
    pn <- readr::read_csv(pa_n_f, show_col_types = FALSE)
    file.copy(pa_n_f, file.path(TAB$dir, "inputs"), overwrite = TRUE)
    pw <- pn |> tidyr::pivot_wider(names_from = PA01, values_from = n, names_prefix = "pa") |> dplyr::arrange(year_pa)
    cnt_note <- sprintf("Training points by year (presence / absence): %s; total %s / %s.",
                        paste(sprintf("%d %s / %s", pw$year_pa, num(pw$pa1, 0), num(pw$pa0, 0)), collapse = ", "),
                        num(sum(pw$pa1), 0), num(sum(pw$pa0), 0))
    check("Presence CV: training points by year add up to n", sum(pn$n) == pcv$n[1], sprintf("%d vs %d", sum(pn$n), pcv$n[1]))
  }
  save_table("pa_cv", "Five-fold cross-validation of the seagrass presence model at three probability thresholds.",
             S2t, c("Random forest in GEE (smileRandomForest: 500 trees, 10 variables per split, minimum leaf 3, bag fraction 0.8, seed 42) on the 64 GSE bands and ACA depth; training points 2018-2023; random folds (seed 42). The folds and settings are those of the presence model that defined the extent.",
                    "Presences: field quadrats with cover > 0; absences: quadrats with zero cover and non-seagrass areas of the Allen Coral Atlas.",
                    "The extent uses 0.7 (MaxProb over the eight annual probability layers).",
                    cnt_note))
} else cat(sprintf("  -> missing: %s (copy P4_gee_pa_cv.csv and P4_gee_pa_cv_counts.csv from Drive P4v2_GEE)\n", pa_cv_f))

# =============================================================================
# S3: AGC samples inside the extent
# =============================================================================
step("S3")
sie <- src("P4_samples_in_extent.csv")
if (!is.null(sie)) {
  S3 <- sie |> dplyr::mutate(level = factor(level, levels = c("all", "ecoregion", "year")),
                             key = ifelse(level == "ecoregion", -n, 0)) |>
    dplyr::arrange(level, key, group) |>
    dplyr::transmute(Group = dplyr::case_when(level == "all" ~ "All AGC samples",
                                              level == "year" ~ paste("Survey year", group), TRUE ~ group),
                     `AGC samples (n)` = num(n, 0), `In final extent (%)` = num(in_extent_pct, 1),
                     `MaxProb < 0.7 (%)` = num(fail_prob_pct, 1), `Depth outside 0-5 m (%)` = num(fail_depth_pct, 1),
                     !!sprintf("%s (%%)", TAB$trf_label) := num(class16_pct, 1),
                     !!sprintf("Mean AGC inside (%s)", GCM2) := num(mean_agc_in, 1),
                     !!sprintf("Mean AGC outside (%s)", GCM2) := num(mean_agc_out, 1))
  save_table("in_extent", "AGC samples that fall inside the final extent, and the masks that exclude the others.", S3,
             c("A sample can fail more than one mask, so the mask columns do not add up to the share outside.",
               sprintf("Final extent = core + add-on: shallow (ACA depth 0-5 m), MaxProb >= 0.7 (maximum over the eight annual layers), not %s (ACA geomorphic map).", TAB$trf_label),
               "Every AGC sample is used in the models and the validation, inside the extent or not."))
}

# =============================================================================
# S4: extent under the geomorphic mask options and the KKP 2025 area
# =============================================================================
step("S4")
esn <- src("P4_extent_sensitivity_national.csv"); ese <- src("P4_extent_sensitivity_ecoregion.csv")
if (!is.null(esn) && !is.null(ese) && !is.null(px_nat)) {
  core <- esn$area_km2_at_10m[esn$mask == "current"]
  f30  <- core / sum(ese$current_km2)              # 30 m area -> 10 m map area (as P4_07)
  only16_est <- f30 * sum(ese$only16_km2)
  regional <- f30 * (sum(ese$only16_km2) + sum((ese$nogeo_km2 - ese$only16_km2)[!ese$ecoregion %in% TAB$turbid_eco]))
  nogeo <- f30 * sum(ese$nogeo_km2)
  dens <- px_nat$mean_gC_m2                         # annual national mean density of the map
  opt <- tibble::tibble(
    Option = c(sprintf("Core extent: ACA geomorphic data present, %s removed", TAB$trf_label),
               sprintf("Final extent: core + add-on (%s removed everywhere)", TAB$trf_label),
               sprintf("%s removed only in %s", TAB$trf_label, paste(TAB$turbid_eco, collapse = ", ")),
               "No geomorphic mask",
               "KKP 2025 national seagrass map"),
    Basis = c("mapped (10 m)", "mapped (10 m)", "estimated from 30 m", "estimated from 30 m", "official map"),
    km2 = c(core, AREA, regional, nogeo, TAB$kkp_km2))
  S4 <- opt |> dplyr::transmute(
    Option, Basis,
    !!sprintf("Extent (%s)", KM2) := num(km2, 0),
    `Share of KKP 2025 (%)` = num(100 * km2 / TAB$kkp_km2, 1),
    `Tier 1, pooled mean (kt C)` = num(NAT_MEAN * km2 / 1e3, 1),
    `Map density x extent, range 2017-2024 (kt C)` = span(min(dens) * km2 / 1e3, max(dens) * km2 / 1e3, 1))
  derive("extent_core_km2", core, "core extent (P4_extent_sensitivity_national, current, at 10 m)")
  derive("extent_trf_regional_km2", regional, "Terrestrial Reef Flat removed only in TAB$turbid_eco, scaled to 10 m")
  derive("extent_nogeo_km2", nogeo, "no geomorphic mask, scaled to 10 m")
  derive("tier1_kkp_tC", NAT_MEAN * TAB$kkp_km2, "national field mean x KKP 2025 area")
  derive("map_kkp_min_tC", min(dens) * TAB$kkp_km2, "lowest annual national map density x KKP 2025 area")
  derive("map_kkp_max_tC", max(dens) * TAB$kkp_km2, "highest annual national map density x KKP 2025 area")
  save_table("extent_kkp",
             "Seagrass extent under the geomorphic mask options, compared with the KKP 2025 national map, and the national AGC each extent implies.",
             S4, c("All options share MaxProb >= 0.7 and ACA depth 0-5 m; they differ only in the ACA geomorphic mask.",
                   sprintf("Estimated rows: areas of the 30 m sensitivity run (P4_G05) scaled to the 10 m map by the core ratio %.4f. The same scaling gives %s %s for the final extent, %.1f %% from the mapped %s %s.",
                           f30, num(only16_est, 0), KM2, 100 * (only16_est - AREA) / AREA, num(AREA, 0), KM2),
                   sprintf("Tier 1, pooled mean: national field mean (%.2f %s) x extent.", NAT_MEAN, GCM2),
                   sprintf("Map density x extent: annual national mean density of the AGC map (%.1f-%.1f %s) x extent; for the final extent this equals the pixel-based totals.",
                           min(dens), max(dens), GCM2)))
}

# =============================================================================
# S5: Sentinel-2 availability at the AGC samples
# =============================================================================
step("S5")
s2a <- src("P4_s2_availability.csv")
if (!is.null(s2a)) {
  s2a <- s2a |> dplyr::mutate(year = as.character(year))
  yr <- s2a |> dplyr::filter(ecoregion == "All", year != "All") |>
    dplyr::left_join(TAB$s2_scenes |> dplyr::mutate(year = as.character(year)), by = "year") |>
    dplyr::arrange(year)
  er <- s2a |> dplyr::filter(ecoregion != "All", year %in% as.character(CFG$s2_years)) |>
    dplyr::group_by(ecoregion) |>
    dplyr::summarise(n_agc = sum(n_agc), n_has_s2 = sum(n_has_s2), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(n_agc))
  # medians of clear observations per ecoregion from the samples, when P4_samples.csv has S2_nobs
  if ("S2_nobs" %in% names(agc)) {
    mo <- agc |> dplyr::filter(year %in% CFG$s2_years) |> dplyr::group_by(ecoregion) |>
      dplyr::summarise(median_nobs = stats::median(S2_nobs, na.rm = TRUE), .groups = "drop")
    er <- er |> dplyr::left_join(mo, by = "ecoregion")
  } else er$median_nobs <- NA_real_
  all_s2 <- er |> dplyr::summarise(n_agc = sum(n_agc), n_has_s2 = sum(n_has_s2))
  S5 <- dplyr::bind_rows(
    yr |> dplyr::transmute(Group = paste("Survey year", year), `AGC samples (n)` = num(n_agc, 0),
                           `L2A scenes` = num(l2a_scenes, 0), `L1C scenes` = num(l1c_scenes, 0),
                           `In the Sentinel-2 set (n)` = num(n_has_s2, 0),
                           `In the Sentinel-2 set (%)` = num(pct_has_s2, 1),
                           `Clear observations, median (n)` = num(median_nobs, 0)),
    tibble::tibble(Group = sprintf("Surveys %d-%d", min(CFG$s2_years), max(CFG$s2_years)),
                   `AGC samples (n)` = num(all_s2$n_agc, 0), `L2A scenes` = "", `L1C scenes` = "",
                   `In the Sentinel-2 set (n)` = num(all_s2$n_has_s2, 0),
                   `In the Sentinel-2 set (%)` = num(100 * all_s2$n_has_s2 / all_s2$n_agc, 1),
                   `Clear observations, median (n)` = ""),
    er |> dplyr::transmute(Group = paste0("  ", ecoregion), `AGC samples (n)` = num(n_agc, 0),
                           `L2A scenes` = "", `L1C scenes` = "",
                           `In the Sentinel-2 set (n)` = num(n_has_s2, 0),
                           `In the Sentinel-2 set (%)` = num(100 * n_has_s2 / n_agc, 1),
                           `Clear observations, median (n)` = num(median_nobs, 0)))
  check("Sentinel-2 set in P4_s2_availability = samples with has_s2", identical(as.integer(all_s2$n_has_s2), as.integer(N_S2)),
        sprintf("%s vs %s", all_s2$n_has_s2, N_S2))
  derive("s2_l2a_scenes_2017", TAB$s2_scenes$l2a_scenes[TAB$s2_scenes$year == 2017], "P4_G11 console")
  derive("s2_l1c_scenes_2017", TAB$s2_scenes$l1c_scenes[TAB$s2_scenes$year == 2017], "P4_G11 console")
  derive("s2_l2a_scenes_2018", TAB$s2_scenes$l2a_scenes[TAB$s2_scenes$year == 2018], "P4_G11 console")
  derive("s2_l1c_scenes_2018", TAB$s2_scenes$l1c_scenes[TAB$s2_scenes$year == 2018], "P4_G11 console")
  save_table("s2_avail",
             "Sentinel-2 L2A surface reflectance at the AGC samples: scenes per survey year, samples with a clear pixel, and the Sentinel-2 set by ecoregion.",
             S5, c("Scenes: Sentinel-2 scenes over the AGC samples of the survey year in the Earth Engine collections COPERNICUS/S2_SR_HARMONIZED (L2A) and COPERNICUS/S2_HARMONIZED (L1C), printed by P4_G11.",
                   sprintf("Annual mosaic: calendar year, every scene, Cloud Score+ cs_cdf >= 0.6 per pixel, one observation per pass; percentiles 10, 25, 50, 75 and 90 of B2, B3, B4, B8 and B11 (25 predictors). A sample enters the Sentinel-2 set with >= %d clear observation in its survey year.", CFG$s2_min_obs),
                   sprintf("Surveys of %s are left out of the Sentinel-2 comparison because the L2A archive in Earth Engine holds few scenes for these years.",
                           paste(setdiff(sort(unique(agc$year)), CFG$s2_years), collapse = " and ")),
                   "Clear observations: median number of clear passes per sample in its survey year (by ecoregion: survey years of the Sentinel-2 set)."))
}

# =============================================================================
# S6: random forest settings of the four AGC models
# =============================================================================
step("S6")
tg <- dplyr::bind_rows(lapply(c("AGC", "AGC_GSE_SUB", "AGC_S2", "AGC_S2MED"), function(m) {
  f <- sprintf("P4_tuning_%s.csv", m)
  if (have(f)) src(f) |> dplyr::mutate(model = m) else NULL
}))
if (!is.null(bp)) {
  if (nrow(tg) > 0) {
    tgr <- tg |> dplyr::group_by(model) |>
      dplyr::summarise(n_set = dplyr::n(), rmse_min = min(RMSE), rmse_max = max(RMSE),
                       mtry_set = paste(sort(unique(mtry)), collapse = ", "), .groups = "drop")
    for (i in seq_len(nrow(tgr))) {
      derive(sprintf("tuning_%s_n_settings", tgr$model[i]), tgr$n_set[i], "rows of the P4_02 tuning grid")
      derive(sprintf("tuning_%s_rmse_range", tgr$model[i]), tgr$rmse_max[i] - tgr$rmse_min[i], "max - min RMSE over the grid")
    }
    bp6 <- bp |> dplyr::left_join(tgr, by = "model")
  } else bp6 <- bp |> dplyr::mutate(n_set = NA_real_, rmse_min = NA_real_, rmse_max = NA_real_, mtry_set = NA_character_)
  bp6 <- bp6 |> dplyr::mutate(model = factor(model, levels = names(MODEL_NAME))) |> dplyr::arrange(model)
  S6 <- bp6 |> dplyr::transmute(
    Model = MODEL_NAME[as.character(model)],
    `Predictors (n)` = num(n_predictors, 0), `Training samples (n)` = num(n, 0), Trees = num(num_trees, 0),
    `mtry grid` = dplyr::coalesce(mtry_set, "–"),
    `Variables per split (mtry)` = num(mtry, 0), `Minimum node size` = num(min_node_size, 0),
    `Sample fraction` = num(sample_fraction, 1),
    `Tuning CV RMSE` = num(RMSE, 2), `Tuning CV MAE` = num(MAE, 2),
    !!paste("Tuning CV", R2) := num(R2, 3), `Tuning CV MEC` = num(MEC, 3), `Tuning CV bias (%)` = num(bias_pct, 2),
    `Settings in the grid (n)` = num(n_set, 0),
    `RMSE range over the grid` = span(rmse_min, rmse_max, 2))
  save_table("rf_settings",
             "Random forest settings of the four AGC models: GSE on all AGC samples, and GSE, Sentinel-2 percentiles and Sentinel-2 median on the Sentinel-2 set.",
             S6, c(sprintf("Settings chosen by %d-fold CV grouped by 10 m pixel over mtry (grid in the table), minimum node size {%s} and sample fraction {%s}; the folds are the same for the four models.",
                           CFG$tune_folds, paste(CFG$grid_node, collapse = ", "), paste(CFG$grid_frac, collapse = ", ")),
                   "The mtry grid adds floor(sqrt(p)) and floor(p / 3) where the number of predictors p is smaller than the largest grid value.",
                   "GSE: 64 bands of the annual Google Satellite Embedding. Sentinel-2: B2, B3, B4, B8 and B11 of the L2A annual mosaic, percentiles 10-90 (25 predictors) or the median only (5).",
                   "The GSE settings of all AGC samples are also used in GEE for the national maps (smileRandomForest: numberOfTrees, variablesPerSplit, minLeafPopulation, bagFraction)."))
}

# =============================================================================
# S7: national AGC by year and method
# =============================================================================
step("S7")
if (!is.null(nat3)) {
  w <- function(m) nat3 |> dplyr::filter(method == m) |> dplyr::select(year, total_kt, ci_low_kt, ci_up_kt)
  px <- w("Pixel-based GSE"); t1y <- w("Tier 1 (survey-year mean)"); t1p <- w("Tier 1 (pooled mean)")
  t2y <- w("Tier 2 (ecoregion strata)")
  ny <- agc |> dplyr::count(year, name = "n")
  neff <- if (!is.null(t1g)) t1g |> dplyr::filter(year != "pooled") |>
    dplyr::transmute(year = as.integer(year), n_eff) else tibble::tibble(year = integer(0), n_eff = numeric(0))
  b <- tibble::tibble(year = TAB$years) |>
    dplyr::left_join(ny, by = "year") |> dplyr::left_join(neff, by = "year") |>
    dplyr::left_join(px |> dplyr::rename(px = total_kt, pxl = ci_low_kt, pxu = ci_up_kt), by = "year") |>
    dplyr::left_join(t1y |> dplyr::rename(ty = total_kt, tyl = ci_low_kt, tyu = ci_up_kt), by = "year") |>
    dplyr::left_join(t1p |> dplyr::rename(tp = total_kt, tpl = ci_low_kt, tpu = ci_up_kt), by = "year") |>
    dplyr::left_join(t2y |> dplyr::rename(t2 = total_kt, t2l = ci_low_kt, t2u = ci_up_kt), by = "year")
  S7 <- b |> dplyr::transmute(
    Year = as.character(year),
    `AGC samples (n)` = ifelse(is.na(n), "0", num(n, 0)),
    `Pixel-based GSE` = paste0(rng(px, pxl, pxu, 1), ifelse(year %in% ny$year, "", " †")),
    `Tier 1, survey-year mean` = ifelse(is.na(ty), "–",
                                        ifelse(is.na(tyl), paste0(num(ty, 1), " ‡"), rng(ty, tyl, tyu, 1))),
    `Tier 1, pooled mean` = rng(tp, tpl, tpu, 1),
    `Tier 2, ecoregion strata` = rng(t2, t2l, t2u, 1))
  sm <- function(x, f) num(f(x, na.rm = TRUE), 1)
  S7 <- dplyr::bind_rows(S7,
    tibble::tibble(Year = "Mean 2017–2024", `AGC samples (n)` = num(nrow(agc), 0),
                   `Pixel-based GSE` = sm(b$px, mean), `Tier 1, survey-year mean` = "–",
                   `Tier 1, pooled mean` = rng(b$tp[1], b$tpl[1], b$tpu[1], 1),
                   `Tier 2, ecoregion strata` = rng(b$t2[1], b$t2l[1], b$t2u[1], 1)),
    tibble::tibble(Year = "Range of annual values", `AGC samples (n)` = "",
                   `Pixel-based GSE` = span(min(b$px, na.rm = TRUE), max(b$px, na.rm = TRUE), 1),
                   `Tier 1, survey-year mean` = span(min(b$ty, na.rm = TRUE), max(b$ty, na.rm = TRUE), 1),
                   `Tier 1, pooled mean` = "constant", `Tier 2, ecoregion strata` = "constant"))
  derive("t1_year_min_tC", 1e3 * min(b$ty, na.rm = TRUE), "lowest Tier 1 survey-year total")
  derive("t1_year_max_tC", 1e3 * max(b$ty, na.rm = TRUE), "highest Tier 1 survey-year total")
  notes7 <- c(sprintf("All methods use the same fixed extent of %s %s (2017-2024).", num(AREA, 2), KM2),
              "Pixel-based CI: 1.96 x cluster-robust SE of the mean out-of-fold residual (cluster design) x extent (Table S10).",
              "Tier 1 and Tier 2 CIs use n / DEFF as the effective sample size; Tier 2 is stratified by ecoregion (Table S9).",
              "Tier 1 and Tier 2 pooled over the survey years are constant because the extent is fixed.",
              "† No field data in 2024: the CI assumes the 2017-2023 bias SE.",
              if (any(!is.na(b$ty) & is.na(b$tyl))) sprintf("‡ Effective sample size below 2 (n = %s): no CI.",
                                                          paste(b$n[!is.na(b$ty) & is.na(b$tyl)], collapse = ", ")) else "")
  if (!is.null(t2n)) notes7 <- c(notes7, sprintf("Tier 2 with one survey year left out of the strata: %s kt C.",
                                                 span(t2n$loyo_min_tC[1] / 1e3, t2n$loyo_max_tC[1] / 1e3, 1)))
  save_table("national_year", "National seagrass AGC (kt C) by year and estimation method, with 95 % confidence intervals.",
             S7, notes7)
}

# =============================================================================
# S8: Tier 1 by survey year and the survey composition
# =============================================================================
step("S8")
yco <- src("P4_year_composition.csv")
if (!is.null(t1g)) {
  comp <- if (!is.null(yco)) yco |> dplyr::group_by(year) |> dplyr::arrange(dplyr::desc(share), .by_group = TRUE) |>
    dplyr::summarise(main = paste(sprintf("%s (%.0f %%)", ecoregion[seq_len(min(2, dplyr::n()))],
                                          100 * share[seq_len(min(2, dplyr::n()))]), collapse = "; "),
                     .groups = "drop") |> dplyr::mutate(year = as.character(year)) else
      tibble::tibble(year = character(0), main = character(0))
  S8 <- t1g |> dplyr::mutate(year = as.character(year)) |> dplyr::left_join(comp, by = "year") |>
    dplyr::transmute(`Survey year` = ifelse(year == "pooled", "Pooled", year), `AGC samples (n)` = num(n, 0),
                     `Effective n` = num(n_eff, 1),
                     !!sprintf("Mean ± SD (%s)", GCM2) := pm(mean_gC_m2, sd_gC_m2, 2),
                     !!sprintf("95 %% CI half-width (%s)", GCM2) := num(ci_half_gC_m2, 2),
                     `Main ecoregions (share of the year's samples)` = ifelse(is.na(main), "all", main))
  save_table("tier1_year",
             "Field AGC by survey year (the Tier 1 factors) and the ecoregions that dominate each year's samples.",
             S8, c("Survey year and region are confounded: differences between years mostly reflect where the surveys took place.",
                   "Effective n = n / DEFF (Table S10). No CI where the effective sample size is below 2."))
}

# =============================================================================
# S9: Tier 2 strata
# =============================================================================
step("S9")
t2s <- src("P4_tier2_eco_strata.csv")
if (!is.null(t2e)) {
  lr <- if (!is.null(lyr)) lyr |> dplyr::filter(stratum %in% TRUE, ecoregion != "National") |>
    dplyr::group_by(ecoregion) |>
    dplyr::summarise(lo = min(mean_gC_m2), hi = max(mean_gC_m2), n_left = dplyr::n(), .groups = "drop") else
      tibble::tibble(ecoregion = character(0), lo = numeric(0), hi = numeric(0), n_left = integer(0))
  nat_lr <- if (!is.null(lyr)) lyr |> dplyr::filter(ecoregion == "National") else NULL
  ny9 <- if (!is.null(t2s)) t2s |> dplyr::select(ecoregion, n_years) else tibble::tibble(ecoregion = character(0), n_years = integer(0))
  if (!is.null(t2s)) {
    k <- t2e |> dplyr::filter(stratum %in% TRUE) |> dplyr::inner_join(t2s |> dplyr::select(ecoregion, m2 = mean_gC_m2), by = "ecoregion")
    check("Tier 2 stratum means in P4_tier2_by_ecoregion = P4_tier2_eco_strata", all(close_pct(k$mean_gC_m2, k$m2)),
          sprintf("%d strata", nrow(k)))
  }
  S9 <- t2e |> dplyr::left_join(ny9, by = "ecoregion") |> dplyr::left_join(lr, by = "ecoregion") |>
    dplyr::arrange(dplyr::desc(stratum %in% TRUE), dplyr::desc(n)) |>
    dplyr::transmute(
      Ecoregion = ecoregion, !!sprintf("Extent (%s)", KM2) := num(area_km2, 1),
      `AGC samples (n)` = num(n, 0), `Survey years (n)` = num(n_years, 0),
      !!sprintf("Field mean ± SD (%s)", GCM2) := pm(mean_gC_m2, sd_gC_m2, 1),
      `Effective n (n / DEFF)` = num(n_eff, 1),
      Stratum = ifelse(stratum %in% TRUE, "yes", "no: national mean"),
      !!sprintf("Density used (%s)", GCM2) := num(dens_gC_m2, 1),
      `Tier 2 total [95 % CI] (kt C)` = rng(t2_total_tC / 1e3, (t2_total_tC - t2_ci_half_tC) / 1e3,
                                            (t2_total_tC + t2_ci_half_tC) / 1e3, 2),
      !!sprintf("Leave-one-year-out mean range (%s)", GCM2) := span(lo, hi, 1))
  if (!is.null(t2n)) S9 <- dplyr::bind_rows(S9, tibble::tibble(
    Ecoregion = "Indonesia", !!sprintf("Extent (%s)", KM2) := num(sum(t2e$area_km2), 1),
    `AGC samples (n)` = num(sum(t2e$n), 0), `Survey years (n)` = num(dplyr::n_distinct(agc$year), 0),
    !!sprintf("Field mean ± SD (%s)", GCM2) := pm(NAT_MEAN, sd(agc$AGC_pred), 1),
    `Effective n (n / DEFF)` = num(sum(t2e$n_eff), 1),
    Stratum = sprintf("%d strata", t2n$n_strata[1]),
    !!sprintf("Density used (%s)", GCM2) := num(t2n$total_tC[1] / sum(t2e$area_km2), 1),
    `Tier 2 total [95 % CI] (kt C)` = rng(t2n$total_tC[1] / 1e3, t2n$ci_low_tC[1] / 1e3, t2n$ci_up_tC[1] / 1e3, 1),
    !!sprintf("Leave-one-year-out mean range (%s)", GCM2) := if (!is.null(nat_lr)) span(min(nat_lr$mean_gC_m2), max(nat_lr$mean_gC_m2), 1) else "–"))
  if (!is.null(t2n)) {
    derive("t2_area_strata_pct", t2n$area_strata_pct[1], "share of the extent in Tier 2 strata")
    derive("t2_area_fallback_km2", t2n$area_fallback_km2[1], "extent of the ecoregions with the national mean")
  }
  save_table("tier2_strata",
             "Tier 2 strata: field AGC by marine ecoregion (pooled over survey years), the density used for each ecoregion, and the Tier 2 totals.",
             S9, c(sprintf("Stratum: ecoregion with at least %d AGC samples; the other ecoregions use the national field mean, and their extent forms one pooled term of the CI.", CFG$min_n_eco),
                   "95 % CI: 1.96 x SD / sqrt(n / DEFF) x extent per stratum; national CI from the sum of the stratum variances and the pooled term.",
                   "Leave-one-year-out: stratum mean recomputed without one survey year at a time, over the left-out years where the ecoregion stays a stratum.",
                   "Strata are pooled over survey years and the extent is fixed, so Tier 2 is the same in every year."))
}

# =============================================================================
# S10: clustering, uncertainty parameters and checks of the totals
# =============================================================================
step("S10")
up <- src("P4_uncertainty_params.csv"); csm <- src("P4_cluster_summary.csv")
if (!is.null(up)) {
  lab <- c(ICC = "Intra-cluster correlation of out-of-fold residuals",
           adjusted_mean_cluster_size = "Adjusted mean cluster size (n0)",
           DEFF = "Design effect, 1 + (n0 - 1) x ICC", N_agc = "AGC samples", N_eff = "Effective sample size, n / DEFF",
           SD_model_gC_m2 = sprintf("SD_model: RMSE of the cluster design (%s)", GCM2),
           SD_data_gC_m2 = sprintf("SD_data: RMS of (AGC_up - AGC_low) / 3.92 (%s)", GCM2),
           SD_total_gC_m2 = sprintf("SD_total: sqrt(SD_model² + SD_data²) (%s)", GCM2),
           mean_resid_gC_m2 = sprintf("Mean out-of-fold residual, obs - pred, cluster design (%s)", GCM2),
           se_mean_resid_gC_m2 = sprintf("Cluster-robust SE of the mean residual (%s)", GCM2),
           bias_halfwidth_national_tC = "Pixel-based 95 % CI half-width, 1.96 x SE x extent (t C)")
  dg <- c(ICC = 3, adjusted_mean_cluster_size = 1, DEFF = 2, N_agc = 0, N_eff = 0, SD_model_gC_m2 = 2,
          SD_data_gC_m2 = 2, SD_total_gC_m2 = 2, mean_resid_gC_m2 = 2, se_mean_resid_gC_m2 = 2,
          bias_halfwidth_national_tC = 0)
  up <- up |> dplyr::filter(parameter %in% names(lab)) |> dplyr::distinct(parameter, .keep_all = TRUE)
  S10 <- up |> dplyr::transmute(Parameter = unname(lab[parameter]),
                                Value = mapply(function(v, p) num(v, dg[[p]]), value, parameter))
  if ("se_mean_resid_gC_m2" %in% up$parameter && !is.na(AREA)) {
    hw <- 1.96 * up$value[up$parameter == "se_mean_resid_gC_m2"] * AREA
    S10$Value[S10$Parameter == lab[["bias_halfwidth_national_tC"]]] <- num(hw, 0)
    derive("pixel_ci_half_tC", hw, "1.96 x cluster-robust SE of the mean residual x final extent")
  }
  if (!is.null(csm)) S10 <- dplyr::bind_rows(
    tibble::tibble(Parameter = c("DBSCAN eps (m) / minPts", "Meadow clusters (n)", "AGC samples in clusters (n)",
                                 "Unclustered AGC samples (n)", "Largest cluster (n samples)"),
                   Value = c(sprintf("%s / %s", num(csm$eps_m, 0), csm$minPts), num(csm$n_clusters, 0),
                             num(csm$n_clustered, 0), num(csm$n_noise, 0), num(csm$largest_cluster, 0))), S10)
  # ICC of the field AGC values themselves, and repeat surveys of the 100 m blocks
  clu <- agc |> dplyr::filter(grepl("^c", g_cluster))
  ic <- icc_oneway(clu$AGC_pred, clu$g_cluster)
  b100 <- agc |> dplyr::group_by(g_b100) |> dplyr::summarise(ny = dplyr::n_distinct(year), .groups = "drop")
  derive("icc_field_agc", ic$icc, "one-way ANOVA ICC of AGC_pred by g_cluster, clustered AGC samples (P4_samples.csv)")
  derive("deff_field_agc", ic$deff, "1 + (n0 - 1) x ICC of field AGC")
  derive("blocks100_with_agc", nrow(b100), "100 m blocks (g_b100) with AGC samples")
  derive("blocks100_multi_year", sum(b100$ny > 1), "100 m blocks with AGC samples from more than one survey year")
  S10 <- dplyr::bind_rows(S10, tibble::tibble(
    Parameter = c("ICC of field AGC (sensitivity)", "Design effect with the ICC of field AGC",
                  "100 m blocks with AGC samples (n)", "100 m blocks with AGC samples from more than one survey year (n)"),
    Value = c(num(ic$icc, 3), num(ic$deff, 2), num(nrow(b100), 0), num(sum(b100$ny > 1), 0))))
  if (!is.null(awb)) {
    a <- awb |> dplyr::arrange(dplyr::desc(grepl("30", ecoregions)))
    S10 <- dplyr::bind_rows(S10, purrr::map_dfr(seq_len(nrow(a)), function(i) tibble::tibble(
      Parameter = sprintf("Ecoregion held out, area-weighted bias of the pixel-based map, %s (%d ecoregions, %s %% of extent): %s / kt C nationally",
                          a$ecoregions[i], a$n_ecoregions[i], num(a$area_covered_pct[i], 1), GCM2),
      Value = sprintf("%s / %s", num(a$bias_pixel_gC_m2[i], 2), num(a$bias_pixel_national_kt[i], 1)))))
    a1 <- a[1, ]
    S10 <- dplyr::bind_rows(S10, tibble::tibble(
      Parameter = sprintf("Same for Tier 1 (%s): %s / kt C", a1$ecoregions, GCM2),
      Value = sprintf("%s / %s", num(a1$bias_tier1_gC_m2, 2), num(a1$bias_tier1_national_kt, 1))))
    derive("lolo_aw_bias_pixel_gC_m2", a1$bias_pixel_gC_m2, "area-weighted bias, ecoregion held out, pixel-based (n >= 30)")
    derive("lolo_aw_bias_pixel_kt", a1$bias_pixel_national_kt, "the same x extent, kt C")
    derive("lolo_aw_bias_tier1_kt", a1$bias_tier1_national_kt, "area-weighted bias, ecoregion held out, Tier 1, kt C")
  }
  if (!is.null(t1b)) {
    S10 <- dplyr::bind_rows(S10, purrr::map_dfr(seq_len(nrow(t1b)), function(i) tibble::tibble(
      Parameter = sprintf("Tier 1 pooled from %s: mean (%s) / total [95 %% CI] (kt C)", t1b$response[i], GCM2),
      Value = sprintf("%s / %s", num(t1b$mean_gC_m2[i], 2),
                      rng(t1b$total_tC[i] / 1e3, t1b$ci_low_tC[i] / 1e3, t1b$ci_up_tC[i] / 1e3, 1)))))
  }
  save_table("uncertainty",
             "Clustering of the AGC samples, the parameters of the confidence intervals, and checks of the totals.", S10,
             c(sprintf("CI half-width recomputed here with the final extent (%s %s).", num(AREA, 2), KM2),
               "The mean residual is not different from zero (|mean| < 1.96 x SE), so the totals are not bias-corrected.",
               "Area-weighted bias: mean (predicted - field) of each ecoregion when it is held out of training, weighted by the ecoregion extent.",
               "AGC_low and AGC_up: lower and upper AGC of each sample from the cover-to-AGC models; they bound the field-data uncertainty of Tier 1."))
}

# =============================================================================
# S11: sample-level accuracy by design, method and sample set
# =============================================================================
step("S11")
if (!is.null(ms)) {
  sh <- if (!is.null(dcd)) setNames(dcd$area_pct, dcd$design) else setNames(rep(NA_real_, 5), DESIGN_LEVELS)
  sa <- if (!is.null(dcd) && "agc_pct_mean" %in% names(dcd)) setNames(dcd$agc_pct_mean, dcd$design) else
    setNames(rep(NA_real_, 5), DESIGN_LEVELS)
  d11 <- ms |>
    dplyr::mutate(sample_set = factor(sample_set, levels = names(SET_NAME)),
                  design = factor(design, levels = DESIGN_LEVELS), method = factor(method, levels = METHODS)) |>
    dplyr::arrange(sample_set, design, method) |>
    dplyr::mutate(method = as.character(method)) |>
    collapse_same(c("sample_set", "design"), c("RMSE_mean", "MEC_mean", "bias_pct_mean"))
  first_in <- function(x) c(TRUE, x[-1] != x[-length(x)])
  S11 <- d11 |>
    dplyr::mutate(new_set = first_in(as.character(sample_set)), new_des = first_in(paste(sample_set, design)),
                  dz = as.character(design)) |>
    dplyr::transmute(`Sample set` = ifelse(new_set, sprintf("%s (n = %s)", SET_NAME[as.character(sample_set)], num(n, 0)), ""),
                     Design = ifelse(new_des, sprintf("%s (%s)", DESIGN_NAME[dz], DESIGN_RANGE[dz]), ""),
                     `Median distance to training` = ifelse(new_des, dist_lab(median_dist_m), ""),
                     `Share of extent (%)` = ifelse(new_des, share_lab(sh[dz]), ""),
                     `Share of pixel-based AGC (%)` = ifelse(new_des, share_lab(sa[dz]), ""),
                     Method = method_label,
                     !!sprintf("RMSE (%s)", GCM2) := rng(RMSE_mean, RMSE_q025, RMSE_q975, 1),
                     !!sprintf("MAE (%s)", GCM2) := rng(MAE_mean, MAE_q025, MAE_q975, 1),
                     !!R2 := rng(R2_mean, R2_q025, R2_q975, 3),
                     MEC = rng(MEC_mean, MEC_q025, MEC_q975, 3),
                     `Bias (%)` = rng(bias_pct_mean, bias_pct_q025, bias_pct_q975, 1))
  nrep <- ms |> dplyr::distinct(design, n_repeats)
  save_table("sample_level",
             "Out-of-fold accuracy of AGC at the sample level by validation design, estimation method and sample set.",
             S11, c(sprintf("Mean over %s repeats of %s-fold CV (ecoregion design: %s seeds, one fold per ecoregion), with the 2.5-97.5 percentile range in brackets.",
                            nrep$n_repeats[nrep$design == "pixel"][1], CFG$k_folds, nrep$n_repeats[nrep$design == "ecoregion"][1]),
                    "Distance range: distance from an extent pixel to the nearest AGC sample that this design represents; share of extent and of the mean annual pixel-based AGC in that range (Table S14).",
                    "Median distance to training: median distance from a test sample to the nearest training sample in that design (all AGC samples).",
                    "Methods equal by definition in a design share one row. Sentinel-2 set: Table S5."))
}

# =============================================================================
# S12: site-level accuracy by ecoregion
# =============================================================================
step("S12")
sbe <- src("P4_inventory_site_by_ecoregion.csv")
if (!is.null(sbe)) {
  ordr <- sbe |> dplyr::filter(sample_set == "all") |> dplyr::distinct(ecoregion, n_sites) |>
    dplyr::arrange(dplyr::desc(n_sites)) |> dplyr::pull(ecoregion)
  d12 <- sbe |>
    dplyr::mutate(sample_set = factor(sample_set, levels = names(SET_NAME)),
                  ecoregion = factor(ecoregion, levels = unique(c(ordr, unique(sbe$ecoregion)))),
                  method = factor(method, levels = METHODS)) |>
    dplyr::arrange(sample_set, ecoregion, method) |> dplyr::mutate(method = as.character(method))
  first_in <- function(x) c(TRUE, x[-1] != x[-length(x)])
  S12 <- d12 |>
    dplyr::mutate(new_set = first_in(as.character(sample_set)), new_eco = first_in(paste(sample_set, ecoregion))) |>
    dplyr::transmute(`Sample set` = ifelse(new_set, SET_NAME[as.character(sample_set)], ""),
                     Ecoregion = ifelse(new_eco, paste0(as.character(ecoregion), ifelse(reported, "", " *")), ""),
                     `Sites (n)` = ifelse(new_eco, num(n_sites, 0), ""),
                     Method = method,
                     !!sprintf("RMSE (%s)", GCM2) := num(RMSE, 1), !!sprintf("MAE (%s)", GCM2) := num(MAE, 1),
                     !!R2 := num(R2, 2), MEC = num(MEC, 2), `Bias (%)` = num(bias_pct, 1), Slope = num(slope, 2))
  save_table("site_eco",
             "Accuracy of the site means (100 m blocks held out) by marine ecoregion and estimation method.",
             S12, c(sprintf("Ecoregions with at least %d sites; * fewer sites, shown for completeness.", CFG$min_sites_eco),
                    "Metrics over the out-of-fold site means of the ecoregion in each CV repeat, averaged over the repeats; definitions as in Table 3.",
                    "MEC within one ecoregion compares the site means with the mean of the same sites, so a constant offset (bias) pushes MEC far below zero."))
}

# =============================================================================
# S13: accuracy for each ecoregion when it is held out
# =============================================================================
step("S13")
me <- src("P4_metrics_ecoregion.csv"); mew <- src("P4_metrics_ecoregion_equal_weight.csv")
lolo <- src("P4_lolo_region_means.csv")
if (!is.null(me)) {
  ordr <- me |> dplyr::distinct(ecoregion, n) |> dplyr::arrange(dplyr::desc(n)) |> dplyr::pull(ecoregion)
  d13 <- me |> dplyr::mutate(ecoregion = factor(ecoregion, levels = ordr), method = factor(method, levels = METHODS)) |>
    dplyr::arrange(ecoregion, method) |> dplyr::mutate(method = as.character(method)) |>
    collapse_same("ecoregion", c("RMSE", "MAE", "bias_pct"))
  fm <- if (!is.null(lolo)) setNames(lolo$field, lolo$ecoregion) else NULL
  first_in <- function(x) c(TRUE, x[-1] != x[-length(x)])
  S13 <- d13 |> dplyr::mutate(new_eco = first_in(as.character(ecoregion))) |>
    dplyr::transmute(Ecoregion = ifelse(new_eco, paste0(as.character(ecoregion), ifelse(reported, "", " *")), ""),
                     `AGC samples (n)` = ifelse(new_eco, num(n, 0), ""),
                     `Median distance to training (km)` = ifelse(new_eco, num(median_dist_km, 0), ""),
                     !!sprintf("Field mean (%s)", GCM2) := ifelse(new_eco & !is.null(fm), num(fm[as.character(ecoregion)], 1), ""),
                     Method = method_label,
                     !!sprintf("RMSE (%s)", GCM2) := num(RMSE, 1), !!sprintf("MAE (%s)", GCM2) := num(MAE, 1),
                     !!R2 := num(R2, 3), MEC = num(MEC, 3), `Bias (%)` = num(bias_pct, 1))
  notes <- c(sprintf("* n < %d: shown for completeness, not used in the summaries.", CFG$min_n_eco),
             "With the ecoregion held out, Tier 2 and the local mean fall back to the training national mean, so they equal Tier 1.",
             "R² is undefined (–) where a method predicts one value for the whole ecoregion (Tier 1).")
  if (!is.null(mew)) {
    mew <- mew |> dplyr::mutate(method = factor(method, levels = METHODS)) |> dplyr::arrange(method) |>
      dplyr::mutate(method = as.character(method)) |> collapse_same(character(0), c("MEC_median", "abs_bias_pct_median"))
    notes <- c(notes, paste0(
      "Across ecoregions (n >= ", CFG$min_n_eco, "): median MEC ",
      paste(sprintf("%s %s", mew$method_label, num(mew$MEC_median, 3)), collapse = ", "),
      "; median |bias| ", paste(sprintf("%s %s %%", mew$method_label, num(mew$abs_bias_pct_median, 1)), collapse = ", "), "."))
  }
  save_table("lolo", "Accuracy for each ecoregion when it is held out of training (ecoregion design, pooled over seeds).",
             S13, notes)
}

# =============================================================================
# S14: distance to AGC samples and between-year CV by ecoregion
# =============================================================================
step("S14")
dce <- src("P4_distance_coverage_eco.csv"); cvs <- src("P4_cv_summary.csv")
if (!is.null(dce)) {
  cvw <- if (!is.null(cvs)) cvs |>
    dplyr::mutate(grp = dplyr::case_when(cv_label == "<10" ~ "lt10", cv_label == "10-20" ~ "c1020", TRUE ~ "ge20")) |>
    dplyr::group_by(ecoregion, grp) |> dplyr::summarise(p = sum(area_pct), .groups = "drop") |>
    tidyr::pivot_wider(names_from = grp, values_from = p, values_fill = 0) else
      tibble::tibble(ecoregion = character(0), lt10 = numeric(0), c1020 = numeric(0), ge20 = numeric(0))
  esh <- if (!is.null(eco)) eco |> dplyr::select(ecoregion, n_agc, t2_stratum) else
    tibble::tibble(ecoregion = character(0), n_agc = numeric(0), t2_stratum = logical(0))
  d14 <- dce |> dplyr::mutate(ecoregion = ifelse(ecoregion == "All", "Indonesia", ecoregion)) |>
    dplyr::left_join(cvw |> dplyr::mutate(ecoregion = ifelse(ecoregion == "All", "Indonesia", ecoregion)), by = "ecoregion") |>
    dplyr::left_join(esh, by = "ecoregion") |>
    dplyr::mutate(beyond5_km2 = area_km2 * beyond_5km_pct / 100) |>
    dplyr::arrange(ecoregion == "Indonesia", dplyr::coalesce(eco_id, 0) == 0 & ecoregion != "Indonesia",
                   dplyr::desc(beyond5_km2))
  S14 <- d14 |>
    dplyr::transmute(Ecoregion = ecoregion, !!sprintf("Extent (%s)", KM2) := num(area_km2, 1),
                     `AGC samples (n)` = ifelse(ecoregion == "Indonesia", num(nrow(agc), 0), num(dplyr::coalesce(n_agc, 0), 0)),
                     `Within 1 km of a sample (%)` = num(within_1km_pct, 1),
                     `Beyond 5 km (%)` = num(beyond_5km_pct, 1),
                     !!sprintf("Beyond 5 km (%s)", KM2) := num(beyond5_km2, 1),
                     `Beyond 50 km (%)` = num(beyond_50km_pct, 1),
                     `Median distance class` = median_class,
                     `CV < 10 % (% of extent)` = num(lt10, 1), `CV 10-20 %` = num(c1020, 1), `CV >= 20 %` = num(ge20, 1))
  top <- d14 |> dplyr::filter(ecoregion != "Indonesia", dplyr::coalesce(eco_id, 0) > 0)
  nat14 <- d14 |> dplyr::filter(ecoregion == "Indonesia")
  if (nrow(nat14) == 1) {
    derive("extent_beyond5km_pct", nat14$beyond_5km_pct, "share of the extent > 5 km from an AGC sample")
    derive("extent_beyond50km_pct", nat14$beyond_50km_pct, "share of the extent > 50 km from an AGC sample")
    derive("extent_within1km_pct", nat14$within_1km_pct, "share of the extent within 1 km of an AGC sample")
    if (nrow(top) > 0) derive("top_eco_beyond5km_share_pct", 100 * top$beyond5_km2[1] / nat14$beyond5_km2,
                              sprintf("share of the national extent beyond 5 km that lies in %s", top$ecoregion[1]))
  }
  save_table("dist_cv",
             "Distance from the extent to the nearest AGC sample, and the between-year coefficient of variation (CV) of mapped AGC, by ecoregion (ordered by extent beyond 5 km).",
             S14, c("Distance: from each extent pixel to the nearest AGC sample of any year (P4_G03, 11 classes).",
                    "CV: SD / mean of the 8 annual AGC maps (2017-2024) per pixel; shares of the extent."))
}

# =============================================================================
# Write the combined Markdown file, the checks, and copies of the inputs
# =============================================================================
step("WRITE")
md_table <- function(d) {
  d <- as.data.frame(d); d[] <- lapply(d, function(v) gsub("\\|", "/", as.character(v)))
  c(paste0("| ", paste(names(d), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
    apply(d, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
md <- c("# Paper 4: tables", "",
        sprintf("Generated by P4_09_tables.R on %s from %s.", format(Sys.time(), "%Y-%m-%d %H:%M"), CFG$out_dir), "",
        "## Consistency checks", "", paste0("- ", CHECKS), "")
tab_order <- function(id) vapply(id, function(i) if (substr(i, 1, 1) == "S") 100 + as.numeric(substring(i, 2)) else
  as.numeric(i), numeric(1))
TABLES <- TABLES[order(tab_order(names(TABLES)))]
for (id in names(TABLES)) {
  t <- TABLES[[id]]
  md <- c(md, paste0("## ", t$title), "", sprintf("File: `%s`", t$file), "", md_table(t$data), "",
          if (length(t$notes) > 0) paste0("- ", t$notes), "")
}
writeLines(enc2utf8(md), file.path(TAB$dir, "P4_tables_all.md"), useBytes = TRUE)
writeLines(enc2utf8(CHECKS), file.path(TAB$dir, "P4_09_checks.txt"), useBytes = TRUE)
readr::write_csv(DERIVED, file.path(TAB$dir, "P4_09_derived_numbers.csv"))
cat(sprintf("Derived numbers for the text: %d (P4_09_derived_numbers.csv)\n", nrow(DERIVED)))
inp <- setdiff(used, "P4_samples.csv")         # P4_samples.csv is large and holds the field sample locations
ok_copy <- file.copy(out_file(inp), file.path(TAB$dir, "inputs"), overwrite = TRUE)
cat(sprintf("\n%d tables written to %s\n", length(TABLES), TAB$dir))
cat(sprintf("Copied %d source CSVs to %s\n", sum(ok_copy), file.path(TAB$dir, "inputs")))
n_cek <- sum(grepl("^\\[CEK", CHECKS))
cat(sprintf("Consistency checks: %d OK, %d CEK\n", sum(grepl("^\\[OK", CHECKS)), n_cek))
if (n_cek > 0) cat("Some checks failed: send the folder P4_tables_v3 as it is; do not edit the tables yet.\n")
cat("Zip the folder P4_tables_v3 and send it.\n")
cat("\nP4_09_tables.R COMPLETE\n")
