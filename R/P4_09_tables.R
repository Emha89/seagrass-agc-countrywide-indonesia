# =============================================================================
# P4_09_tables.R
# Paper 4 (revised design), added 7 Nov 2026: all tables of the manuscript and
# the Supplement, built from the R results (P4_01, P4_04, P4_05, P4_07) and the
# GEE exports. Nothing is recomputed from raw data: every number comes from a
# result CSV, so a table and its source can be checked line by line.
#
# Update 9 Nov 2026: tables renumbered in the order of their
# first citation, written to a NEW folder P4_tables_v2 (P4_tables is not touched).
# Other changes: the ACA geomorphic class removed from the extent is labelled
# "Terrestrial Reef Flat"; new Table S5 (presence model CV, P4_G09); Table 4 adds
# the share of pixel-based AGC per design; Table S3 adds the ICC of field AGC and
# the 100 m blocks surveyed in more than one year; Table S7 adds the RMSE range
# over the tuning grid; Table S9 adds the cover bias; P4_09_derived_numbers.csv
# keeps the numbers computed here that the text quotes (read by P4_10).
#
# Main tables (was -> now)
#   T1  Field data and extent by ecoregion                                (T1)
#   T2  Extent under the geomorphic mask options and the KKP 2025 area     (T5)
#   T3  AGC by ecoregion and leave-one-ecoregion-out bias                  (T4)
#   T4  Validation by design and method: RMSE, MAE, R2, MEC, bias, with the
#       share of extent and of pixel-based AGC that each design represents (T2)
#   T5  National AGC 2017-2024 by method (95 % CI)                         (T3)
# Supplement
#   S1  AGC samples by ecoregion and survey year                           (S1)
#   S2  Tier 1 by survey year and the survey composition                   (S9)
#   S3  Sampling structure and uncertainty parameters                      (S5)
#   S4  The 22 ROIs: field, final map and out-of-fold means                (S10)
#   S5  Presence model cross-validation at thresholds 0.5, 0.6, 0.7        (new)
#   S6  AGC samples inside the extent and the masks that removed the rest  (S8)
#   S7  Random forest settings (AGC and cover models)                      (S2)
#   S8  Tier 2 strata: class means and class areas by year                 (S4)
#   S9  Cover model validation by design (continuous and class accuracy)   (S3)
#   S10 Validation by ecoregion (leave-one-ecoregion-out)                  (S6)
#   S11 Distance to AGC samples and between-year CV by ecoregion           (S7)
# The code sections keep the order of the 7 Nov version; P4_tables_all.md lists
# the tables in manuscript order.
#
# Output (new folder, nothing else is touched): <out_dir>/P4_tables_v2/
#   P4_T1_...csv ... P4_S11_...csv  one file per table (UTF-8 with BOM, opens in Excel)
#   P4_tables_all.md                all tables with their notes, for drafting
#   P4_09_checks.txt                consistency checks between the source files
#   P4_09_derived_numbers.csv       numbers computed here and quoted in the text
#   P4_S4_roi_sites.csv             S4 site names: every location code of the samples
#                                   in each ROI (field label file), with counts and the
#                                   display name used in S4 (TAB$site_names, 8 Nov 2026)
#   inputs/                         copies of the source CSVs
# Rerunning P4_09 rewrites only the files in P4_tables_v2.
#
# Extra inputs of this version: <gee_dir>/P4_gee_pa_cv.csv and
# P4_gee_pa_cv_counts.csv (P4_G09; copy them from Drive P4v2_GEE to result/P4v2/gee),
# P4_tuning_AGC.csv and P4_tuning_COVER.csv (P4_02), P4_oof_cov_<design>.csv (P4_03).
# A missing input only skips the part that needs it.
#
# Number format: values with a 95 % range are written "mean [low, high]".
# For T4, S9 and S10 the range is the 2.5-97.5 percentile over the CV repeats
# (10 repeats; 3 seeds for the ecoregion design); for totals it is the 95 % CI.
# =============================================================================

# setwd("path/to/seagrass-agc-countrywide-indonesia/R")   # folder that holds these scripts
source("P4_00_config_functions.R")
cat("P4_09_tables.R | version 2026-11-09b\n")

# -----------------------------------------------------------------------------
# TABLE CONFIG -- edit here only
# -----------------------------------------------------------------------------
TAB <- list(
  dir        = out_file("P4_tables_v2"),        # new folder (9 Nov 2026)
  gee_dir    = if (!is.null(CFG$gee_dir)) CFG$gee_dir else out_file("gee"),
  trf_label  = "Terrestrial Reef Flat",       # ACA geomorphic class removed from the extent (code 16)
  kkp_km2    = 6601.56,                     # KKP 2025 national seagrass map (official area)
  turbid_eco = c("Sunda Shelf/Java Sea"),   # T2: option with Terrestrial Reef Flat removed only here
  years      = 2017:2024,
  tol_pct    = 0.1,                         # consistency checks: allowed difference (%)
  # 22 ROIs (Supplementary Table S4). Site name, ecoregion and centroid
  # come from the data (8 Nov 2026): the ROI polygons (CFG$roi_path, else
  # <out_dir>/gee/P4_roi_polygons.geojson, as P4_07) and the location names of the AGC
  # samples inside each ROI. master_raw_YYYY.csv does not keep the location column
  # (01_build_master_raw.R keeps only model columns), so the names are read from the
  # field label file that 01_build_master_raw.R reads (label_path there), joined by gee_id.
  roi22 = sprintf("ROI-%03d", c(95, 70, 272, 115, 83, 210, 192, 317, 153, 214, 320,
                                41, 149, 3, 156, 57, 96, 159, 148, 103, 53, 102)),
  location_path  = "../data/RF_AGC_summary_v4_clean_for_QGIS_merged.csv",   # field label file (gee_id, location; not shared)
  location_field = "location",   # if absent, the first column named like loc*/site/lokasi is used
  # Display names of the location codes for the manuscript (confirmed by the authors, 8 Nov 2026).
  # A code not listed here is written in title case with "_" as a space, and S4 checks report it.
  site_names = c(kepulauan_seribu = "Seribu Islands", bajo = "Labuan Bajo", tapteng = "Central Tapanuli",
                 komodo = "Komodo", manado = "Manado", ayau = "Ayau", bali = "Bali", rote = "Rote",
                 derawan = "Derawan", belitung = "Belitung", karimunjawa = "Karimunjawa", bintan = "Bintan"),
  roi_min_n  = 10                           # ROIs below this n are marked in S4
)

DESIGN_LEVELS <- c("pixel", "b100", "b1000", "cluster", "ecoregion")
DESIGN_NAME <- c(pixel = "Pixel", b100 = "100 m block", b1000 = "1 km block",
                 cluster = "Meadow cluster", ecoregion = "Ecoregion")
DESIGN_RANGE <- c(pixel = "\u2264 10 m", b100 = "10\u2013100 m", b1000 = "100 m\u20131 km",
                  cluster = "1\u201320 km", ecoregion = "> 20 km")
METHODS <- c("Pixel-based", "Tier 1", "Tier 2")
KM2 <- "km\u00b2"; GCM2 <- "gC m\u207b\u00b2"; R2 <- "R\u00b2"

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
  ifelse(is.na(x), "\u2013", formatC(round(x, d) + 0, format = "f", digits = d, big.mark = big))
}
share_lab <- function(p) ifelse(!is.na(p) & p > 0 & p < 0.05, "< 0.1", num(p, 1))
rng <- function(m, lo, hi, d = 1) ifelse(is.na(m), "\u2013",
                                         ifelse(is.na(lo) | is.na(hi), num(m, d),
                                                sprintf("%s [%s, %s]", num(m, d), num(lo, d), num(hi, d))))
pm  <- function(m, s, d = 1) ifelse(is.na(m), "\u2013", sprintf("%s \u00b1 %s", num(m, d), num(s, d)))
dist_lab <- function(m) ifelse(m < 1000, sprintf("%.0f m", m),
                               ifelse(m < 1e5, sprintf("%.1f km", m / 1000), sprintf("%.0f km", m / 1000)))

TABLES <- list()                               # name -> list(title, data, notes)
CHECKS <- character(0)
check <- function(label, ok, detail = "") {
  line <- sprintf("[%s] %s%s", if (isTRUE(ok)) "OK  " else "CEK ", label,
                  if (nzchar(detail)) paste0(": ", detail) else "")
  cat(line, "\n"); CHECKS <<- c(CHECKS, line)
}
close_pct <- function(a, b) abs(a - b) <= TAB$tol_pct / 100 * max(abs(b), 1e-9)
# Numbers computed in this script and quoted in the text (9 Nov 2026), written to
# P4_09_derived_numbers.csv and read by P4_10_text_numbers.R
DERIVED <- tibble::tibble(name = character(0), value = numeric(0), how = character(0))
derive <- function(name, value, how) {
  DERIVED <<- dplyr::bind_rows(DERIVED, tibble::tibble(name = name, value = as.numeric(value), how = how))
}

save_table <- function(id, file, title, d, notes = character(0)) {
  readr::write_excel_csv(d, file.path(TAB$dir, file), na = "")
  TABLES[[id]] <<- list(title = title, data = d, notes = notes, file = file)
  cat(sprintf("\n==== %s ====\n%s\n", id, title))
  print(as.data.frame(d), row.names = FALSE, right = FALSE)
  if (length(notes) > 0) cat(paste0("  ", notes, collapse = "\n"), "\n")
}
step <- function(txt) cat(sprintf("\n######## %s ########\n", txt))

# -----------------------------------------------------------------------------
# Inputs used by several tables
# -----------------------------------------------------------------------------
samples <- read_samples(); used <- union(used, "P4_samples.csv")
agc <- samples |> dplyr::filter(is_agc)
cov <- samples |> dplyr::filter(is_cov)
px_nat <- src("P4_pixel_national.csv")
nat3   <- src("P4_national_three_methods.csv")
eco    <- src("P4_ecoregion_summary.csv")
AREA <- if (!is.null(px_nat)) mean(px_nat$area_km2) else NA_real_
NAT_MEAN <- mean(agc$AGC_pred)
cat(sprintf("AGC samples %d, cover samples %d, extent %.2f %s, national field mean %.3f %s\n",
            nrow(agc), nrow(cov), AREA, KM2, NAT_MEAN, GCM2))

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
  check("Ecoregion-by-year areas add up to the national extent",
        all(close_pct(s$a, s$area_km2)), sprintf("largest difference %.3f %%", max(abs(s$a - s$area_km2) / s$area_km2 * 100)))
}
if (!is.null(nat3)) {
  t1p <- nat3 |> dplyr::filter(method == "Tier 1 (pooled mean)")
  check("Tier 1 pooled total = national field mean x extent",
        close_pct(t1p$total_tC[1], NAT_MEAN * AREA), sprintf("%.0f vs %.0f t C", t1p$total_tC[1], NAT_MEAN * AREA))
  pxn <- nat3 |> dplyr::filter(method == "Pixel-based")
  check("Pixel-based totals match P4_pixel_national", all(close_pct(sort(pxn$total_tC), sort(px_nat$total_tC))))
}
t2n <- src("P4_tier2_national_map.csv")
if (!is.null(t2n) && !is.null(nat3)) {
  t2a <- nat3 |> dplyr::filter(method == "Tier 2 (pooled class means)") |>
    dplyr::inner_join(t2n |> dplyr::select(year, t2 = total_tC), by = "year")
  check("Tier 2 totals in P4_tier2_national_map match P4_national_three_methods",
        all(close_pct(t2a$t2, t2a$total_tC)), sprintf("largest difference %.3f %%",
                                                     max(abs(t2a$t2 - t2a$total_tC) / t2a$total_tC * 100)))
}
cls <- src("gee_cover_class_area_by_year.csv")
if (!is.null(cls)) {
  s <- cls |> dplyr::group_by(year) |> dplyr::summarise(a = sum(area_km2), .groups = "drop")
  check("Cover class areas add up to the extent every year", all(close_pct(s$a, AREA)),
        sprintf("%.2f-%.2f %s", min(s$a), max(s$a), KM2))
}
dcd <- src("P4_distance_coverage_design.csv")
if (!is.null(dcd)) check("Distance shares by design add up to 100 %", close_pct(sum(dcd$area_pct), 100),
                         sprintf("%.3f %%", sum(dcd$area_pct)))
t1g <- src("P4_tier1_national_gee.csv")
if (!is.null(t1g)) {
  pooled <- t1g |> dplyr::filter(year == "pooled")
  check("Tier 1 pooled mean in P4_tier1_national_gee = national field mean",
        close_pct(pooled$mean_gC_m2, NAT_MEAN), sprintf("%.3f vs %.3f", pooled$mean_gC_m2, NAT_MEAN))
  check("Tier 1 pooled total in P4_tier1_national_gee uses the final extent",
        close_pct(pooled$total_tC, NAT_MEAN * AREA),
        sprintf("%.0f vs %.0f t C (a smaller value means the file is older than the add-on run)",
                pooled$total_tC, NAT_MEAN * AREA))
}

# =============================================================================
# T1: field data and extent by ecoregion
# =============================================================================
step("T1")
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
  ncov <- cov |> dplyr::count(ecoregion, name = "n_cov")
  t1 <- eco |>
    dplyr::left_join(yrs, by = "ecoregion") |> dplyr::left_join(clu, by = "ecoregion") |>
    dplyr::left_join(ncov, by = "ecoregion") |>
    dplyr::arrange(eco_id == 0, dplyr::desc(area_km2))
  tot_row <- tibble::tibble(
    Ecoregion = "Indonesia", a = sum(t1$area_km2), sa = 100, n = nrow(agc), sn = 100,
    years = paste(sort(unique(agc$year)), collapse = ", "),
    clusters = dplyr::n_distinct(agc$g_cluster[grepl("^c", agc$g_cluster)]),
    noise = sum(!grepl("^c", agc$g_cluster)),
    fm = mean(agc$AGC_pred), fs = sd(agc$AGC_pred), nc = nrow(cov))
  T1 <- dplyr::bind_rows(
    t1 |> dplyr::transmute(Ecoregion = ecoregion, a = area_km2, sa = share_area_pct, n = n_agc,
                           sn = share_samples_pct, years, clusters, noise,
                           fm = field_mean_gC_m2, fs = field_sd_gC_m2, nc = n_cov),
    tot_row) |>
    dplyr::transmute(Ecoregion,
                     !!sprintf("Extent (%s)", KM2) := num(a, 1),
                     `Extent share (%)` = num(sa, 1),
                     `AGC samples (n)` = num(dplyr::coalesce(n, 0), 0),
                     `Sample share (%)` = num(dplyr::coalesce(sn, 0), 1),
                     `Survey years` = ifelse(is.na(years) | years == "", "\u2013", years),
                     `Meadow clusters (n)` = num(dplyr::coalesce(clusters, 0L), 0),
                     `Unclustered samples (n)` = num(dplyr::coalesce(noise, 0L), 0),
                     !!sprintf("Field AGC, mean \u00b1 SD (%s)", GCM2) := pm(fm, fs, 1),
                     `Cover samples (n)` = num(dplyr::coalesce(nc, 0L), 0))
  save_table("T1", "P4_T1_ecoregion_data.csv",
             "Table 1. Seagrass extent, field samples and field AGC by marine ecoregion (MEOW).",
             T1, c(sprintf("Extent: final map extent (core + add-on), %s %s.", num(AREA, 2), KM2),
                   sprintf("Meadow clusters: DBSCAN clusters of AGC samples (eps %s m, minPts %s)%s",
                           CFG$dbscan_eps, CFG$dbscan_minPts,
                           if (sum(t1$clusters, na.rm = TRUE) == tot_row$clusters) "; every cluster lies in one ecoregion." else
                             "; a cluster that crosses an ecoregion boundary is counted in each ecoregion, so the column does not add up to the national count."),
                   "Unclustered samples: DBSCAN noise points, each treated as its own group in the cluster design.",
                   "Cover samples: samples with seagrass cover (tSPC) used for the cover model and the Tier 2 strata."))
}

# =============================================================================
# T4: validation by design and method
# =============================================================================
step("T4")
ms <- src("P4_metrics_summary.csv")
if (!is.null(ms)) {
  sh <- if (!is.null(dcd)) setNames(dcd$area_pct, dcd$design) else setNames(rep(NA_real_, 5), DESIGN_LEVELS)
  sa <- if (!is.null(dcd) && "agc_pct_mean" %in% names(dcd)) setNames(dcd$agc_pct_mean, dcd$design) else
    setNames(rep(NA_real_, 5), DESIGN_LEVELS)
  T4 <- ms |>
    dplyr::mutate(design = factor(design, levels = DESIGN_LEVELS), method = factor(method, levels = METHODS)) |>
    dplyr::arrange(design, method) |>
    dplyr::transmute(Design = ifelse(method == "Pixel-based", sprintf("%s (%s)", DESIGN_NAME[as.character(design)],
                                                                        DESIGN_RANGE[as.character(design)]), ""),
                     `Median distance to training` = ifelse(method == "Pixel-based", dist_lab(median_dist_m), ""),
                     `Share of extent (%)` = ifelse(method == "Pixel-based", share_lab(sh[as.character(design)]), ""),
                     `Share of pixel-based AGC (%)` = ifelse(method == "Pixel-based", share_lab(sa[as.character(design)]), ""),
                     Method = as.character(method),
                     !!sprintf("RMSE (%s)", GCM2) := rng(RMSE_mean, RMSE_q025, RMSE_q975, 1),
                     !!sprintf("MAE (%s)", GCM2) := rng(MAE_mean, MAE_q025, MAE_q975, 1),
                     !!R2 := rng(R2_mean, R2_q025, R2_q975, 3),
                     MEC = rng(MEC_mean, MEC_q025, MEC_q975, 3),
                     `Bias (%)` = rng(bias_pct_mean, bias_pct_q025, bias_pct_q975, 1))
  nrep <- ms |> dplyr::distinct(design, n_repeats)
  save_table("T4", "P4_T4_validation.csv",
             "Table 4. Out-of-fold accuracy of AGC by validation design and estimation method.",
             T4, c(sprintf("n = %s AGC samples in every design; mean over %s repeats of %s-fold CV (ecoregion design: %s seeds, one fold per ecoregion), with the 2.5-97.5 percentile range in brackets.",
                           num(unique(ms$n)[1], 0), nrep$n_repeats[nrep$design == "pixel"], CFG$k_folds,
                           nrep$n_repeats[nrep$design == "ecoregion"]),
                   "Distance range: distance from an extent pixel to the nearest AGC sample that this design represents; share of extent: part of the national extent in that range (Table S11, Fig. 5); share of pixel-based AGC: part of the mean annual pixel-based total in that range.",
                   "Median distance to training: median distance from a test sample to the nearest training sample in that design.",
                   "Tier 1: national field mean of the training folds; Tier 2: training-fold mean of the predicted cover class.",
                   "MEC: model efficiency coefficient, 1 - SSE / sum((obs - mean(obs))^2); bias: (mean pred - mean obs) / mean obs x 100."))
}

# =============================================================================
# T5: national AGC by year and method
# =============================================================================
step("T5")
if (!is.null(nat3)) {
  w <- function(m) nat3 |> dplyr::filter(method == m) |> dplyr::select(year, total_kt, ci_low_kt, ci_up_kt)
  px <- w("Pixel-based"); t1y <- w("Tier 1 (survey-year mean)"); t1p <- w("Tier 1 (pooled mean)")
  t2 <- w("Tier 2 (pooled class means)")
  ny <- agc |> dplyr::count(year, name = "n")
  loyo <- if (!is.null(t2n)) t2n |> dplyr::transmute(year, lmin = loyo_min_tC / 1e3, lmax = loyo_max_tC / 1e3) else
    tibble::tibble(year = integer(0), lmin = numeric(0), lmax = numeric(0))
  neff <- if (!is.null(t1g)) t1g |> dplyr::filter(year != "pooled") |>
    dplyr::transmute(year = as.integer(year), n_eff) else tibble::tibble(year = integer(0), n_eff = numeric(0))
  b <- tibble::tibble(year = TAB$years) |>
    dplyr::left_join(ny, by = "year") |> dplyr::left_join(neff, by = "year") |>
    dplyr::left_join(px |> dplyr::rename(px = total_kt, pxl = ci_low_kt, pxu = ci_up_kt), by = "year") |>
    dplyr::left_join(t1y |> dplyr::rename(ty = total_kt, tyl = ci_low_kt, tyu = ci_up_kt), by = "year") |>
    dplyr::left_join(t1p |> dplyr::rename(tp = total_kt, tpl = ci_low_kt, tpu = ci_up_kt), by = "year") |>
    dplyr::left_join(t2 |> dplyr::rename(t2 = total_kt, t2l = ci_low_kt, t2u = ci_up_kt), by = "year") |>
    dplyr::left_join(loyo, by = "year")
  T5 <- b |> dplyr::transmute(
    Year = as.character(year),
    `AGC samples (n)` = ifelse(is.na(n), "0", num(n, 0)),
    `Pixel-based` = paste0(rng(px, pxl, pxu, 1), ifelse(year %in% ny$year, "", " \u2020")),
    `Tier 1, survey-year mean` = ifelse(is.na(ty), "\u2013",
                                        ifelse(is.na(tyl), paste0(num(ty, 1), " \u2021"), rng(ty, tyl, tyu, 1))),
    `Tier 1, pooled mean` = rng(tp, tpl, tpu, 1),
    `Tier 2, pooled class means` = rng(t2, t2l, t2u, 1),
    `Tier 2, leave-one-year-out range` = ifelse(is.na(lmin), "\u2013", sprintf("%s\u2013%s", num(lmin, 1), num(lmax, 1))))
  sm <- function(x, f) num(f(x, na.rm = TRUE), 1)
  T5 <- dplyr::bind_rows(T5,
    tibble::tibble(Year = "Mean 2017\u20132024", `AGC samples (n)` = num(nrow(agc), 0),
                   `Pixel-based` = sm(b$px, mean), `Tier 1, survey-year mean` = "\u2013",
                   `Tier 1, pooled mean` = rng(b$tp[1], b$tpl[1], b$tpu[1], 1),
                   `Tier 2, pooled class means` = sm(b$t2, mean), `Tier 2, leave-one-year-out range` = "\u2013"),
    tibble::tibble(Year = "Range of annual values", `AGC samples (n)` = "",
                   `Pixel-based` = sprintf("%s\u2013%s", sm(b$px, min), sm(b$px, max)),
                   `Tier 1, survey-year mean` = sprintf("%s\u2013%s", sm(b$ty, min), sm(b$ty, max)),
                   `Tier 1, pooled mean` = "constant",
                   `Tier 2, pooled class means` = sprintf("%s\u2013%s", sm(b$t2, min), sm(b$t2, max)),
                   `Tier 2, leave-one-year-out range` = sprintf("%s\u2013%s", sm(b$lmin, min), sm(b$lmax, max))))
  save_table("T5", "P4_T5_national_totals.csv",
             "Table 5. National seagrass AGC (kt C) by year and estimation method, with 95 % confidence intervals.",
             T5, c(sprintf("All methods use the same fixed extent of %s %s (2017-2024).", num(AREA, 2), KM2),
                   "Pixel-based CI: 1.96 x cluster-robust SE of the mean out-of-fold residual (cluster design) x extent (Table S3).",
                   "Tier 1 and Tier 2 CIs use n / DEFF as the effective sample size (Table S3).",
                   "\u2020 No field data in 2024: the CI assumes the 2017-2023 bias SE.",
                   sprintf("\u2021 Effective sample size below 2 (n = %s): no CI.",
                           paste(b$n[!is.na(b$ty) & is.na(b$tyl)], collapse = ", ")),
                   "Leave-one-year-out range: Tier 2 total with class means recomputed without one survey year at a time (sensitivity to survey composition)."))
}

# =============================================================================
# T3: AGC by ecoregion and leave-one-ecoregion-out bias
# =============================================================================
step("T3")
lolo <- src("P4_lolo_region_means.csv")
if (!is.null(eco)) {
  mm <- if (!is.null(eby)) eby |> dplyr::group_by(ecoregion) |>
    dplyr::summarise(pmin = min(px_total_tC), pmax = max(px_total_tC), .groups = "drop") else
      tibble::tibble(ecoregion = character(0), pmin = numeric(0), pmax = numeric(0))
  lb <- if (!is.null(lolo)) lolo |> dplyr::filter(reported) |>
    dplyr::select(ecoregion, bp = bias_pct_pixel, b1 = bias_pct_tier1, b2 = bias_pct_tier2) else
      tibble::tibble(ecoregion = character(0), bp = numeric(0), b1 = numeric(0), b2 = numeric(0))
  t4 <- eco |> dplyr::left_join(mm, by = "ecoregion") |> dplyr::left_join(lb, by = "ecoregion") |>
    dplyr::arrange(eco_id == 0, dplyr::desc(area_km2))
  T3 <- dplyr::bind_rows(
    t4 |> dplyr::transmute(
      Ecoregion = ecoregion,
      !!sprintf("Extent (%s)", KM2) := num(area_km2, 1),
      `AGC samples (n)` = num(dplyr::coalesce(n_agc, 0), 0),
      !!sprintf("Field mean (%s)", GCM2) := num(field_mean_gC_m2, 1),
      !!sprintf("Map mean (%s)", GCM2) := num(px_mean_gC_m2, 1),
      `Pixel-based, mean [min, max] (kt C)` = rng(px_total_tC / 1e3, pmin / 1e3, pmax / 1e3, 2),
      `Tier 1, ecoregion mean [95 % CI] (kt C)` = ifelse(is.na(t1_eco_total_tC), "\u2013",
        rng(t1_eco_total_tC / 1e3, (t1_eco_total_tC - t1_eco_ci_half_tC) / 1e3,
            (t1_eco_total_tC + t1_eco_ci_half_tC) / 1e3, 2)),
      `Tier 1, national mean (kt C)` = num(t1_national_factor_tC / 1e3, 2),
      `Tier 2 (kt C)` = num(t2_total_tC / 1e3, 2),
      `LOLO bias, pixel-based (%)` = num(bp, 1),
      `LOLO bias, Tier 1 (%)` = num(b1, 1),
      `LOLO bias, Tier 2 (%)` = num(b2, 1)),
    tibble::tibble(Ecoregion = "Indonesia",
                   !!sprintf("Extent (%s)", KM2) := num(sum(eco$area_km2), 1),
                   `AGC samples (n)` = num(sum(eco$n_agc, na.rm = TRUE), 0),
                   !!sprintf("Field mean (%s)", GCM2) := num(NAT_MEAN, 1),
                   !!sprintf("Map mean (%s)", GCM2) := num(sum(eco$px_total_tC) / sum(eco$area_km2), 1),
                   `Pixel-based, mean [min, max] (kt C)` = if (!is.null(px_nat))
                     rng(mean(px_nat$total_tC) / 1e3, min(px_nat$total_tC) / 1e3, max(px_nat$total_tC) / 1e3, 2) else "\u2013",
                   `Tier 1, ecoregion mean [95 % CI] (kt C)` = "\u2013",
                   `Tier 1, national mean (kt C)` = num(sum(eco$t1_national_factor_tC) / 1e3, 2),
                   `Tier 2 (kt C)` = num(sum(eco$t2_total_tC) / 1e3, 2),
                   `LOLO bias, pixel-based (%)` = "", `LOLO bias, Tier 1 (%)` = "", `LOLO bias, Tier 2 (%)` = ""))
  save_table("T3", "P4_T3_ecoregion_totals.csv",
             "Table 3. Seagrass AGC by marine ecoregion: means, totals by method, and bias when the ecoregion is left out.",
             T3, c("Pixel-based: mean of the annual totals 2017-2024, with the lowest and highest annual total; map mean = mean total / extent.",
                   sprintf("Tier 1, ecoregion mean: ecoregion field mean x extent, only where n >= %d; the 95 %% CI uses the national DEFF (approximation).", CFG$min_n_eco),
                   "Tier 1, national mean: national field mean x ecoregion extent. Tier 2: pooled national class means x class areas, mean of 2017-2024.",
                   sprintf("LOLO bias: (predicted - field mean) / field mean x 100 for the samples of the ecoregion when it is left out of training (ecoregion design); only where n >= %d.", CFG$min_n_eco)))
}

# =============================================================================
# T2: extent under the geomorphic mask options and the KKP 2025 area
# =============================================================================
step("T2")
esn <- src("P4_extent_sensitivity_national.csv"); ese <- src("P4_extent_sensitivity_ecoregion.csv")
if (!is.null(esn) && !is.null(ese) && !is.null(px_nat)) {
  core <- esn$area_km2_at_10m[esn$mask == "current"]
  f30  <- core / sum(ese$current_km2)              # 30 m area -> 10 m map area (as P4_07 STEP 12)
  only16_est <- f30 * sum(ese$only16_km2)
  regional <- f30 * (sum(ese$only16_km2) + sum((ese$nogeo_km2 - ese$only16_km2)[!ese$ecoregion %in% TAB$turbid_eco]))
  nogeo <- f30 * sum(ese$nogeo_km2)
  dens <- px_nat$mean_gC_m2                         # annual national mean density of the map
  opt <- tibble::tibble(
    Option = c(sprintf("Core extent: ACA geomorphic data present, %s removed", TAB$trf_label),
               sprintf("Final extent: core + add-on (only %s removed)", TAB$trf_label),
               sprintf("%s removed only in %s", TAB$trf_label, paste(TAB$turbid_eco, collapse = ", ")),
               "No geomorphic mask",
               "KKP 2025 national seagrass map"),
    Basis = c("mapped (10 m)", "mapped (10 m)", "estimated from 30 m", "estimated from 30 m", "official map"),
    km2 = c(core, AREA, regional, nogeo, TAB$kkp_km2))
  T2 <- opt |> dplyr::transmute(
    Option, Basis,
    !!sprintf("Extent (%s)", KM2) := num(km2, 0),
    `Share of KKP 2025 (%)` = num(100 * km2 / TAB$kkp_km2, 1),
    `Tier 1, pooled mean (kt C)` = num(NAT_MEAN * km2 / 1e3, 1),
    `Map density x extent, range 2017-2024 (kt C)` = sprintf("%s\u2013%s", num(min(dens) * km2 / 1e3, 1),
                                                             num(max(dens) * km2 / 1e3, 1)))
  derive("extent_core_km2", core, "core extent (P4_extent_sensitivity_national, current, at 10 m)")
  derive("extent_trf_regional_km2", regional, "Terrestrial Reef Flat removed only in TAB$turbid_eco, scaled to 10 m")
  derive("extent_nogeo_km2", nogeo, "no geomorphic mask, scaled to 10 m")
  derive("tier1_kkp_tC", NAT_MEAN * TAB$kkp_km2, "national field mean x KKP 2025 area")
  derive("map_kkp_min_tC", min(dens) * TAB$kkp_km2, "lowest annual national map density x KKP 2025 area")
  derive("map_kkp_max_tC", max(dens) * TAB$kkp_km2, "highest annual national map density x KKP 2025 area")
  save_table("T2", "P4_T2_extent_options_kkp.csv",
             "Table 2. Seagrass extent under the geomorphic mask options, compared with the KKP 2025 national map, and the national AGC each extent implies.",
             T2, c("All options share MaxProb >= 0.7 and ACA depth 0-5 m; they differ only in the ACA geomorphic mask.",
                   sprintf("Estimated rows: areas of the 30 m sensitivity run (P4_G05) scaled to the 10 m map by the core ratio %.4f. The same scaling gives %s %s for the final extent, %.1f %% from the mapped %s %s.",
                           f30, num(only16_est, 0), KM2, 100 * (only16_est - AREA) / AREA, num(AREA, 0), KM2),
                   sprintf("Tier 1, pooled mean: national field mean (%.2f %s) x extent.", NAT_MEAN, GCM2),
                   sprintf("Map density x extent: annual national mean density of the AGC map (%.1f-%.1f %s) x extent; for the final extent this equals the pixel-based totals.",
                           min(dens), max(dens), GCM2)))
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
  save_table("S1", "P4_S1_samples_ecoregion_year.csv",
             "Table S1. AGC field samples by marine ecoregion and survey year.", S1,
             "No AGC samples were collected in 2024.")
}

# =============================================================================
# S7: random forest settings
# =============================================================================
step("S7")
bp <- src("P4_best_params.csv")
# Tuning grid (9 Nov 2026): RMSE range over all settings, from the P4_02 grid files
tg <- dplyr::bind_rows(lapply(c("P4_tuning_AGC.csv", "P4_tuning_COVER.csv"), function(f) if (have(f)) src(f) else NULL))
if (!is.null(bp)) {
  if (nrow(tg) > 0) {
    tgr <- tg |> dplyr::mutate(model = ifelse(toupper(model) == "AGC", "AGC", "COVER")) |> dplyr::group_by(model) |>
      dplyr::summarise(n_set = dplyr::n(), rmse_min = min(RMSE), rmse_max = max(RMSE),
                       mec_min = min(MEC), mec_max = max(MEC), .groups = "drop")
    for (i in seq_len(nrow(tgr))) {
      derive(sprintf("tuning_%s_n_settings", tgr$model[i]), tgr$n_set[i], "rows of the P4_02 tuning grid")
      derive(sprintf("tuning_%s_rmse_range", tgr$model[i]), tgr$rmse_max[i] - tgr$rmse_min[i], "max - min RMSE over the grid")
    }
    bp <- bp |> dplyr::mutate(key = ifelse(toupper(model) == "AGC", "AGC", "COVER")) |>
      dplyr::left_join(tgr, by = c("key" = "model"))
  } else bp <- bp |> dplyr::mutate(n_set = NA_real_, rmse_min = NA_real_, rmse_max = NA_real_)
  S7 <- bp |> dplyr::transmute(
    Model = ifelse(model == "AGC", "AGC (gC m\u207b\u00b2)", "Seagrass cover (%)"),
    `Training samples (n)` = num(n, 0), Trees = num(num_trees, 0),
    `Variables per split (mtry)` = num(mtry, 0), `Minimum node size` = num(min_node_size, 0),
    `Sample fraction` = num(sample_fraction, 1),
    `Tuning CV RMSE` = num(RMSE, 2), `Tuning CV MAE` = num(MAE, 2),
    !!paste("Tuning CV", R2) := num(R2, 3), `Tuning CV MEC` = num(MEC, 3), `Tuning CV bias (%)` = num(bias_pct, 2),
    `Settings in the grid (n)` = num(n_set, 0),
    `RMSE range over the grid` = ifelse(is.na(rmse_min), "\u2013", sprintf("%s\u2013%s", num(rmse_min, 2), num(rmse_max, 2))))
  save_table("S7", "P4_S7_model_settings.csv",
             "Table S7. Random forest settings of the AGC and cover models (64 GSE bands as predictors).", S7,
             c(sprintf("Settings chosen by %d-fold CV grouped by 10 m pixel over mtry {%s}, minimum node size {%s} and sample fraction {%s}.",
                       CFG$tune_folds, paste(CFG$grid_mtry, collapse = ", "), paste(CFG$grid_node, collapse = ", "),
                       paste(CFG$grid_frac, collapse = ", ")),
               "The same settings are used in GEE (smileRandomForest: numberOfTrees, variablesPerSplit, minLeafPopulation, bagFraction).",
               "Cover model: trained on samples with cover > 0."))
}

# =============================================================================
# S9: cover model validation
# =============================================================================
step("S9")
cm <- src("P4_cover_metrics.csv"); ca <- src("P4_cover_class_accuracy.csv")
# Bias of the cover model (9 Nov 2026): P4_cover_metrics.csv has no bias column, so it is
# computed here from the out-of-fold cover files exactly as P4_04 STEP 6 computes the other
# metrics (cover > 0, metrics per repeat, mean over repeats). RMSE is recomputed as a check.
cov_bias <- purrr::map_dfr(DESIGN_LEVELS, function(d) {
  f <- sprintf("P4_oof_cov_%s.csv", d)
  if (!have(f)) return(NULL)
  # read directly (large file, not copied to inputs/)
  readr::read_csv(out_file(f), show_col_types = FALSE) |> dplyr::filter(obs > 0) |> dplyr::group_by(rep) |>
    dplyr::summarise(reg_metrics(obs, pred), .groups = "drop") |>
    dplyr::summarise(design = d, bias_pct = mean(bias_pct), RMSE_chk = mean(RMSE))
})
if (!is.null(cm) && nrow(cov_bias) > 0) {
  k <- dplyr::inner_join(cm, cov_bias, by = "design")
  check("Cover RMSE recomputed from P4_oof_cov_* = P4_cover_metrics", all(close_pct(k$RMSE_chk, k$RMSE)),
        sprintf("largest difference %.4f %%", max(abs(k$RMSE_chk - k$RMSE) / k$RMSE * 100)))
}
if (!is.null(cm) && !is.null(ca)) {
  if (nrow(cov_bias) > 0) cm <- cm |> dplyr::left_join(cov_bias |> dplyr::select(design, bias_pct), by = "design") else
    cm$bias_pct <- NA_real_
  S9 <- cm |> dplyr::left_join(ca, by = "design") |>
    dplyr::mutate(design = factor(design, levels = DESIGN_LEVELS)) |> dplyr::arrange(design) |>
    dplyr::transmute(Design = DESIGN_NAME[as.character(design)], `Samples (n)` = num(n, 0),
                     `RMSE (% cover)` = num(RMSE, 1), `MAE (% cover)` = num(MAE, 1),
                     !!R2 := num(R2, 3), MEC = num(MEC, 3), `Bias (%)` = num(bias_pct, 1),
                     `Overall accuracy` = num(overall_accuracy, 3), Kappa = num(kappa, 3),
                     `Producer's accuracy by class` = producer_acc, `User's accuracy by class` = user_acc)
  save_table("S9", "P4_S9_cover_validation.csv",
             "Table S9. Out-of-fold accuracy of the cover model by validation design: continuous cover and the three cover classes (KepMenLH 200/2004).",
             S9, c("Continuous metrics: mean over the CV repeats (cover > 0); class accuracy: pooled over the repeats. Classes: < 30 %, 30-59.9 %, >= 60 % cover.",
                   "The cover map sets only the class areas of Tier 2; the class mean AGC comes from field samples (Table S8)."))
}

# =============================================================================
# S8: Tier 2 strata
# =============================================================================
step("S8")
cmm <- src("P4_tier2_class_means_map.csv"); lym <- src("P4_tier2_loyo_map.csv")
if (!is.null(cmm) && !is.null(cls)) {
  lr <- if (!is.null(lym)) lym |> dplyr::group_by(cov_class) |>
    dplyr::summarise(lo = min(mean_gC_m2), hi = max(mean_gC_m2), .groups = "drop") else
      tibble::tibble(cov_class = character(0), lo = numeric(0), hi = numeric(0))
  ar <- cls |> dplyr::mutate(year = as.character(year)) |>
    tidyr::pivot_wider(names_from = year, values_from = area_km2)
  S8 <- cmm |> dplyr::left_join(lr, by = "cov_class") |>
    dplyr::left_join(ar, by = c("cov_class" = "class")) |>
    dplyr::mutate(cov_class = factor(cov_class, levels = CFG$cover_labels)) |> dplyr::arrange(cov_class)
  yc <- intersect(as.character(TAB$years), names(S8))
  S8 <- S8 |> dplyr::transmute(
    `Cover class (%)` = as.character(cov_class), `AGC samples (n)` = num(n, 0),
    !!sprintf("Mean \u00b1 SD (%s)", GCM2) := pm(mean_gC_m2, sd_gC_m2, 1),
    `Effective n (n / DEFF)` = num(n_eff, 0),
    !!sprintf("Leave-one-year-out mean range (%s)", GCM2) := ifelse(is.na(lo), "\u2013",
                                                                     sprintf("%s\u2013%s", num(lo, 1), num(hi, 1))),
    dplyr::across(dplyr::all_of(yc), ~ num(.x, 1), .names = "Area {.col} (km\u00b2)"))
  save_table("S8", "P4_S8_tier2_strata.csv",
             "Table S8. Tier 2 strata: field AGC by mapped cover class (pooled over survey years) and the class areas of each year.",
             S8, c("Each AGC sample takes the class of the GEE cover map at its pixel in its survey year.",
                   "Class means are pooled over survey years, so annual change in Tier 2 comes only from change in class area."))
}

# =============================================================================
# S3: uncertainty parameters
# =============================================================================
step("S3")
up <- src("P4_uncertainty_params.csv"); csm <- src("P4_cluster_summary.csv")
if (!is.null(up)) {
  lab <- c(ICC = "Intra-cluster correlation of out-of-fold residuals",
           adjusted_mean_cluster_size = "Adjusted mean cluster size (n0)",
           DEFF = "Design effect, 1 + (n0 - 1) x ICC", N_agc = "AGC samples", N_eff = "Effective sample size, n / DEFF",
           SD_model_gC_m2 = sprintf("SD_model: RMSE of the cluster design (%s)", GCM2),
           SD_data_gC_m2 = sprintf("SD_data: RMS of (AGC_up - AGC_low) / 3.92 (%s)", GCM2),
           SD_total_gC_m2 = sprintf("SD_total: sqrt(SD_model\u00b2 + SD_data\u00b2) (%s)", GCM2),
           mean_resid_gC_m2 = sprintf("Mean out-of-fold residual, obs - pred, cluster design (%s)", GCM2),
           se_mean_resid_gC_m2 = sprintf("Cluster-robust SE of the mean residual (%s)", GCM2),
           bias_halfwidth_national_tC = "Pixel-based 95 % CI half-width, 1.96 x SE x extent (t C)")
  dg <- c(ICC = 3, adjusted_mean_cluster_size = 1, DEFF = 2, N_agc = 0, N_eff = 0, SD_model_gC_m2 = 2,
          SD_data_gC_m2 = 2, SD_total_gC_m2 = 2, mean_resid_gC_m2 = 2, se_mean_resid_gC_m2 = 2,
          bias_halfwidth_national_tC = 0)
  S3 <- up |> dplyr::filter(parameter %in% names(lab)) |>
    dplyr::transmute(Parameter = unname(lab[parameter]),
                     Value = mapply(function(v, p) num(v, dg[[p]]), value, parameter))
  if ("bias_halfwidth_national_tC" %in% up$parameter && !is.na(AREA)) {
    hw <- 1.96 * up$value[up$parameter == "se_mean_resid_gC_m2"] * AREA
    S3$Value[S3$Parameter == lab[["bias_halfwidth_national_tC"]]] <- num(hw, 0)
  }
  if (!is.null(csm)) S3 <- dplyr::bind_rows(
    tibble::tibble(Parameter = c("DBSCAN eps (m) / minPts", "Meadow clusters (n)", "AGC samples in clusters (n)",
                                 "Unclustered AGC samples (n)", "Largest cluster (n samples)"),
                   Value = c(sprintf("%s / %s", num(csm$eps_m, 0), csm$minPts), num(csm$n_clusters, 0),
                             num(csm$n_clustered, 0), num(csm$n_noise, 0), num(csm$largest_cluster, 0))), S3)
  # Sensitivity and sampling structure (9 Nov 2026): ICC of the field AGC values themselves
  # (clustered samples, same one-way ANOVA as P4_05) and the 100 m blocks with AGC samples
  # from more than one survey year (repeat surveys)
  clu <- agc |> dplyr::filter(grepl("^c", g_cluster))
  ic <- icc_oneway(clu$AGC_pred, clu$g_cluster)
  b100 <- agc |> dplyr::group_by(g_b100) |> dplyr::summarise(ny = dplyr::n_distinct(year), .groups = "drop")
  derive("icc_field_agc", ic$icc, "one-way ANOVA ICC of AGC_pred by g_cluster, clustered AGC samples (P4_samples.csv)")
  derive("deff_field_agc", ic$deff, "1 + (n0 - 1) x ICC of field AGC")
  derive("blocks100_with_agc", nrow(b100), "100 m blocks (g_b100) with AGC samples")
  derive("blocks100_multi_year", sum(b100$ny > 1), "100 m blocks with AGC samples from more than one survey year")
  S3 <- dplyr::bind_rows(S3, tibble::tibble(
    Parameter = c("ICC of field AGC (sensitivity)", "Design effect with the ICC of field AGC",
                  "100 m blocks with AGC samples (n)", "100 m blocks with AGC samples from more than one survey year (n)"),
    Value = c(num(ic$icc, 3), num(ic$deff, 2), num(nrow(b100), 0), num(sum(b100$ny > 1), 0))))
  if ("se_mean_resid_gC_m2" %in% up$parameter && !is.na(AREA))
    derive("pixel_ci_half_tC", 1.96 * up$value[up$parameter == "se_mean_resid_gC_m2"] * AREA,
           "1.96 x cluster-robust SE of the mean residual x final extent")
  save_table("S3", "P4_S3_uncertainty_parameters.csv",
             "Table S3. Clustering of the AGC samples, the parameters of the confidence intervals, and the repeat surveys.", S3,
             c(sprintf("CI half-width recomputed here with the final extent (%s %s); P4_05 used CFG$extent_km2.", num(AREA, 2), KM2),
               "The mean residual is not different from zero (|mean| < 1.96 x SE), so the totals are not bias-corrected."))
}

# =============================================================================
# S10: validation by ecoregion (leave-one-ecoregion-out)
# =============================================================================
step("S10")
me <- src("P4_metrics_ecoregion.csv"); mew <- src("P4_metrics_ecoregion_equal_weight.csv")
if (!is.null(me)) {
  ordr <- me |> dplyr::distinct(ecoregion, n) |> dplyr::arrange(dplyr::desc(n)) |> dplyr::pull(ecoregion)
  S10 <- me |> dplyr::mutate(ecoregion = factor(ecoregion, levels = ordr), method = factor(method, levels = METHODS)) |>
    dplyr::arrange(ecoregion, method) |>
    dplyr::transmute(Ecoregion = ifelse(method == "Pixel-based",
                                        paste0(as.character(ecoregion), ifelse(reported, "", " *")), ""),
                     `AGC samples (n)` = ifelse(method == "Pixel-based", num(n, 0), ""),
                     `Median distance to training (km)` = ifelse(method == "Pixel-based", num(median_dist_km, 0), ""),
                     Method = as.character(method),
                     !!sprintf("RMSE (%s)", GCM2) := num(RMSE, 1), !!sprintf("MAE (%s)", GCM2) := num(MAE, 1),
                     !!R2 := num(R2, 3), MEC = num(MEC, 3), `Bias (%)` = num(bias_pct, 1))
  notes <- c(sprintf("* n < %d: shown for completeness, not used in the summaries.", CFG$min_n_eco),
             "R\u00b2 is undefined (\u2013) where a method predicts one value for the whole ecoregion (Tier 1).")
  if (!is.null(mew)) notes <- c(notes, paste0(
    "Across ecoregions (n >= ", CFG$min_n_eco, "): median MEC ",
    paste(sprintf("%s %s", mew$method, num(mew$MEC_median, 3)), collapse = ", "),
    "; median |bias| ", paste(sprintf("%s %s %%", mew$method, num(mew$abs_bias_pct_median, 1)), collapse = ", "), "."))
  save_table("S10", "P4_S10_validation_ecoregion.csv",
             "Table S10. Accuracy for each ecoregion when it is left out of training (ecoregion design, pooled over seeds).",
             S10, notes)
}

# =============================================================================
# S11: distance to AGC samples and between-year CV by ecoregion
# =============================================================================
step("S11")
dce <- src("P4_distance_coverage_eco.csv"); cvs <- src("P4_cv_summary.csv")
if (!is.null(dce)) {
  cvw <- if (!is.null(cvs)) cvs |>
    dplyr::mutate(grp = dplyr::case_when(cv_label == "<10" ~ "lt10", cv_label == "10-20" ~ "c1020", TRUE ~ "ge20")) |>
    dplyr::group_by(ecoregion, grp) |> dplyr::summarise(p = sum(area_pct), .groups = "drop") |>
    tidyr::pivot_wider(names_from = grp, values_from = p, values_fill = 0) else
      tibble::tibble(ecoregion = character(0), lt10 = numeric(0), c1020 = numeric(0), ge20 = numeric(0))
  S11 <- dce |> dplyr::mutate(ecoregion = ifelse(ecoregion == "All", "Indonesia", ecoregion)) |>
    dplyr::left_join(cvw |> dplyr::mutate(ecoregion = ifelse(ecoregion == "All", "Indonesia", ecoregion)), by = "ecoregion") |>
    dplyr::arrange(ecoregion == "Indonesia", eco_id == 0, dplyr::desc(area_km2)) |>
    dplyr::transmute(Ecoregion = ecoregion, !!sprintf("Extent (%s)", KM2) := num(area_km2, 1),
                     `Within 1 km of a sample (%)` = num(within_1km_pct, 1),
                     `Beyond 5 km (%)` = num(beyond_5km_pct, 1), `Beyond 50 km (%)` = num(beyond_50km_pct, 1),
                     `Median distance class` = median_class,
                     `CV < 10 % (% of extent)` = num(lt10, 1), `CV 10-20 %` = num(c1020, 1), `CV >= 20 %` = num(ge20, 1))
  save_table("S11", "P4_S11_distance_cv_ecoregion.csv",
             "Table S11. Distance from the extent to the nearest AGC sample, and the between-year coefficient of variation (CV) of mapped AGC, by ecoregion.",
             S11, c("Distance: from each extent pixel to the nearest AGC sample of any year (P4_G03, 11 classes).",
                   "CV: SD / mean of the 8 annual AGC maps (2017-2024) per pixel; shares of the extent."))
}

# =============================================================================
# S6: AGC samples inside the extent
# =============================================================================
step("S6")
sie <- src("P4_samples_in_extent.csv")
if (!is.null(sie)) {
  S6 <- sie |> dplyr::mutate(level = factor(level, levels = c("all", "ecoregion", "year")),
                             key = ifelse(level == "ecoregion", -n, 0)) |>
    dplyr::arrange(level, key, group) |>
    dplyr::transmute(Group = dplyr::case_when(level == "all" ~ "All AGC samples",
                                              level == "year" ~ paste("Survey year", group), TRUE ~ group),
                     `AGC samples (n)` = num(n, 0), `In final extent (%)` = num(in_extent_pct, 1),
                     `MaxProb < 0.7 (%)` = num(fail_prob_pct, 1), `Depth outside 0-5 m (%)` = num(fail_depth_pct, 1),
                     !!sprintf("%s (%%)", TAB$trf_label) := num(class16_pct, 1),
                     !!sprintf("Mean AGC inside (%s)", GCM2) := num(mean_agc_in, 1),
                     !!sprintf("Mean AGC outside (%s)", GCM2) := num(mean_agc_out, 1))
  save_table("S6", "P4_S6_samples_in_extent.csv",
             "Table S6. AGC samples that fall inside the final extent, and the masks that exclude the others.", S6,
             c("A sample can fail more than one mask, so the mask columns do not add up to the share outside.",
               "Final extent = core + add-on: shallow (ACA depth 0-5 m), MaxProb >= 0.7 in the sample year, not Terrestrial Reef Flat (ACA geomorphic map)."))
}

# =============================================================================
# S2: Tier 1 by survey year and the survey composition
# =============================================================================
step("S2")
yco <- src("P4_year_composition.csv")
if (!is.null(t1g)) {
  comp <- if (!is.null(yco)) yco |> dplyr::group_by(year) |> dplyr::arrange(dplyr::desc(share), .by_group = TRUE) |>
    dplyr::summarise(main = paste(sprintf("%s (%.0f %%)", ecoregion[seq_len(min(2, dplyr::n()))],
                                          100 * share[seq_len(min(2, dplyr::n()))]), collapse = "; "),
                     .groups = "drop") |> dplyr::mutate(year = as.character(year)) else
      tibble::tibble(year = character(0), main = character(0))
  S2 <- t1g |> dplyr::mutate(year = as.character(year)) |> dplyr::left_join(comp, by = "year") |>
    dplyr::transmute(`Survey year` = ifelse(year == "pooled", "Pooled", year), `AGC samples (n)` = num(n, 0),
                     `Effective n` = num(n_eff, 1),
                     !!sprintf("Mean \u00b1 SD (%s)", GCM2) := pm(mean_gC_m2, sd_gC_m2, 2),
                     !!sprintf("95 %% CI half-width (%s)", GCM2) := num(ci_half_gC_m2, 2),
                     `Main ecoregions (share of the year's samples)` = ifelse(is.na(main), "all", main))
  save_table("S2", "P4_S2_tier1_by_year.csv",
             "Table S2. Field AGC by survey year (the Tier 1 factors) and the ecoregions that dominate each year's samples.",
             S2, c("Survey year and region are confounded: differences between years mostly reflect where the surveys took place.",
                   "No CI where the effective sample size is below 2."))
}

# =============================================================================
# S4: the 22 ROIs
# =============================================================================
step("S4")
rc <- src("P4_roi_comparison.csv"); sr <- src("P4_shrinkage_roi.csv")
# Site names and ecoregion from the data: AGC samples inside each ROI polygon (st_within,
# as P4_07) joined to the location names of the field label file by gee_id
roi_poly <- if (!is.na(CFG$roi_path)) CFG$roi_path else
  file.path(if (!is.null(CFG$gee_dir)) CFG$gee_dir else out_file("gee"), "P4_roi_polygons.geojson")
roi_info <- tibble::tibble(loc_id = TAB$roi22, site = NA_character_, lon = NA_real_, lat = NA_real_,
                           ecoregion = NA_character_)
if (file.exists(roi_poly)) {
  sf::sf_use_s2(FALSE)
  poly <- sf::st_transform(sf::st_read(roi_poly, quiet = TRUE), 4326)
  poly$loc_id <- as.character(poly[[CFG$roi_id_field]])
  poly <- poly[poly$loc_id %in% TAB$roi22, "loc_id"]
  cen <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(poly))))
  pts <- sf::st_as_sf(agc, coords = c("xcoord", "ycoord"), crs = 4326, remove = FALSE)
  inr <- suppressMessages(sf::st_join(pts, poly, join = sf::st_within)) |> sf::st_drop_geometry() |>
    dplyr::filter(!is.na(loc_id))
  locs <- NULL; fld <- character(0)
  if (file.exists(TAB$location_path)) {
    hdr <- names(readr::read_csv(TAB$location_path, n_max = 0, show_col_types = FALSE))
    fld <- if (TAB$location_field %in% hdr) TAB$location_field else
      utils::head(setdiff(grep("^loc|site|lokasi", hdr, ignore.case = TRUE, value = TRUE), c("loc_id", "gee_id")), 1)
    if (!identical(fld, TAB$location_field))
      cat(sprintf("Columns in %s:\n  %s\n", basename(TAB$location_path), paste(hdr, collapse = ", ")))
    if (length(fld) == 1 && "gee_id" %in% hdr)
      locs <- readr::read_csv(TAB$location_path, col_select = dplyr::all_of(c("gee_id", fld)),
                              col_types = readr::cols(.default = "c")) |>
        dplyr::rename(location = dplyr::all_of(fld)) |> dplyr::distinct(gee_id, .keep_all = TRUE)
  } else cat(sprintf("  -> location file not found: %s\n", TAB$location_path))
  check(sprintf("Location names read from %s (column '%s')", basename(TAB$location_path),
                if (length(fld) == 1) fld else "none"),
        !is.null(locs) && nrow(locs) > 0, sprintf("%d rows", if (is.null(locs)) 0L else nrow(locs)))
  if (!is.null(locs)) inr <- inr |> dplyr::left_join(locs, by = "gee_id") else inr$location <- NA_character_
  inr$location <- trimws(inr$location)
  sites <- inr |> dplyr::filter(!is.na(location), location != "") |>
    dplyr::count(loc_id, location, name = "n") |> dplyr::arrange(loc_id, dplyr::desc(n))
  # Display name: TAB$site_names, else title case of the code
  unmapped <- setdiff(unique(sites$location), names(TAB$site_names))
  sites$site_name <- ifelse(sites$location %in% names(TAB$site_names),
                            unname(TAB$site_names[sites$location]),
                            tools::toTitleCase(gsub("_", " ", sites$location)))
  check("Every location code of the 22 ROIs has a display name in TAB$site_names", length(unmapped) == 0,
        if (length(unmapped) > 0) sprintf("not listed (title case used): %s", paste(unmapped, collapse = ", ")) else
          sprintf("%d codes", length(unique(sites$location))))
  readr::write_excel_csv(sites, file.path(TAB$dir, "P4_S4_roi_sites.csv"))
  top <- sites |> dplyr::group_by(loc_id) |>
    dplyr::summarise(site = paste0(dplyr::first(site_name),
                                   ifelse(dplyr::n() > 1, sprintf(" (+%d)", dplyr::n() - 1), "")), .groups = "drop")
  eco_mode <- inr |> dplyr::count(loc_id, ecoregion) |> dplyr::group_by(loc_id) |>
    dplyr::slice_max(n, n = 1, with_ties = FALSE) |> dplyr::ungroup() |> dplyr::select(loc_id, ecoregion)
  roi_info <- tibble::tibble(loc_id = poly$loc_id, lon = cen[, 1], lat = cen[, 2]) |>
    dplyr::left_join(top, by = "loc_id") |> dplyr::left_join(eco_mode, by = "loc_id")
  check("ROI polygons found for all 22 ROIs", all(TAB$roi22 %in% roi_info$loc_id),
        sprintf("%d of %d (%s)", sum(TAB$roi22 %in% roi_info$loc_id), length(TAB$roi22), basename(roi_poly)))
  check("Every ROI has a location name", all(!is.na(roi_info$site)),
        if (any(is.na(roi_info$site))) sprintf("missing: %s", paste(roi_info$loc_id[is.na(roi_info$site)], collapse = ", ")) else "")
  if (!is.null(rc)) {
    nchk <- inr |> dplyr::count(loc_id, year, name = "n_join") |>
      dplyr::inner_join(rc |> dplyr::select(loc_id, dom_year, n_year), by = c("loc_id", "year" = "dom_year"))
    check("Samples per ROI in the dominant year match P4_roi_comparison", all(nchk$n_join == nchk$n_year) &&
            nrow(nchk) == length(TAB$roi22), sprintf("%d ROIs compared", nrow(nchk)))
  }
} else check("ROI polygons found", FALSE, sprintf("%s not found: S4 without site names", roi_poly))
if (!is.null(rc) && !is.null(sr)) {
  r <- rc |> dplyr::filter(loc_id %in% TAB$roi22) |>
    dplyr::left_join(sr |> dplyr::select(loc_id, year, n_oof = n, mean_pred_oof, pct_diff), by = c("loc_id", "dom_year" = "year")) |>
    dplyr::left_join(roi_info, by = "loc_id") |> dplyr::arrange(field_mean_gC_m2)
  check("All 22 ROIs found in P4_roi_comparison", nrow(r) == length(TAB$roi22), sprintf("%d of %d", nrow(r), length(TAB$roi22)))
  check("Out-of-fold ROI rows match the dominant year and n", all(!is.na(r$mean_pred_oof)) && all(r$n_oof == r$n_year, na.rm = TRUE))
  S4 <- r |> dplyr::transmute(
    ROI = paste0(loc_id, ifelse(n_year < TAB$roi_min_n, " *", "")), Site = dplyr::coalesce(site, "\u2013"),
    Ecoregion = dplyr::coalesce(ecoregion, "\u2013"), Longitude = num(lon, 3), Latitude = num(lat, 3),
    `Survey year` = as.character(dom_year), `AGC samples (n)` = num(n_year, 0),
    !!sprintf("Extent in ROI (%s)", KM2) := num(area_km2, 3),
    !!sprintf("Field mean \u00b1 SD (%s)", GCM2) := pm(field_mean_gC_m2, field_sd_gC_m2, 1),
    !!sprintf("Final map mean (%s)", GCM2) := num(px_mean_gC_m2, 1),
    !!sprintf("Out-of-fold mean (%s)", GCM2) := num(mean_pred_oof, 1),
    `Tier 1, local (t C)` = num(t1_local_tC, 2), `Pixel-based (t C)` = num(px_total_tC, 2),
    `Difference, final map (%)` = num(bias_px_pct, 1), `Difference, out-of-fold (%)` = num(pct_diff, 1))
  fit <- function(f, m) { k <- !is.na(f) & !is.na(m)
    c(slope = unname(coef(lm(m[k] ~ f[k]))[2]), r2 = cor(f[k], m[k])^2,
      mec = 1 - sum((f[k] - m[k])^2) / sum((f[k] - mean(f[k]))^2), rmse = sqrt(mean((f[k] - m[k])^2))) }
  big <- r$n_year >= TAB$roi_min_n
  fi <- fit(r$field_mean_gC_m2[big], r$px_mean_gC_m2[big]); fo <- fit(r$field_mean_gC_m2[big], r$mean_pred_oof[big])
  derive("roi_n_ge_min", sum(big), sprintf("ROIs of the 22 with n >= %d in the dominant year", TAB$roi_min_n))
  for (k in names(fi)) derive(paste0("roi_final_", k), fi[[k]], "final map mean regressed on field mean, ROIs n >= min")
  for (k in names(fo)) derive(paste0("roi_oof_", k), fo[[k]], "out-of-fold mean regressed on field mean, ROIs n >= min")
  save_table("S4", "P4_S4_roi22.csv",
             "Table S4. The 22 ROIs: field mean AGC in the survey year with the most samples, the final map (ROI samples in training) and the out-of-fold prediction (meadow cluster held out).",
             S4, c(sprintf("* n < %d.", TAB$roi_min_n),
                    "Extent in ROI: final extent inside the ROI polygon. Tier 1, local: field mean x extent in ROI. Difference: (map - field) / field x 100.",
                    sprintf("ROIs with n >= %d (%d): final map slope %.2f, R\u00b2 %.2f, MEC %.2f, RMSE %.1f; out-of-fold slope %.2f, R\u00b2 %.2f, MEC %.2f, RMSE %.1f (map mean regressed on field mean).",
                            TAB$roi_min_n, sum(big), fi[["slope"]], fi[["r2"]], fi[["mec"]], fi[["rmse"]],
                            fo[["slope"]], fo[["r2"]], fo[["mec"]], fo[["rmse"]]),
                    sprintf("Site: location of most AGC samples in the ROI (%s, display names from TAB$site_names); (+k) = k further names (P4_S4_roi_sites.csv). Ecoregion: of most samples in the ROI. Longitude, latitude: ROI polygon centroid.",
                            basename(TAB$location_path))))
}

# =============================================================================
# S5: presence model cross-validation (new, 9 Nov 2026; GEE script P4_G09)
# =============================================================================
step("S5")
pa_cv_f <- file.path(TAB$gee_dir, "P4_gee_pa_cv.csv"); pa_n_f <- file.path(TAB$gee_dir, "P4_gee_pa_cv_counts.csv")
if (file.exists(pa_cv_f)) {
  pcv <- readr::read_csv(pa_cv_f, show_col_types = FALSE)
  file.copy(pa_cv_f, file.path(TAB$dir, "inputs"), overwrite = TRUE)
  check("Presence CV: TP + FP + FN + TN = n at every threshold", all(pcv$TP + pcv$FP + pcv$FN + pcv$TN == pcv$n),
        sprintf("n = %s", paste(unique(pcv$n), collapse = ", ")))
  S5 <- pcv |> dplyr::arrange(threshold) |> dplyr::transmute(
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
  save_table("S5", "P4_S5_presence_cv.csv",
             "Table S5. Five-fold cross-validation of the seagrass presence model at three probability thresholds.",
             S5, c("Random forest in GEE (smileRandomForest: 500 trees, 10 variables per split, minimum leaf 3, bag fraction 0.8, seed 42) on the 64 GSE bands and ACA depth; training points 2018-2023; random folds (seed 42). The folds and settings are those of the presence model that defined the extent.",
                   "Presences: field quadrats with cover > 0; absences: quadrats with zero cover and non-seagrass areas of the Allen Coral Atlas.",
                   "The extent uses 0.7 (MaxProb over the eight annual probability layers).",
                   cnt_note))
} else cat(sprintf("  -> missing: %s (copy P4_gee_pa_cv.csv and P4_gee_pa_cv_counts.csv from Drive P4v2_GEE)\n", pa_cv_f))

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
tab_order <- function(id) ifelse(substr(id, 1, 1) == "T", 0, 100) + as.integer(substring(id, 2))
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
inp <- setdiff(used, "P4_samples.csv")         # P4_samples.csv is large and is not copied
ok_copy <- file.copy(out_file(inp), file.path(TAB$dir, "inputs"), overwrite = TRUE)
cat(sprintf("\n%d tables written to %s\n", length(TABLES), TAB$dir))
cat(sprintf("Copied %d source CSVs to %s\n", sum(ok_copy), file.path(TAB$dir, "inputs")))
n_cek <- sum(grepl("^\\[CEK", CHECKS))
cat(sprintf("Consistency checks: %d OK, %d CEK\n", sum(grepl("^\\[OK", CHECKS)), n_cek))
if (n_cek > 0) cat("Some checks failed: see the [CEK] lines in P4_09_checks.txt before using the tables.\n")
cat("\nP4_09_tables.R COMPLETE\n")
