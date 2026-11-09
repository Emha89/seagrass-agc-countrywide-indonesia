# =============================================================================
# P4_07_gee_results.R
# Paper 4 revised design: joins the GEE exports (P4_G01-G04) with the R
# validation results into the estimation tables (national, ecoregion, ROI) and
# the coverage and temporal checks.
#
# Copy the CSV/GeoJSON files from the Drive folder P4v2_GEE into
# <out_dir>/gee (or set CFG$gee_dir). A missing GEE file only skips its step,
# so the script can be run as the exports arrive.
#
# Reads (GEE): P4_gee_totals_YYYY.csv, P4_gee_points_fitted.csv,
#   P4_gee_roi_totals.csv, P4_roi_polygons.geojson, P4_gee_samples_extent.csv,
#   P4_gee_pa_points.csv, P4_gee_distance_area.csv, P4_gee_distance_pct.csv,
#   P4_gee_cv_summary.csv, P4_gee_gse_similarity_eco.csv,
#   P4_gee_gse_sim_sample.csv, P4_gee_gse_sim_points.csv, P4_gee_importance.csv
# Reads (R): P4_samples.csv, gee_upload/P4_gee_training.csv, P4_meow_lookup.csv,
#   P4_uncertainty_params.csv, P4_metrics_summary.csv, P4_tier2_class_means.csv
# Writes:
#   STEP 1  P4_pixel_national.csv              pixel-based totals 2017-2024, CI from
#                                              the cluster-robust SE of the mean residual
#           P4_pixel_ci_sensitivity.csv        same totals with the 09b-style CI
#           gee_cover_class_area_by_year.csv   class areas (input of P4_05 STEP 7)
#   STEP 2  P4_tier1_national_gee.csv          Tier 1 with the GEE extent area
#   STEP 3  P4_tier2_class_means_map.csv       class means, classes from the GEE cover map
#           P4_tier2_loyo_map.csv, P4_tier2_national_map.csv
#   STEP 4  P4_national_three_methods.csv
#   STEP 5  P4_ecoregion_by_year.csv, P4_ecoregion_summary.csv
#   STEP 6  P4_roi_comparison.csv
#   STEP 7  P4_samples_in_extent.csv
#   STEP 8  P4_pa_threshold_by_year.csv
#   STEP 9  P4_distance_coverage_class.csv, P4_distance_coverage_design.csv,
#           P4_distance_coverage_eco.csv, P4_distance_percentiles.csv
#   STEP 10 P4_cv_summary.csv
#   STEP 11 P4_gse_similarity_eco.csv (optional GEE task), P4_gse_similarity_eco_sample.csv,
#           P4_gse_similarity_pairs.csv,
#           P4_gse_similarity_points_lag.csv
# Add-on (2 Nov 2026): when the P4_G07 files (*_addon*.csv) are in gee_dir, the
# core sums of P4_G02/G03/G04 and the add-on sums are added (see "Core + add-on").
# Units: gC m-2 x km2 = t C. All plots are printed; ggsave lines are commented.
# =============================================================================

# setwd("path/to/seagrass-agc-countrywide-indonesia/R")   # folder that holds these scripts
source("P4_00_config_functions.R")

gee_dir   <- if (!is.null(CFG$gee_dir)) CFG$gee_dir else out_file("gee")
NA_CODE   <- -9999
YEARS_MAP <- 2017:2024
DIST_LABELS <- c("<=10 m", "10-50 m", "50-100 m", "100-250 m", "250-500 m", "500 m-1 km",
                 "1-2 km", "2-5 km", "5-20 km", "20-50 km", ">50 km")
# validation design whose median distance (P4_04) lies in each distance class
DIST_DESIGN <- c("pixel", "b100", "b100", "b1000", "b1000", "b1000",
                 "cluster", "cluster", "cluster", "ecoregion", "ecoregion")
CV_LABELS   <- c("<10", "10-20", "20-30", "30-50", "50-100", ">=100")

gee_path <- function(name) file.path(gee_dir, name)
has_gee  <- function(name) {
  ok <- file.exists(gee_path(name))
  if (!ok) cat(sprintf("  -> skipped: %s not found in %s\n", name, gee_dir))
  ok
}
read_gee <- function(name) readr::read_csv(gee_path(name), show_col_types = FALSE, guess_max = 100000)
step <- function(txt) cat(sprintf("\n==== %s ====\n", txt))

# GEE exports are joined on rid. This check catches an outdated training asset
# in GEE (rid numbering from an older P4_06 run): gee_id must match at every rid.
check_rid <- function(d, name) {
  if (!all(c("rid", "gee_id") %in% names(d))) return(invisible(TRUE))
  ref <- samples$gee_id[match(d$rid, samples$rid)]
  bad <- sum(is.na(ref) | ref != as.character(d$gee_id))
  if (bad > 0) stop(sprintf("%s: %d rows whose rid does not match gee_id. The GEE training asset is older than P4_samples.csv: re-upload P4_gee_training and rerun the GEE scripts.", name, bad))
  invisible(TRUE)
}

wquantile <- function(x, w, p) {
  ok <- is.finite(x) & is.finite(w); x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; cw <- cumsum(w[o]) / sum(w)
  vapply(p, function(pp) x[which(cw >= pp)[1]], numeric(1))
}

# ---- R-side inputs -----------------------------------------------------------
samples <- read_samples() |> dplyr::mutate(rid = dplyr::row_number())
up <- readr::read_csv(out_file("gee_upload/P4_gee_training.csv"), show_col_types = FALSE,
                      col_select = c("rid", "gee_id", "year"), col_types = readr::cols(gee_id = "c"))
chk <- dplyr::inner_join(samples[, c("rid", "gee_id", "year")], up, by = "rid")
if (nrow(chk) != nrow(samples) || any(chk$gee_id.x != chk$gee_id.y))
  stop("P4_samples.csv and gee_upload/P4_gee_training.csv do not match row by row: rerun P4_06.")
agc <- samples |> dplyr::filter(is_agc)

lookup   <- readr::read_csv(out_file("P4_meow_lookup.csv"), show_col_types = FALSE)
eco_name <- function(id) ifelse(id == 0, "Outside MEOW polygons",
                                lookup$ECOREGION[match(id, lookup$eco_id)])

unc  <- readr::read_csv(out_file("P4_uncertainty_params.csv"), show_col_types = FALSE)
upar <- setNames(unc$value, unc$parameter)
DEFF <- upar[["DEFF"]]
cat(sprintf("DEFF %.2f | SE of mean residual %.2f gC m-2 | SD_total %.2f gC m-2\n",
            DEFF, upar[["se_mean_resid_gC_m2"]], upar[["SD_total_gC_m2"]]))

# =============================================================================
# Core + add-on (2 Nov 2026). P4_G01 maps the core extent; P4_G06/G07 map and sum
# the add-on area (shallow pixels without ACA geomorphic data, dropped by the 09b
# mask by mistake). When the add-on files are present, every GEE sum below is
# core + add-on, i.e. the intended extent: shallow, MaxProb >= 0.7, only class 16
# removed. Without them, the core extent is used and a NOTE is printed.
# =============================================================================
addon_files <- list.files(gee_dir, pattern = "_addon(_\\d{4})?\\.csv$")
ADDON <- length(addon_files) > 0
cat(if (ADDON) "\nEXTENT: core (P4_G01) + add-on (P4_G06/G07) = only ACA class 16 removed\n" else
      "\nEXTENT: core only (P4_G01). NOTE: shallow pixels without ACA geomorphic data are missing; run P4_G06/G07.\n")

# Read a core table and, when present, its add-on table, then add the sums.
# keys = columns that identify a row; all other numeric columns are summed.
# optional = TRUE: a missing add-on table gives a NOTE and the core table only
# (used for the CV summary, a descriptive check), instead of stopping.
read_core_addon <- function(core_name, addon_name, keys, optional = FALSE) {
  d <- read_gee(core_name)
  if (ADDON) {
    if (!file.exists(gee_path(addon_name))) {
      if (!optional)
        stop(sprintf("Add-on files found but %s is missing: copy every P4_G07 output.", addon_name))
      cat(sprintf("  NOTE: %s not found: this table covers the core extent only.\n", addon_name))
    } else {
      d <- dplyr::bind_rows(d, read_gee(addon_name))
    }
  }
  num <- setdiff(names(d)[vapply(d, is.numeric, logical(1))], keys)
  d |> dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(num), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
}

# =============================================================================
# STEP 1: Pixel-based national totals and class areas
# =============================================================================
step("STEP 1: pixel-based national totals")
tot_files <- list.files(gee_dir, pattern = "^P4_gee_totals_\\d{4}\\.csv$", full.names = TRUE)
tot <- NULL
if (length(tot_files) == 0) {
  cat("  -> skipped: no P4_gee_totals_YYYY.csv (STEPS 1-5 need them)\n")
} else {
  tot_raw <- purrr::map_dfr(tot_files, readr::read_csv, show_col_types = FALSE)
  # core extent area (P4_G01/G02 only); STEP 12 scales the P4_G05 mask ratios to it
  AREA_CORE_KM2 <- mean(tapply(tot_raw$area_m2, tot_raw$year, sum)) / 1e6
  # core area per ecoregion: weights of the GSE similarity sample (STEP 11),
  # which P4_G04 drew inside the core extent only
  ECO_AREA_CORE <- tot_raw |> dplyr::filter(year == min(year)) |>
    dplyr::mutate(eco_id = as.integer(round(eco_id))) |>
    dplyr::group_by(eco_id) |> dplyr::summarise(area = sum(area_m2), .groups = "drop")
  if (ADDON) {
    add_files <- list.files(gee_dir, pattern = "^P4_gee_totals_addon_\\d{4}\\.csv$", full.names = TRUE)
    y_core <- sort(unique(tot_raw$year))
    y_add  <- as.integer(sub(".*_(\\d{4})\\.csv$", "\\1", basename(add_files)))
    if (!setequal(y_core, y_add))
      stop("Add-on totals do not cover the same years as the core totals: copy all P4_gee_totals_addon_YYYY.")
    tot_raw <- dplyr::bind_rows(tot_raw, purrr::map_dfr(add_files, readr::read_csv, show_col_types = FALSE))
  }
  tot <- tot_raw |>
    dplyr::mutate(eco_id = as.integer(round(eco_id)), cover_class = as.integer(round(cover_class))) |>
    dplyr::group_by(year, eco_id, cover_class) |>
    dplyr::summarise(area_m2 = sum(area_m2), agc_g = sum(agc_g), n_px = sum(n_px), .groups = "drop") |>
    dplyr::mutate(year = as.integer(year), eco_id = as.integer(round(eco_id)),
                  cover_class = as.integer(round(cover_class)),
                  class = factor(CFG$cover_labels[cover_class], levels = CFG$cover_labels),
                  ecoregion = eco_name(eco_id))
  miss_y <- setdiff(YEARS_MAP, unique(tot$year))
  if (length(miss_y) > 0) cat("  NOTE: years without a totals file:", paste(miss_y, collapse = ", "), "\n")

  nat_px <- tot |>
    dplyr::group_by(year) |>
    dplyr::summarise(area_km2 = sum(area_m2) / 1e6, n_px = sum(n_px),
                     total_tC = sum(agc_g) / 1e6, .groups = "drop") |>
    dplyr::mutate(mean_gC_m2 = total_tC / area_km2,
                  ci_half_tC = 1.96 * upar[["se_mean_resid_gC_m2"]] * area_km2,
                  ci_low_tC = total_tC - ci_half_tC, ci_up_tC = total_tC + ci_half_tC,
                  field_data = year %in% unique(agc$year),
                  ci_note = ifelse(field_data, "",
                                   "no field data this year; CI assumes the 2017-2023 bias SE"))
  print_tbl(nat_px |> dplyr::select(-ci_note), 1)
  readr::write_csv(nat_px, out_file("P4_pixel_national.csv"))

  rel_area <- diff(range(nat_px$area_km2)) / mean(nat_px$area_km2)
  cat(sprintf("Extent area: %.2f km2 (relative range across years %.1e; expected 0, fixed extent)\n",
              mean(nat_px$area_km2), rel_area))
  if (ADDON) cat(sprintf("  = core %.2f km2 (P4_G02) + add-on %.2f km2 (P4_G07)\n",
                         AREA_CORE_KM2, mean(nat_px$area_km2) - AREA_CORE_KM2))
  cat(sprintf("Extent area in CFG$extent_km2 (used by P4_05): %.2f km2\n", CFG$extent_km2))
  if (nrow(nat_px |> dplyr::filter(!field_data)) > 0)
    cat("Years without field data (2024) use the 2017-2023 bias SE: shown, labelled without validation.\n",
        "Support for this assumption: the 2023-2024 GSE similarity in STEP 11.\n")
  AREA_KM2 <- mean(nat_px$area_km2)

  ci_sens <- nat_px |>
    dplyr::transmute(year, total_tC,
                     ci_half_bias_tC = ci_half_tC,
                     ci_half_sdtotal_neff_tC = 1.96 * upar[["SD_total_gC_m2"]] * area_km2 /
                       sqrt(upar[["N_eff"]]))
  readr::write_csv(ci_sens, out_file("P4_pixel_ci_sensitivity.csv"))

  cls_area <- tot |>
    dplyr::group_by(year, class, .drop = FALSE) |>
    dplyr::summarise(area_km2 = sum(area_m2) / 1e6, .groups = "drop") |>
    dplyr::mutate(class = as.character(class))
  readr::write_csv(cls_area, out_file("gee_cover_class_area_by_year.csv"))
  cat("Cover class areas (km2):\n")
  print(as.data.frame(tidyr::pivot_wider(cls_area, names_from = class, values_from = area_km2)),
        row.names = FALSE)
}

# =============================================================================
# STEP 2: Tier 1 national totals with the GEE extent area
# =============================================================================
step("STEP 2: Tier 1 national (GEE extent area)")
tier1 <- NULL
if (!is.null(tot)) {
  tier1 <- dplyr::bind_rows(
    agc |> dplyr::group_by(year) |>
      dplyr::group_modify(~ tier1_summary(.x$AGC_pred, DEFF, AREA_KM2)) |>
      dplyr::ungroup() |> dplyr::mutate(year = as.character(year)),
    tier1_summary(agc$AGC_pred, DEFF, AREA_KM2) |> dplyr::mutate(year = "pooled")
  )
  print_tbl(tier1, 1)
  readr::write_csv(tier1, out_file("P4_tier1_national_gee.csv"))
}

# =============================================================================
# STEP 3: Tier 2 with classes from the GEE cover map
# (class at each AGC sample = map class at that pixel, as in stratified estimation)
# =============================================================================
step("STEP 3: Tier 2 (classes from the GEE cover map)")
class_means <- function(d) {
  d |>
    dplyr::group_by(cov_class, .drop = FALSE) |>
    dplyr::summarise(n = dplyr::n(), mean_gC_m2 = mean(AGC_pred), sd_gC_m2 = sd(AGC_pred),
                     .groups = "drop") |>
    dplyr::mutate(n_eff = n / DEFF)
}
tier2_total <- function(areas, cm, overall_mean) {
  d <- dplyr::left_join(areas, cm, by = c("class" = "cov_class"))
  empty <- is.na(d$mean_gC_m2) | d$n == 0
  d$mean_gC_m2[empty] <- overall_mean
  d$var_term <- ifelse(empty | d$n_eff < 2, NA_real_, d$area_km2^2 * d$sd_gC_m2^2 / d$n_eff)
  tibble::tibble(total_tC = sum(d$mean_gC_m2 * d$area_km2),
                 ci_half_tC = 1.96 * sqrt(sum(d$var_term, na.rm = TRUE)),
                 classes_without_samples = paste(d$class[empty], collapse = ", "))
}

t2 <- NULL; cm_map <- NULL; agc_map <- NULL
if (has_gee("P4_gee_points_fitted.csv")) {
  fit <- read_gee("P4_gee_points_fitted.csv")
  check_rid(fit, "P4_gee_points_fitted.csv")
  fit <- fit |> dplyr::select(rid, AGC_map, tSPC_map, cov_class_map)
  agc_map <- agc |>
    dplyr::left_join(fit, by = "rid") |>
    dplyr::mutate(cov_class = factor(CFG$cover_labels[cov_class_map], levels = CFG$cover_labels))
  n_miss <- sum(is.na(agc_map$cov_class_map))
  if (n_miss > 0) cat(sprintf("  NOTE: %d AGC samples without a fitted value from GEE\n", n_miss))
  agc_map <- agc_map |> dplyr::filter(!is.na(cov_class))

  ins <- agc_map |> dplyr::summarise(n = dplyr::n(), r = cor(AGC_pred, AGC_map),
                                     mean_obs = mean(AGC_pred), mean_fit = mean(AGC_map))
  cat("GEE fitted AGC at the training points (in-sample, sanity check only):\n")
  print_tbl(ins, 3)

  cm_map <- class_means(agc_map)
  cat("Class means, classes from the GEE cover map:\n")
  print_tbl(cm_map, 2)
  readr::write_csv(cm_map, out_file("P4_tier2_class_means_map.csv"))
  if (file.exists(out_file("P4_tier2_class_means.csv"))) {
    prov <- readr::read_csv(out_file("P4_tier2_class_means.csv"), show_col_types = FALSE) |>
      dplyr::select(cov_class, n_R = n, mean_R = mean_gC_m2)
    cat("Compared with the provisional R classes (P4_05):\n")
    print_tbl(dplyr::left_join(cm_map |> dplyr::mutate(cov_class = as.character(cov_class)) |>
                                 dplyr::select(cov_class, n_map = n, mean_map = mean_gC_m2),
                               prov |> dplyr::mutate(cov_class = as.character(cov_class)),
                               by = "cov_class"), 2)
  }

  loyo_map <- purrr::map_dfr(sort(unique(agc_map$year)), function(y) {
    class_means(agc_map |> dplyr::filter(year != y)) |> dplyr::mutate(left_out_year = y, .before = 1)
  })
  readr::write_csv(loyo_map, out_file("P4_tier2_loyo_map.csv"))

  if (!is.null(tot)) {
    cm_c <- cm_map |> dplyr::mutate(cov_class = as.character(cov_class))
    t2 <- cls_area |>
      dplyr::group_by(year) |>
      dplyr::group_modify(~ tier2_total(.x, cm_c, mean(agc_map$AGC_pred))) |>
      dplyr::ungroup()
    t2_loyo <- purrr::map_dfr(unique(loyo_map$left_out_year), function(y) {
      cm_y <- loyo_map |> dplyr::filter(left_out_year == y) |>
        dplyr::mutate(cov_class = as.character(cov_class))
      cls_area |> dplyr::group_by(year) |>
        dplyr::group_modify(~ tier2_total(.x, cm_y, mean(agc_map$AGC_pred[agc_map$year != y]))) |>
        dplyr::ungroup()
    }) |>
      dplyr::group_by(year) |>
      dplyr::summarise(loyo_min_tC = min(total_tC), loyo_max_tC = max(total_tC), .groups = "drop")
    t2 <- dplyr::left_join(t2, t2_loyo, by = "year") |>
      dplyr::mutate(ci_low_tC = total_tC - ci_half_tC, ci_up_tC = total_tC + ci_half_tC)
    print_tbl(t2, 0)
    readr::write_csv(t2, out_file("P4_tier2_national_map.csv"))
  }
}

# =============================================================================
# STEP 4: Three methods side by side (national)
# =============================================================================
step("STEP 4: national totals, three methods")
if (!is.null(tot)) {
  t1_pooled <- tier1 |> dplyr::filter(year == "pooled")
  nat3 <- dplyr::bind_rows(
    nat_px |> dplyr::transmute(year, method = "Pixel-based", total_tC, ci_low_tC, ci_up_tC),
    tier1 |> dplyr::filter(year != "pooled") |>
      dplyr::transmute(year = as.integer(year), method = "Tier 1 (survey-year mean)",
                       total_tC, ci_low_tC, ci_up_tC),
    tibble::tibble(year = sort(unique(nat_px$year)), method = "Tier 1 (pooled mean)",
                   total_tC = t1_pooled$total_tC, ci_low_tC = t1_pooled$ci_low_tC,
                   ci_up_tC = t1_pooled$ci_up_tC),
    if (!is.null(t2)) t2 |> dplyr::transmute(year, method = "Tier 2 (pooled class means)",
                                             total_tC, ci_low_tC, ci_up_tC)
  ) |>
    dplyr::mutate(total_kt = total_tC / 1e3, ci_low_kt = ci_low_tC / 1e3, ci_up_kt = ci_up_tC / 1e3)
  readr::write_csv(nat3, out_file("P4_national_three_methods.csv"))
  print(as.data.frame(nat3 |> dplyr::select(year, method, total_kt) |>
                        tidyr::pivot_wider(names_from = method, values_from = total_kt) |>
                        dplyr::arrange(year)), row.names = FALSE, digits = 4)

  p_nat <- ggplot(nat3, aes(x = year, y = total_kt, colour = method, fill = method)) +
    geom_ribbon(aes(ymin = ci_low_kt, ymax = ci_up_kt), alpha = 0.12, colour = NA,
                data = nat3 |> dplyr::filter(method != "Tier 1 (survey-year mean)")) +
    geom_line(data = nat3 |> dplyr::filter(method != "Tier 1 (survey-year mean)")) +
    geom_pointrange(aes(ymin = ci_low_kt, ymax = ci_up_kt),
                    data = nat3 |> dplyr::filter(method == "Tier 1 (survey-year mean)"),
                    position = position_nudge(x = 0.1), size = 0.3) +
    scale_x_continuous(breaks = YEARS_MAP) +
    labs(title = "National seagrass AGC by method", x = NULL, y = "AGC (kt C)",
         colour = NULL, fill = NULL) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  print(p_nat)
  # ggsave(out_file("P4_fig_national_three_methods.png"), p_nat, width = 8, height = 5, dpi = 300)
}

# =============================================================================
# STEP 5: Ecoregion level (middle reporting unit)
# =============================================================================
step("STEP 5: ecoregions")
if (!is.null(tot)) {
  eco_year <- tot |>
    dplyr::group_by(year, eco_id, ecoregion) |>
    dplyr::summarise(area_km2 = sum(area_m2) / 1e6, px_total_tC = sum(agc_g) / 1e6, .groups = "drop") |>
    dplyr::mutate(px_mean_gC_m2 = px_total_tC / area_km2)
  if (!is.null(cm_map)) {
    t2_eco <- tot |>
      dplyr::group_by(year, eco_id, class) |>
      dplyr::summarise(area_km2 = sum(area_m2) / 1e6, .groups = "drop") |>
      dplyr::left_join(cm_map |> dplyr::select(class = cov_class, mean_gC_m2), by = "class") |>
      dplyr::group_by(year, eco_id) |>
      dplyr::summarise(t2_total_tC = sum(area_km2 * mean_gC_m2, na.rm = TRUE), .groups = "drop")
    eco_year <- dplyr::left_join(eco_year, t2_eco, by = c("year", "eco_id"))
  }
  readr::write_csv(eco_year, out_file("P4_ecoregion_by_year.csv"))

  eco_field <- agc |>
    dplyr::group_by(ecoregion) |>
    dplyr::summarise(n_agc = dplyr::n(), n_years = dplyr::n_distinct(year),
                     field_mean_gC_m2 = mean(AGC_pred), field_sd_gC_m2 = sd(AGC_pred), .groups = "drop")
  nat_mean <- mean(agc$AGC_pred)

  eco_sum <- eco_year |>
    dplyr::group_by(eco_id, ecoregion) |>
    # min and max first: summarise() evaluates in order, so after px_total_tC = mean(...)
    # a later min(px_total_tC) would see the mean (bug found 7 Nov 2026: min = max = mean)
    dplyr::summarise(area_km2 = mean(area_km2),
                     px_total_min_tC = min(px_total_tC), px_total_max_tC = max(px_total_tC),
                     px_total_tC = mean(px_total_tC),
                     t2_total_tC = if ("t2_total_tC" %in% names(eco_year)) mean(t2_total_tC) else NA_real_,
                     .groups = "drop") |>
    dplyr::relocate(px_total_tC, .before = px_total_min_tC) |>
    dplyr::left_join(eco_field, by = "ecoregion") |>
    dplyr::mutate(n_agc = dplyr::coalesce(n_agc, 0L),
                  px_mean_gC_m2 = px_total_tC / area_km2,
                  reported = n_agc >= CFG$min_n_eco,
                  t1_eco_total_tC = ifelse(reported, field_mean_gC_m2 * area_km2, NA_real_),
                  t1_eco_ci_half_tC = ifelse(reported & n_agc / DEFF >= 2,
                                             1.96 * field_sd_gC_m2 / sqrt(n_agc / DEFF) * area_km2,
                                             NA_real_),
                  t1_national_factor_tC = nat_mean * area_km2,
                  share_area_pct = 100 * area_km2 / sum(area_km2),
                  share_samples_pct = 100 * n_agc / sum(n_agc)) |>
    dplyr::arrange(dplyr::desc(area_km2))
  cat("Tier 1 ecoregion CI uses the national DEFF (an approximation).\n")
  print(as.data.frame(eco_sum |> dplyr::select(ecoregion, area_km2, share_area_pct, n_agc,
                                               share_samples_pct, px_mean_gC_m2, field_mean_gC_m2,
                                               px_total_tC, t1_eco_total_tC, t2_total_tC)),
        row.names = FALSE, digits = 3)
  readr::write_csv(eco_sum, out_file("P4_ecoregion_summary.csv"))

  p_share <- eco_sum |>
    dplyr::select(ecoregion, `Extent area` = share_area_pct, `AGC samples` = share_samples_pct) |>
    tidyr::pivot_longer(-ecoregion, names_to = "what", values_to = "pct") |>
    ggplot(aes(x = reorder(ecoregion, pct), y = pct, fill = what)) +
    geom_col(position = "dodge") + coord_flip() +
    labs(title = "Where the extent is vs where the samples are", x = NULL, y = "Share (%)", fill = NULL) +
    theme_minimal(base_size = 11)
  print(p_share)
  # ggsave(out_file("P4_fig_ecoregion_share.png"), p_share, width = 8, height = 5, dpi = 300)
}

# =============================================================================
# STEP 6: ROIs -- local field mean x area vs pixel-based vs the national factor
# Year = the survey year with most samples in the ROI; area = national extent
# inside the ROI.
# =============================================================================
step("STEP 6: ROI comparison")
roi_poly <- if (!is.na(CFG$roi_path)) CFG$roi_path else gee_path("P4_roi_polygons.geojson")
if (has_gee("P4_gee_roi_totals.csv") && file.exists(roi_poly)) {
  idf <- CFG$roi_id_field
  roi_tot <- read_core_addon("P4_gee_roi_totals.csv", "P4_gee_roi_totals_addon.csv", c(idf, "year")) |>
    dplyr::rename(loc_id = dplyr::all_of(idf)) |>
    dplyr::mutate(loc_id = as.character(loc_id), year = as.integer(year))

  sf::sf_use_s2(FALSE)
  roi <- sf::st_transform(sf::st_read(roi_poly, quiet = TRUE), 4326)
  pts <- sf::st_as_sf(agc, coords = c("xcoord", "ycoord"), crs = 4326, remove = FALSE)
  in_roi <- suppressMessages(sf::st_join(pts, roi[, idf], join = sf::st_within)) |>
    sf::st_drop_geometry() |>
    dplyr::rename(loc_id = dplyr::all_of(idf)) |>
    dplyr::filter(!is.na(loc_id)) |>
    dplyr::mutate(loc_id = as.character(loc_id))

  roi_cmp <- in_roi |>
    dplyr::group_by(loc_id) |>
    dplyr::mutate(n_all = dplyr::n(), years = paste(sort(unique(year)), collapse = ";")) |>
    dplyr::count(loc_id, year, n_all, years, name = "n_year") |>
    dplyr::slice_max(n_year, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::rename(dom_year = year) |>
    dplyr::left_join(in_roi |> dplyr::group_by(loc_id, year) |>
                       dplyr::summarise(field_mean_gC_m2 = mean(AGC_pred), field_sd_gC_m2 = sd(AGC_pred),
                                        .groups = "drop"),
                     by = c("loc_id", "dom_year" = "year")) |>
    dplyr::left_join(roi_tot, by = c("loc_id", "dom_year" = "year")) |>
    dplyr::mutate(area_km2 = area_m2 / 1e6,
                  px_total_tC = agc_g / 1e6,
                  px_mean_gC_m2 = px_total_tC / area_km2,
                  t1_local_tC = field_mean_gC_m2 * area_km2,
                  t1_national_tC = mean(agc$AGC_pred) * area_km2,
                  bias_px_pct = 100 * (px_mean_gC_m2 - field_mean_gC_m2) / field_mean_gC_m2,
                  bias_national_pct = 100 * (mean(agc$AGC_pred) - field_mean_gC_m2) / field_mean_gC_m2,
                  px_closer = abs(bias_px_pct) < abs(bias_national_pct),
                  small_n = n_year < 10) |>
    dplyr::filter(!is.na(area_km2), area_km2 > 0) |>
    dplyr::arrange(field_mean_gC_m2)
  print(as.data.frame(roi_cmp |> dplyr::select(loc_id, dom_year, n_year, area_km2, field_mean_gC_m2,
                                               px_mean_gC_m2, bias_px_pct, bias_national_pct)),
        row.names = FALSE, digits = 3)
  cat(sprintf("ROIs: %d | pixel-based closer than the national factor in %d | median |bias| pixel %.0f%%, national factor %.0f%%\n",
              nrow(roi_cmp), sum(roi_cmp$px_closer), median(abs(roi_cmp$bias_px_pct)),
              median(abs(roi_cmp$bias_national_pct))))
  cat("Pixel-based ROI values here are fitted where the ROI samples trained the model;\n",
      "for the out-of-fold check set CFG$roi_path to this GeoJSON and rerun P4_04 (shrinkage, ROI).\n")
  readr::write_csv(roi_cmp, out_file("P4_roi_comparison.csv"))

  p_roi <- ggplot(roi_cmp, aes(x = field_mean_gC_m2, y = px_mean_gC_m2)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_hline(yintercept = mean(agc$AGC_pred), colour = "firebrick", linetype = 3) +
    geom_point(aes(size = n_year), alpha = 0.7) +
    labs(title = "ROI mean AGC: field vs pixel-based (dominant survey year)",
         subtitle = "dashed = 1:1, dotted = national Tier 1 factor",
         x = "Field mean (gC m-2)", y = "Pixel-based mean over the ROI extent (gC m-2)", size = "n") +
    theme_minimal(base_size = 11)
  print(p_roi)
  # ggsave(out_file("P4_fig_roi_field_vs_pixel.png"), p_roi, width = 6, height = 5, dpi = 300)
} else if (!file.exists(roi_poly)) {
  cat(sprintf("  -> skipped: ROI polygons not found (%s)\n", roi_poly))
}

# =============================================================================
# STEP 7: Field samples inside the national extent
# =============================================================================
step("STEP 7: samples inside the extent")
if (has_gee("P4_gee_samples_extent.csv")) {
  se <- read_gee("P4_gee_samples_extent.csv")
  check_rid(se, "P4_gee_samples_extent.csv")
  se <- se |> dplyr::select(rid, maxprob, depth_raw, geomorphic, in_extent)
  s7 <- samples |>
    dplyr::left_join(se, by = "rid") |>
    dplyr::mutate(in_extent = dplyr::coalesce(in_extent == 1, FALSE),
                  prob_ok  = dplyr::coalesce(maxprob >= 0.7, FALSE),
                  depth_ok = dplyr::coalesce(depth_raw > 0 & depth_raw <= 500, FALSE),
                  geo_ok   = dplyr::coalesce(!(geomorphic %in% c(16, NA_CODE)), FALSE),
                  # mask variants (see P4_G05): only class 16 removed / no geomorphic mask
                  in_only16 = prob_ok & depth_ok & dplyr::coalesce(geomorphic != 16, FALSE),
                  in_nogeo  = prob_ok & depth_ok,
                  geo_class16 = dplyr::coalesce(geomorphic == 16, FALSE),
                  geo_no_aca  = dplyr::coalesce(geomorphic == NA_CODE, TRUE))
  mism <- sum(s7$in_extent != (s7$prob_ok & s7$depth_ok & s7$geo_ok))
  cat(sprintf("Check: core in_extent differs from the three masks in %d samples (expect ~0)\n", mism))
  s7$in_extent_core <- s7$in_extent
  if (ADDON) s7$in_extent <- s7$in_only16      # core + add-on = only class 16 removed

  sum_ext <- function(d) {
    d |> dplyr::summarise(n = dplyr::n(), in_extent_pct = 100 * mean(in_extent),
                          in_core_pct = 100 * mean(in_extent_core),
                          fail_prob_pct = 100 * mean(!prob_ok), fail_depth_pct = 100 * mean(!depth_ok),
                          fail_geomorphic_pct = 100 * mean(!geo_ok),
                          class16_pct = 100 * mean(geo_class16), no_aca_geo_pct = 100 * mean(geo_no_aca),
                          in_only16_pct = 100 * mean(in_only16), in_nogeo_pct = 100 * mean(in_nogeo),
                          mean_agc_in = mean(AGC_pred[in_extent]),
                          mean_agc_out = mean(AGC_pred[!in_extent]), .groups = "drop")
  }
  a7 <- s7 |> dplyr::filter(is_agc)
  ext_tbl <- dplyr::bind_rows(
    sum_ext(a7) |> dplyr::mutate(level = "all", group = "all", .before = 1),
    a7 |> dplyr::group_by(group = ecoregion) |> sum_ext() |> dplyr::mutate(level = "ecoregion", .before = 1),
    a7 |> dplyr::group_by(group = as.character(year)) |> sum_ext() |> dplyr::mutate(level = "year", .before = 1)
  )
  print(as.data.frame(ext_tbl), row.names = FALSE, digits = 3)
  readr::write_csv(ext_tbl, out_file("P4_samples_in_extent.csv"))
}

# =============================================================================
# STEP 8: Presence probability threshold per year (apparent, in-sample F1)
# =============================================================================
step("STEP 8: PA threshold by year")
if (has_gee("P4_gee_pa_points.csv")) {
  pa <- read_gee("P4_gee_pa_points.csv") |>
    dplyr::filter(prob != NA_CODE, !is.na(prob)) |>
    dplyr::mutate(year = as.integer(year_pa))
  cat("PA cast check (original PA vs PA01):\n"); print(table(pa$PA, pa$PA01, useNA = "ifany"))

  thr_grid <- seq(0.05, 0.95, by = 0.01)
  f1_at <- function(obs, p, thr) {
    pred <- p >= thr
    tp <- sum(pred & obs == 1); fp <- sum(pred & obs == 0); fn <- sum(!pred & obs == 1)
    if (tp == 0) 0 else 2 * tp / (2 * tp + fp + fn)
  }
  thr_one <- function(d) {
    n1 <- sum(d$PA01 == 1); n0 <- sum(d$PA01 == 0)
    if (n1 == 0 || n0 == 0)
      return(tibble::tibble(n_pres = n1, n_abs = n0, best_thr = NA_real_, best_F1 = NA_real_,
                            F1_050 = NA_real_, F1_070 = NA_real_))
    f1 <- vapply(thr_grid, function(t) f1_at(d$PA01, d$prob, t), numeric(1))
    tibble::tibble(n_pres = n1, n_abs = n0, best_thr = thr_grid[which.max(f1)], best_F1 = max(f1),
                   F1_050 = f1_at(d$PA01, d$prob, 0.5), F1_070 = f1_at(d$PA01, d$prob, 0.7))
  }
  pa_thr <- dplyr::bind_rows(
    pa |> dplyr::group_by(year) |> dplyr::group_modify(~ thr_one(.x)) |> dplyr::ungroup() |>
      dplyr::mutate(year = as.character(year)),
    thr_one(pa) |> dplyr::mutate(year = "pooled", .before = 1)
  )
  print_tbl(pa_thr, 3)
  cat("Apparent F1: these points trained the PA model. Use the table for the threshold, not accuracy.\n")
  readr::write_csv(pa_thr, out_file("P4_pa_threshold_by_year.csv"))
}

# =============================================================================
# STEP 9: Distance from the extent to the nearest AGC sample
# =============================================================================
step("STEP 9: distance coverage")
if (has_gee("P4_gee_distance_area.csv")) {
  da <- read_core_addon("P4_gee_distance_area.csv", "P4_gee_distance_area_addon.csv",
                        c("eco_id", "dist_class")) |>
    dplyr::mutate(eco_id = as.integer(round(eco_id)), dist_class = as.integer(round(dist_class)),
                  dist_label = factor(DIST_LABELS[dist_class + 1], levels = DIST_LABELS),
                  design = DIST_DESIGN[dist_class + 1], ecoregion = eco_name(eco_id))
  agc_cols <- grep("^agc_g_\\d{4}$", names(da), value = TRUE)

  by_class <- da |>
    dplyr::group_by(dist_class, dist_label, design) |>
    dplyr::summarise(area_km2 = sum(area_m2) / 1e6,
                     dplyr::across(dplyr::all_of(agc_cols), sum), .groups = "drop") |>
    dplyr::mutate(area_pct = 100 * area_km2 / sum(area_km2), area_cum_pct = cumsum(area_pct))
  agc_share <- sweep(as.matrix(by_class[, agc_cols]), 2, colSums(by_class[, agc_cols]), "/")
  by_class$agc_pct_mean <- 100 * rowMeans(agc_share)
  by_class <- by_class |> dplyr::select(-dplyr::all_of(agc_cols))
  print(as.data.frame(by_class), row.names = FALSE, digits = 3)
  readr::write_csv(by_class, out_file("P4_distance_coverage_class.csv"))

  by_design <- by_class |>
    dplyr::group_by(design) |>
    dplyr::summarise(classes = paste(dist_label, collapse = ", "), area_km2 = sum(area_km2),
                     area_pct = sum(area_pct), agc_pct_mean = sum(agc_pct_mean), .groups = "drop") |>
    dplyr::mutate(design = factor(design, levels = CFG$designs)) |>
    dplyr::arrange(design)
  if (file.exists(out_file("P4_metrics_summary.csv"))) {
    ms <- readr::read_csv(out_file("P4_metrics_summary.csv"), show_col_types = FALSE) |>
      dplyr::filter(method %in% c("Pixel-based", "Tier 1")) |>
      dplyr::select(design, method, median_dist_m, MEC_mean, MEC_q025, MEC_q975) |>
      tidyr::pivot_wider(names_from = method, values_from = c(MEC_mean, MEC_q025, MEC_q975)) |>
      dplyr::mutate(design = factor(design, levels = CFG$designs))
    by_design <- dplyr::left_join(by_design, ms, by = "design")
  }
  print(as.data.frame(by_design), row.names = FALSE, digits = 3)
  readr::write_csv(by_design, out_file("P4_distance_coverage_design.csv"))

  # Distance class holding the p-th area percentile (replaces the slow exact
  # percentiles of P4_gee_distance_pct, which is optional)
  pct_class <- function(area, cls, p) {
    o <- order(cls); cum <- cumsum(area[o]) / sum(area)
    DIST_LABELS[cls[o][which(cum >= p)[1]] + 1]
  }
  eco_dist <- function(d) {
    d |> dplyr::summarise(area_km2 = sum(area_m2) / 1e6,
                          within_1km_pct = 100 * sum(area_m2[dist_class <= 5]) / sum(area_m2),
                          beyond_5km_pct = 100 * sum(area_m2[dist_class >= 8]) / sum(area_m2),
                          beyond_50km_pct = 100 * sum(area_m2[dist_class == 10]) / sum(area_m2),
                          p10_class = pct_class(area_m2, dist_class, 0.10),
                          median_class = pct_class(area_m2, dist_class, 0.50),
                          p90_class = pct_class(area_m2, dist_class, 0.90),
                          .groups = "drop")
  }
  by_eco <- dplyr::bind_rows(
    da |> dplyr::group_by(eco_id, ecoregion) |> eco_dist() |> dplyr::arrange(dplyr::desc(area_km2)),
    da |> eco_dist() |> dplyr::mutate(eco_id = -1L, ecoregion = "All", .before = 1)
  )
  print(as.data.frame(by_eco), row.names = FALSE, digits = 3)
  readr::write_csv(by_eco, out_file("P4_distance_coverage_eco.csv"))

  p_dist <- ggplot(by_class, aes(x = dist_label, y = area_pct, fill = design)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.1f", area_pct)), vjust = -0.3, size = 3) +
    labs(title = "Extent area by distance to the nearest AGC sample",
         x = NULL, y = "Share of extent area (%)", fill = "Validation design") +
    theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p_dist)
  # ggsave(out_file("P4_fig_distance_coverage.png"), p_dist, width = 8, height = 5, dpi = 300)
}
if (has_gee("P4_gee_distance_pct.csv")) {
  dp <- read_gee("P4_gee_distance_pct.csv") |>
    dplyr::mutate(eco_id = as.integer(round(eco_id)),
                  ecoregion = ifelse(eco_id == -1, "All", eco_name(eco_id)))
  cat("Distance percentiles (m; capped at 2 x search radius):\n")
  print(as.data.frame(dp), row.names = FALSE, digits = 4)
  readr::write_csv(dp, out_file("P4_distance_percentiles.csv"))
}

# =============================================================================
# STEP 10: Per-pixel variation of mapped AGC 2017-2024
# =============================================================================
step("STEP 10: AGC CV 2017-2024")
if (has_gee("P4_gee_cv_summary.csv")) {
  cv <- read_core_addon("P4_gee_cv_summary.csv", "P4_gee_cv_summary_addon.csv",
                        c("eco_id", "cv_class"), optional = TRUE) |>
    dplyr::mutate(eco_id = as.integer(round(eco_id)), cv_class = as.integer(round(cv_class)),
                  cv_label = factor(CV_LABELS[cv_class + 1], levels = CV_LABELS),
                  ecoregion = eco_name(eco_id))
  cv_nat <- cv |>
    dplyr::group_by(cv_label) |>
    dplyr::summarise(area_km2 = sum(area_m2) / 1e6, .groups = "drop") |>
    dplyr::mutate(area_pct = 100 * area_km2 / sum(area_km2))
  cv_eco <- cv |>
    dplyr::group_by(eco_id, ecoregion) |>
    dplyr::summarise(area_km2 = sum(area_m2) / 1e6, mean_cv_pct = sum(cv_x_area) / sum(area_m2),
                     cv_ge30_pct = 100 * sum(area_m2[cv_class >= 3]) / sum(area_m2), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(area_km2))
  cat(sprintf("Area-weighted mean CV (national): %.1f %%\n", sum(cv$cv_x_area) / sum(cv$area_m2)))
  print(as.data.frame(cv_nat), row.names = FALSE, digits = 3)
  print(as.data.frame(cv_eco), row.names = FALSE, digits = 3)
  readr::write_csv(dplyr::bind_rows(cv_nat |> dplyr::mutate(ecoregion = "All", .before = 1),
                                    cv |> dplyr::group_by(ecoregion, cv_label) |>
                                      dplyr::summarise(area_km2 = sum(area_m2) / 1e6, .groups = "drop") |>
                                      dplyr::group_by(ecoregion) |>
                                      dplyr::mutate(area_pct = 100 * area_km2 / sum(area_km2)) |>
                                      dplyr::ungroup()),
                   out_file("P4_cv_summary.csv"))
}

# =============================================================================
# STEP 11: GSE inter-annual similarity (dot product of unit vectors)
# =============================================================================
step("STEP 11: GSE similarity")
if (has_gee("P4_gee_gse_similarity_eco.csv")) {
  se_eco <- read_gee("P4_gee_gse_similarity_eco.csv") |>
    dplyr::mutate(eco_id = as.integer(round(eco_id)),
                  ecoregion = ifelse(eco_id == -1, "All", eco_name(eco_id)))
  print(as.data.frame(se_eco |> dplyr::select(ecoregion, band, p5, p25, p50, p75, p95, count)),
        row.names = FALSE, digits = 3)
  readr::write_csv(se_eco, out_file("P4_gse_similarity_eco.csv"))
}
if (has_gee("P4_gee_gse_sim_sample.csv")) {
  ss <- read_gee("P4_gee_gse_sim_sample.csv") |> dplyr::mutate(eco_id = as.integer(round(eco_id)))
  n_drawn <- nrow(ss)
  ss <- ss |> dplyr::filter(!is.na(sim_mean), !is.na(AGC_cv))   # 100 m cell centre outside the extent
  cat(sprintf("Similarity sample: %d points drawn, %d with values at the 10 m pixel\n", n_drawn, nrow(ss)))
  # weight = ecoregion extent area / pixels sampled there (stratified sample -> national).
  # P4_G04 drew the sample inside the core extent, so the weights use the core area
  # (fix 5 Nov 2026: core + add-on area was used before; medians changed < 0.001)
  if (!is.null(tot)) {
    ss <- ss |> dplyr::add_count(eco_id, name = "n_eco") |>
      dplyr::left_join(ECO_AREA_CORE, by = "eco_id") |>
      dplyr::mutate(w = area / n_eco)
    if (ADDON) {
      no_smp <- setdiff(unique(tot$eco_id), unique(ss$eco_id))
      cat("  NOTE: the similarity sample covers the core extent only (P4_G04).",
          if (length(no_smp) > 0) paste("No sample in:", paste(eco_name(no_smp), collapse = ", ")) else "",
          "\n")
    }
  } else {
    ss$w <- 1
    cat("  NOTE: no totals files, similarity percentiles are not area-weighted\n")
  }
  sim_cols <- grep("^sim_\\d{4}_\\d{4}$", names(ss), value = TRUE)
  pairs <- purrr::map_dfr(c(sim_cols, "sim_mean", "sim_min"), function(v) {
    q <- wquantile(ss[[v]], ss$w, c(0.05, 0.25, 0.5, 0.75, 0.95))
    tibble::tibble(pair = v, p05 = q[1], p25 = q[2], p50 = q[3], p75 = q[4], p95 = q[5],
                   below_0.9_pct = 100 * sum(ss$w[ss[[v]] < 0.9], na.rm = TRUE) / sum(ss$w))
  })
  print_tbl(pairs, 3)
  readr::write_csv(pairs, out_file("P4_gse_similarity_pairs.csv"))
  if ("sim_2023_2024" %in% pairs$pair) {
    oth <- pairs |> dplyr::filter(grepl("^sim_\\d{4}_\\d{4}$", pair), pair != "sim_2023_2024")
    cat(sprintf("Check for the 2024 CI: median similarity 2023-2024 = %.3f; other consecutive pairs %.3f-%.3f\n",
                pairs$p50[pairs$pair == "sim_2023_2024"], min(oth$p50), max(oth$p50)))
  }
  eco_sim <- ss |>
    dplyr::group_by(eco_id) |>
    dplyr::summarise(n = dplyr::n(),
                     sim_mean_p05 = quantile(sim_mean, 0.05, na.rm = TRUE),
                     sim_mean_median = median(sim_mean, na.rm = TRUE),
                     sim_mean_p95 = quantile(sim_mean, 0.95, na.rm = TRUE),
                     sim_min_median = median(sim_min, na.rm = TRUE),
                     agc_cv_median = median(AGC_cv, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(ecoregion = eco_name(eco_id), .after = eco_id)
  cat("GSE similarity per ecoregion (stratified sample):\n")
  print(as.data.frame(eco_sim), row.names = FALSE, digits = 3)
  readr::write_csv(eco_sim, out_file("P4_gse_similarity_eco_sample.csv"))
  rho <- suppressWarnings(cor(ss$sim_mean, ss$AGC_cv, method = "spearman", use = "complete.obs"))
  cat(sprintf("Spearman rho, GSE mean similarity vs AGC CV (sampled pixels): %.2f\n", rho))

  p_sim <- ggplot(ss, aes(x = sim_mean, y = AGC_cv)) +
    geom_bin2d(bins = 60) + scale_fill_viridis_c(trans = "log10") +
    labs(title = "Mapped AGC variation vs GSE change between years",
         x = "Mean GSE similarity of consecutive years (1 = identical)",
         y = "AGC CV 2017-2024 (%)", fill = "Pixels") +
    theme_minimal(base_size = 11)
  print(p_sim)
  # ggsave(out_file("P4_fig_gse_similarity_vs_cv.png"), p_sim, width = 6, height = 5, dpi = 300)
}
if (has_gee("P4_gee_gse_sim_points.csv")) {
  sp_raw <- read_gee("P4_gee_gse_sim_points.csv")
  check_rid(sp_raw, "P4_gee_gse_sim_points.csv")
  if ("sim_table" %in% names(sp_raw)) {
    # Identity test: table vector x image vector equals image x image (lag 0)
    # exactly when the uploaded GSE values are the image values
    lag0 <- vapply(seq_len(nrow(sp_raw)), function(i) sp_raw[[paste0("sim_", sp_raw$year[i])]][i], numeric(1))
    ok <- !is.na(sp_raw$sim_table) & sp_raw$sim_table != NA_CODE & !is.na(lag0) & lag0 != NA_CODE
    same <- abs(sp_raw$sim_table[ok] - lag0[ok]) < 1e-6
    cat(sprintf("Training-table GSE identical to the image GSE (survey year): %.1f%% of %d AGC points\n",
                100 * mean(same), sum(ok)))
  }
  sp <- sp_raw |>
    tidyr::pivot_longer(dplyr::matches("^sim_\\d{4}$"), names_to = "map_year", values_to = "sim") |>
    dplyr::mutate(map_year = as.integer(sub("sim_", "", map_year)),
                  lag = map_year - year, abs_lag = abs(lag)) |>
    dplyr::filter(sim != NA_CODE, !is.na(sim))
  lag_tbl <- sp |>
    dplyr::group_by(abs_lag) |>
    dplyr::summarise(n = dplyr::n(), p05 = quantile(sim, 0.05), p25 = quantile(sim, 0.25),
                     median = median(sim), p75 = quantile(sim, 0.75),
                     below_0.9_pct = 100 * mean(sim < 0.9), .groups = "drop")
  print_tbl(lag_tbl, 3)
  if (lag_tbl$median[lag_tbl$abs_lag == 0] < 0.999)
    cat("WARNING: same-year similarity is below 0.999: the table GSE and the image differ.\n")
  readr::write_csv(lag_tbl, out_file("P4_gse_similarity_points_lag.csv"))

  p_lag <- ggplot(sp |> dplyr::filter(abs_lag > 0), aes(x = factor(abs_lag), y = sim)) +
    geom_boxplot(outlier.size = 0.3, fill = "grey90") +
    labs(title = "GSE at the AGC sample pixels: survey year vs other years",
         x = "Years between survey and map year", y = "GSE similarity (dot product)") +
    theme_minimal(base_size = 11)
  print(p_lag)
  # ggsave(out_file("P4_fig_gse_similarity_lag.png"), p_lag, width = 6, height = 4, dpi = 300)
}

if (has_gee("P4_gee_importance.csv")) {
  imp <- read_gee("P4_gee_importance.csv") |>
    dplyr::group_by(model) |> dplyr::slice_max(importance, n = 5) |> dplyr::ungroup()
  cat("Top-5 GSE bands by RF importance (GEE):\n"); print(as.data.frame(imp), row.names = FALSE)
}

# =============================================================================
# STEP 12: Extent under three geomorphic masks (optional, P4_G05)
# current = class 16 and pixels without ACA geomorphic removed (09b, P4_G01)
# only16  = only class 16 removed; nogeo = no geomorphic mask
# =============================================================================
step("STEP 12: extent sensitivity to the geomorphic mask")
es_files <- list.files(gee_dir, pattern = "^P4_gee_extent_sensitivity(_\\d+m)?\\.csv$", full.names = TRUE)
if (length(es_files) == 0) {
  cat("  -> skipped: no P4_gee_extent_sensitivity*.csv\n")
} else {
  es_file <- es_files[order(nchar(basename(es_files)))][1]     # prefer the 10 m file if both exist
  cat("Using", basename(es_file), "\n")
  es <- readr::read_csv(es_file, show_col_types = FALSE) |>
    dplyr::mutate(eco_id = as.integer(round(eco_id)), ecoregion = eco_name(eco_id))
  es_eco <- es |>
    dplyr::transmute(eco_id, ecoregion, current_km2 = a_current_m2 / 1e6,
                     only16_km2 = a_only16_m2 / 1e6, nogeo_km2 = a_nogeo_m2 / 1e6) |>
    dplyr::arrange(dplyr::desc(nogeo_km2))
  nat_es <- es_eco |>
    dplyr::summarise(dplyr::across(dplyr::ends_with("_km2"), sum)) |>
    tidyr::pivot_longer(dplyr::everything(), names_to = "mask", values_to = "area_km2") |>
    dplyr::mutate(mask = sub("_km2", "", mask),
                  ratio_to_current = area_km2 / area_km2[mask == "current"])
  if (!is.null(tot)) {
    # scale to the 10 m extent of P4_G02 (removes the effect of a coarser reduce scale).
    # "current" = the core extent, so the ratios are scaled to the core area, not to
    # core + add-on (fix 5 Nov 2026: before, AREA_KM2 was used, which includes the add-on)
    nat_es <- nat_es |> dplyr::mutate(area_km2_at_10m = ratio_to_current * AREA_CORE_KM2)
    cat(sprintf("Check: current %.2f km2 here vs %.2f km2 core extent at 10 m (P4_G02)\n",
                nat_es$area_km2[nat_es$mask == "current"], AREA_CORE_KM2))
    if (ADDON)
      cat(sprintf("Check: only16 scaled to 10 m %.2f km2 vs %.2f km2 mapped (core + add-on, final extent)\n",
                  nat_es$area_km2_at_10m[nat_es$mask == "only16"], AREA_KM2))
  }
  nat_es <- nat_es |>
    dplyr::mutate(tier1_pooled_tC = mean(agc$AGC_pred) *
                    (if ("area_km2_at_10m" %in% names(nat_es)) area_km2_at_10m else area_km2))
  print(as.data.frame(nat_es), row.names = FALSE, digits = 4)
  print(as.data.frame(es_eco), row.names = FALSE, digits = 4)
  readr::write_csv(nat_es, out_file("P4_extent_sensitivity_national.csv"))
  readr::write_csv(es_eco, out_file("P4_extent_sensitivity_ecoregion.csv"))
}

cat("\nP4_07_gee_results.R COMPLETE\n")
