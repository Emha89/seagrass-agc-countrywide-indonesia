# =============================================================================
# P4_05_tier_uncertainty.R
# Paper 4 revised design: uncertainty terms, Tier 1 national totals, Tier 2
# class means and their sensitivity checks.
#
# Reads : P4_samples.csv, P4_oof_agc_cluster.csv, P4_oof_agc_pixel.csv,
#         P4_best_params.csv
#         optional: gee_cover_class_area_by_year.csv (columns year, class, area_km2)
# Writes: P4_uncertainty_params.csv     ICC, DEFF, N_eff, SD_model, SD_data, bias SE
#         P4_sd_model_by_year.csv
#         P4_tier1_national.csv         per year and pooled, DEFF-corrected CI
#         P4_tier2_class_by_year.csv    class x year means (flag n < CFG$min_n_cell)
#         P4_tier2_class_means.csv      pooled class means (provisional classes)
#         P4_tier2_loyo.csv             class means leaving one survey year out
#         P4_tier2_national.csv         only when the GEE class-area file exists
#         P4_year_composition.csv
#
# Tier 2 classes here come from an R cover model fitted on all cover samples
# and predicted at the AGC sample pixels. They are PROVISIONAL: recompute with
# classes sampled from the GEE cover map, so class means and class areas come
# from the same map.
# =============================================================================

# setwd("path/to/seagrass-agc-countrywide-indonesia/R")   # folder that holds these scripts
source("P4_00_config_functions.R")

samples <- read_samples()
agc <- samples |> dplyr::filter(is_agc)
oof_cl  <- read_oof("agc", designs = "cluster")
oof_pix <- read_oof("agc", designs = "pixel")

# =============================================================================
# STEP 1: Design effect from out-of-fold residuals grouped by DBSCAN cluster
# (clustered points only, as in the thesis definition)
# =============================================================================
icc_from <- function(oof, label) {
  oof |>
    dplyr::filter(grepl("^c", g_cluster)) |>
    dplyr::mutate(resid = obs - pred_pixel) |>
    dplyr::group_by(rep) |>
    dplyr::group_modify(~ icc_oneway(.x$resid, .x$g_cluster)) |>
    dplyr::ungroup() |>
    dplyr::summarise(dplyr::across(c(N, k, n0, icc, deff), mean)) |>
    dplyr::mutate(source = label, .before = 1)
}
icc_tbl <- dplyr::bind_rows(icc_from(oof_cl, "cluster design (used)"),
                            icc_from(oof_pix, "pixel design (sensitivity)"))
print_tbl(icc_tbl)
DEFF <- icc_tbl$deff[1]
cat(sprintf("Design effect used: %.2f | N_eff of all AGC samples: %.0f\n", DEFF, nrow(agc) / DEFF))

# =============================================================================
# STEP 2: SD_model (cluster design RMSE), SD_data (RMS of reference SD)
# =============================================================================
sd_model <- oof_cl |>
  dplyr::group_by(rep) |>
  dplyr::summarise(RMSE = sqrt(mean((obs - pred_pixel)^2)), .groups = "drop") |>
  dplyr::summarise(SD_model = mean(RMSE)) |>
  dplyr::pull(SD_model)

sd_model_year <- oof_cl |>
  dplyr::group_by(rep, year) |>
  dplyr::summarise(n = dplyr::n(), RMSE = sqrt(mean((obs - pred_pixel)^2)), .groups = "drop") |>
  dplyr::group_by(year) |>
  dplyr::summarise(n = mean(n), SD_model = mean(RMSE), .groups = "drop") |>
  dplyr::mutate(use_year_value = n >= 30)
readr::write_csv(sd_model_year, out_file("P4_sd_model_by_year.csv"))
print_tbl(sd_model_year)

sd_data <- sqrt(mean(agc$sd_ref^2))

# Diagnostic: bias of the mean prediction with a cluster-robust SE.
# This is the error that does not shrink with the number of pixels.
bias_tbl <- oof_cl |>
  dplyr::mutate(resid = obs - pred_pixel) |>
  dplyr::group_by(rep) |>
  dplyr::group_modify(~ cluster_se_mean(.x$resid, .x$g_cluster)) |>
  dplyr::ungroup() |>
  dplyr::summarise(dplyr::across(c(mean_resid, se_cluster, n_groups), mean))

unc <- tibble::tribble(
  ~parameter, ~value, ~note,
  "ICC", icc_tbl$icc[1], "out-of-fold residuals, cluster design, clustered points only",
  "adjusted_mean_cluster_size", icc_tbl$n0[1], "n0 of the one-way ANOVA",
  "DEFF", DEFF, "1 + (n0 - 1) x ICC; applied to Tier 1, Tier 2 and pixel-based CIs",
  "N_agc", nrow(agc), "AGC samples",
  "N_eff", nrow(agc) / DEFF, "N_agc / DEFF",
  "SD_model_gC_m2", sd_model, "RMSE of the cluster design (between meadows)",
  "SD_data_gC_m2", sd_data, "root mean square of (AGC_up - AGC_low) / 3.92",
  "SD_total_gC_m2", sqrt(sd_model^2 + sd_data^2), "sqrt(SD_model^2 + SD_data^2)",
  "mean_resid_gC_m2", bias_tbl$mean_resid, "diagnostic: mean (obs - pred), cluster design",
  "se_mean_resid_gC_m2", bias_tbl$se_cluster, "diagnostic: cluster-robust SE of the mean residual",
  "bias_halfwidth_national_tC", 1.96 * bias_tbl$se_cluster * CFG$extent_km2,
  "diagnostic: 1.96 x SE x extent (km2) = +/- on the national pixel-based total"
)
print_tbl(unc)
readr::write_csv(unc, out_file("P4_uncertainty_params.csv"))

# =============================================================================
# STEP 3: Tier 1 national totals (per survey year and pooled), DEFF-corrected CI
# =============================================================================
tier1 <- dplyr::bind_rows(
  agc |> dplyr::group_by(year) |>
    dplyr::group_modify(~ tier1_summary(.x$AGC_pred, DEFF, CFG$extent_km2)) |>
    dplyr::ungroup() |> dplyr::mutate(year = as.character(year)),
  tier1_summary(agc$AGC_pred, DEFF, CFG$extent_km2) |> dplyr::mutate(year = "pooled")
) |>
  dplyr::mutate(ci_note = ifelse(is.na(ci_half_gC_m2), "n_eff < 2, CI not defined", ""))
print_tbl(tier1, 2)
readr::write_csv(tier1, out_file("P4_tier1_national.csv"))

# =============================================================================
# STEP 4: Provisional cover classes at AGC sample pixels (final cover model)
# =============================================================================
p_cov <- read_best_params("COVER")
cov_all <- samples |> dplyr::filter(is_cov)
m_cov_final <- fit_rf(cov_all, CFG$cover_response, p_cov, CFG$seed)
agc$cov_pred  <- predict_rf(m_cov_final, agc)
agc$cov_class <- cover_class(agc$cov_pred)
cat("\nAGC samples per provisional cover class:\n")
print(table(agc$cov_class, useNA = "ifany"))

# =============================================================================
# STEP 5: Class x year table (a check of survey composition, not of time:
# survey year and region are confounded in this dataset)
# =============================================================================
class_year <- agc |>
  dplyr::group_by(cov_class, year) |>
  dplyr::summarise(n = dplyr::n(), mean_gC_m2 = mean(AGC_pred), sd_gC_m2 = sd(AGC_pred),
                   .groups = "drop") |>
  dplyr::mutate(reported = n >= CFG$min_n_cell)
readr::write_csv(class_year, out_file("P4_tier2_class_by_year.csv"))
print_tbl(class_year, 2)

spread <- class_year |>
  dplyr::filter(reported) |>
  dplyr::group_by(cov_class) |>
  dplyr::summarise(n_years = dplyr::n(), min_mean = min(mean_gC_m2), max_mean = max(mean_gC_m2),
                   range = max_mean - min_mean, .groups = "drop")
t1_years <- tier1 |> dplyr::filter(year != "pooled", n >= CFG$min_n_cell)
cat(sprintf("\nRange of Tier 1 annual means (n >= %d): %.1f to %.1f gC m-2\n",
            CFG$min_n_cell, min(t1_years$mean_gC_m2), max(t1_years$mean_gC_m2)))
cat("Range of class means across survey years:\n")
print_tbl(spread, 2)

# =============================================================================
# STEP 6: Pooled class means and leave-one-year-out sensitivity
# =============================================================================
class_means <- function(d) {
  d |>
    dplyr::group_by(cov_class, .drop = FALSE) |>
    dplyr::summarise(n = dplyr::n(), mean_gC_m2 = mean(AGC_pred), sd_gC_m2 = sd(AGC_pred),
                     .groups = "drop") |>
    dplyr::mutate(n_eff = n / DEFF)
}
cm_pooled <- class_means(agc)
readr::write_csv(cm_pooled, out_file("P4_tier2_class_means.csv"))
cat("\nPooled class means (provisional):\n")
print_tbl(cm_pooled, 2)

loyo <- purrr::map_dfr(sort(unique(agc$year)), function(y) {
  class_means(agc |> dplyr::filter(year != y)) |> dplyr::mutate(left_out_year = y, .before = 1)
})
readr::write_csv(loyo, out_file("P4_tier2_loyo.csv"))
print(as.data.frame(loyo |> dplyr::select(left_out_year, cov_class, n, mean_gC_m2) |>
                      tidyr::pivot_wider(names_from = cov_class, values_from = c(n, mean_gC_m2))),
      row.names = FALSE)

# =============================================================================
# STEP 7: Tier 2 national totals (needs class areas from the GEE cover map)
# File: <out_dir>/gee_cover_class_area_by_year.csv with columns year, class, area_km2
# class labels must match CFG$cover_labels
# =============================================================================
tier2_total <- function(areas, cm, overall_mean) {
  d <- dplyr::left_join(areas, cm, by = c("class" = "cov_class"))
  empty <- is.na(d$mean_gC_m2) | d$n == 0
  d$mean_gC_m2[empty] <- overall_mean
  d$var_term <- ifelse(empty | d$n_eff < 2, NA_real_, d$area_km2^2 * d$sd_gC_m2^2 / d$n_eff)
  tibble::tibble(total_tC = sum(d$mean_gC_m2 * d$area_km2),
                 ci_half_tC = 1.96 * sqrt(sum(d$var_term, na.rm = TRUE)),
                 classes_without_samples = paste(d$class[empty], collapse = ", "))
}

area_file <- out_file("gee_cover_class_area_by_year.csv")
if (file.exists(area_file)) {
  areas <- readr::read_csv(area_file, show_col_types = FALSE) |>
    dplyr::mutate(class = as.character(class))
  cm_p <- cm_pooled |> dplyr::mutate(cov_class = as.character(cov_class))
  t2 <- areas |>
    dplyr::group_by(year) |>
    dplyr::group_modify(~ tier2_total(.x, cm_p, mean(agc$AGC_pred))) |>
    dplyr::ungroup()
  t2_loyo <- purrr::map_dfr(unique(loyo$left_out_year), function(y) {
    cm_y <- loyo |> dplyr::filter(left_out_year == y) |> dplyr::mutate(cov_class = as.character(cov_class))
    areas |> dplyr::group_by(year) |>
      dplyr::group_modify(~ tier2_total(.x, cm_y, mean(agc$AGC_pred[agc$year != y]))) |>
      dplyr::ungroup() |> dplyr::mutate(left_out_year = y)
  }) |>
    dplyr::group_by(year) |>
    dplyr::summarise(loyo_min_tC = min(total_tC), loyo_max_tC = max(total_tC), .groups = "drop")
  t2 <- dplyr::left_join(t2, t2_loyo, by = "year") |>
    dplyr::mutate(ci_low_tC = total_tC - ci_half_tC, ci_up_tC = total_tC + ci_half_tC)
  print_tbl(t2, 0)
  readr::write_csv(t2, out_file("P4_tier2_national.csv"))
} else {
  cat(sprintf("\nTier 2 totals skipped: %s not found (export class areas from GEE first).\n",
              basename(area_file)))
}

# =============================================================================
# STEP 8: Survey composition by year and ecoregion (for the time-series figure)
# =============================================================================
comp <- agc |>
  dplyr::count(year, ecoregion, name = "n") |>
  dplyr::group_by(year) |>
  dplyr::mutate(share = n / sum(n)) |>
  dplyr::ungroup()
readr::write_csv(comp, out_file("P4_year_composition.csv"))

print(
  ggplot(comp, aes(x = factor(year), y = share, fill = ecoregion)) +
    geom_col(colour = "white", linewidth = 0.2) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Survey composition behind each Tier 1 annual mean",
         x = "Survey year", y = "Share of AGC samples", fill = "Ecoregion") +
    theme_minimal(base_size = 11)
)

cat("\nP4_05_tier_uncertainty.R COMPLETE\n")
