# =============================================================================
# P4_04_metrics.R
# Paper 4 revised design: accuracy metrics, distance curve, repeat stability,
# error by AGC bin, cover-map accuracy and the shrinkage check.
#
# Reads : P4_oof_agc_<design>.csv, P4_oof_cov_<design>.csv, P4_samples.csv
# Writes: P4_metrics_by_repeat.csv     metrics per design x repeat x method
#         P4_metrics_summary.csv       mean and 2.5-97.5 percentiles + median distance
#         P4_metrics_ecoregion.csv     per-ecoregion metrics (LOLO), n flag
#         P4_lolo_region_means.csv     regional mean AGC per method when the ecoregion is held out
#         P4_repeat_stability.csv      cumulative mean MEC by repeat
#         P4_error_by_agc_bin.csv      absolute error per 10 gC m-2 bin (pixel-based)
#         P4_cover_metrics.csv, P4_cover_confusion.csv, P4_cover_error_by_bin.csv
#         P4_shrinkage_groups.csv, P4_shrinkage_fit.csv
#         P4_shrinkage_roi.csv (only when CFG$roi_path is set)
# All plots are printed; ggsave lines are commented.
# =============================================================================

# setwd("path/to/seagrass-agc-countrywide-indonesia/R")   # folder that holds these scripts
source("P4_00_config_functions.R")

oof <- read_oof("agc")
cov <- read_oof("cov")
samples <- read_samples()

method_cols   <- c("Pixel-based" = "pred_pixel", "Tier 1" = "pred_t1", "Tier 2" = "pred_t2")
design_levels <- CFG$designs
q025 <- function(x) quantile(x, 0.025, na.rm = TRUE)
q975 <- function(x) quantile(x, 0.975, na.rm = TRUE)

long <- oof |>
  tidyr::pivot_longer(dplyr::all_of(unname(method_cols)), names_to = "col", values_to = "pred") |>
  dplyr::mutate(method = names(method_cols)[match(col, method_cols)],
                design = factor(design, levels = design_levels))

# =============================================================================
# STEP 1: Metrics per design x repeat x method (pooled out-of-fold predictions)
# =============================================================================
m_rep <- long |>
  dplyr::group_by(design, rep, method) |>
  dplyr::summarise(reg_metrics(obs, pred), .groups = "drop")

dist_rep <- oof |>
  dplyr::group_by(design, rep) |>
  dplyr::summarise(median_dist_m = median(dist_m), .groups = "drop") |>
  dplyr::mutate(design = factor(design, levels = design_levels))

m_rep <- dplyr::left_join(m_rep, dist_rep, by = c("design", "rep"))
readr::write_csv(m_rep, out_file("P4_metrics_by_repeat.csv"))

m_sum <- m_rep |>
  dplyr::group_by(design, method) |>
  dplyr::summarise(
    n_repeats = dplyr::n(), n = mean(n),
    median_dist_m = mean(median_dist_m),
    dplyr::across(c(RMSE, MAE, R2, MEC, bias_pct),
                  list(mean = ~ mean(.x, na.rm = TRUE), q025 = q025, q975 = q975)),
    .groups = "drop") |>
  dplyr::arrange(design, method)
print_tbl(m_sum |> dplyr::select(design, method, median_dist_m, RMSE_mean, MEC_mean,
                                 MEC_q025, MEC_q975, bias_pct_mean))
readr::write_csv(m_sum, out_file("P4_metrics_summary.csv"))

t2_fb <- oof |>
  dplyr::group_by(design) |>
  dplyr::summarise(pct_tier2_fallback = 100 * mean(t2_fallback), .groups = "drop")
cat("\nTier 2 test samples whose class had no training sample (fell back to Tier 1):\n")
print(as.data.frame(t2_fb), row.names = FALSE)

# =============================================================================
# STEP 2: Leave-one-ecoregion-out, per ecoregion
# =============================================================================
eco_long <- long |> dplyr::filter(design == "ecoregion")
if (nrow(eco_long) > 0) {
  m_eco <- eco_long |>
    dplyr::group_by(rep, ecoregion, method) |>
    dplyr::summarise(reg_metrics(obs, pred), median_dist_km = median(dist_m) / 1000,
                     .groups = "drop") |>
    dplyr::group_by(ecoregion, method) |>
    dplyr::summarise(dplyr::across(c(n, RMSE, MAE, R2, MEC, bias_pct, median_dist_km), mean),
                     .groups = "drop") |>
    dplyr::mutate(reported = n >= CFG$min_n_eco) |>
    dplyr::arrange(dplyr::desc(n), method)
  readr::write_csv(m_eco, out_file("P4_metrics_ecoregion.csv"))
  print_tbl(m_eco |> dplyr::select(ecoregion, method, n, MEC, bias_pct, reported))

  # One region with tiny variance (MEC far below zero) dominates a mean, so medians are reported too
  eq_w <- m_eco |>
    dplyr::filter(reported) |>
    dplyr::group_by(method) |>
    dplyr::summarise(n_ecoregions = dplyr::n(),
                     MEC_equal_weight = mean(MEC), MEC_median = median(MEC),
                     RMSE_equal_weight = mean(RMSE),
                     abs_bias_pct_equal_weight = mean(abs(bias_pct)),
                     abs_bias_pct_median = median(abs(bias_pct)), .groups = "drop")
  cat(sprintf("\nOver ecoregions with n >= %d (equal weight and median):\n", CFG$min_n_eco))
  print_tbl(eq_w)
  readr::write_csv(eq_w, out_file("P4_metrics_ecoregion_equal_weight.csv"))

  # Regional mean AGC: what matters for an area total
  reg_means <- oof |>
    dplyr::filter(design == "ecoregion") |>
    dplyr::group_by(rep, ecoregion) |>
    dplyr::summarise(n = dplyr::n(), field = mean(obs), pixel_based = mean(pred_pixel),
                     tier1 = mean(pred_t1), tier2 = mean(pred_t2), .groups = "drop") |>
    dplyr::group_by(ecoregion) |>
    dplyr::summarise(dplyr::across(c(n, field, pixel_based, tier1, tier2), mean), .groups = "drop") |>
    dplyr::mutate(bias_pct_pixel = 100 * (pixel_based - field) / field,
                  bias_pct_tier1 = 100 * (tier1 - field) / field,
                  bias_pct_tier2 = 100 * (tier2 - field) / field,
                  pixel_closer_than_tier1 = abs(pixel_based - field) < abs(tier1 - field),
                  reported = n >= CFG$min_n_eco) |>
    dplyr::arrange(dplyr::desc(n))
  readr::write_csv(reg_means, out_file("P4_lolo_region_means.csv"))
  cat("\nRegional mean AGC when the ecoregion is held out (gC m-2):\n")
  print_tbl(reg_means, 1)
  rr <- dplyr::filter(reg_means, reported)
  cat(sprintf("Pixel-based closer to the field mean than Tier 1 in %d of %d ecoregions\n",
              sum(rr$pixel_closer_than_tier1), nrow(rr)))
}

# =============================================================================
# STEP 3: Stability of repeats (cumulative mean MEC)
# =============================================================================
stab <- m_rep |>
  dplyr::filter(design != "ecoregion") |>
  dplyr::arrange(design, method, rep) |>
  dplyr::group_by(design, method) |>
  dplyr::mutate(cum_mean_MEC = cumsum(MEC) / dplyr::row_number()) |>
  dplyr::ungroup()
readr::write_csv(stab, out_file("P4_repeat_stability.csv"))

stab_last <- stab |>
  dplyr::group_by(design, method) |>
  dplyr::summarise(last_change = abs(dplyr::last(cum_mean_MEC) - dplyr::nth(cum_mean_MEC, -2)),
                   stable = last_change < 0.01, .groups = "drop")
cat("\nChange in cumulative mean MEC at the last repeat (stable if < 0.01):\n")
print_tbl(stab_last, 4)
if (any(!stab_last$stable))
  cat("Some designs are not yet stable: raise CFG$n_repeats and rerun P4_03 for those designs.\n")

print(
  ggplot(stab, aes(x = rep, y = cum_mean_MEC, colour = method)) +
    geom_line() + geom_point(size = 1) +
    facet_wrap(~ design, nrow = 1) +
    labs(title = "Cumulative mean MEC by number of repeats", x = "Repeat", y = "Cumulative mean MEC",
         colour = NULL) +
    theme_minimal(base_size = 11)
)

# =============================================================================
# STEP 4: MEC against median distance to the nearest training sample
# =============================================================================
p_dist <- ggplot(m_sum, aes(x = median_dist_m, y = MEC_mean, colour = method, group = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line() +
  geom_errorbar(aes(ymin = MEC_q025, ymax = MEC_q975), width = 0.04) +
  geom_point(size = 2) +
  geom_text(data = dplyr::filter(m_sum, method == "Pixel-based"),
            aes(label = design), vjust = -1.1, size = 3, colour = "grey20", show.legend = FALSE) +
  scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  labs(title = "Accuracy against distance from field data",
       x = "Median distance from test to nearest training sample (m, log scale)",
       y = "MEC (mean, 2.5-97.5 percentiles across repeats)", colour = NULL) +
  theme_minimal(base_size = 11)
print(p_dist)
# ggsave(out_file("P4_fig_mec_distance.png"), p_dist, width = 8, height = 5, dpi = 300)

# =============================================================================
# STEP 5: Pixel-based absolute error per 10 gC m-2 bin (out-of-fold)
# =============================================================================
bin_breaks <- seq(0, ceiling(max(oof$obs) / 10) * 10, by = 10)
err_bin <- oof |>
  dplyr::filter(design %in% c("pixel", "cluster")) |>
  dplyr::mutate(abs_err = abs(obs - pred_pixel),
                bin = cut(obs, breaks = bin_breaks, right = FALSE, include.lowest = TRUE)) |>
  dplyr::group_by(design, bin) |>
  dplyr::summarise(n_per_repeat = dplyr::n() / dplyr::n_distinct(rep),
                   median_abs_err = median(abs_err),
                   q25 = quantile(abs_err, 0.25), q75 = quantile(abs_err, 0.75),
                   mean_pred = mean(pred_pixel), .groups = "drop")
readr::write_csv(err_bin, out_file("P4_error_by_agc_bin.csv"))
print_tbl(err_bin, 2)

print(
  ggplot(err_bin, aes(x = bin, y = median_abs_err)) +
    geom_pointrange(aes(ymin = q25, ymax = q75)) +
    facet_wrap(~ design, ncol = 1) +
    labs(title = "Pixel-based absolute error by observed AGC (out-of-fold)",
         x = "Observed AGC bin (gC m-2)", y = "Absolute error (gC m-2), median and IQR") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

# =============================================================================
# STEP 6: Cover model accuracy (out-of-fold) and class confusion
# =============================================================================
if (nrow(cov) > 0) {
  cov <- cov |> dplyr::mutate(design = factor(design, levels = design_levels))
  cov_rep <- cov |>
    dplyr::filter(obs > 0) |>
    dplyr::group_by(design, rep) |>
    dplyr::summarise(reg_metrics(obs, pred), .groups = "drop")
  cov_sum <- cov_rep |>
    dplyr::group_by(design) |>
    dplyr::summarise(dplyr::across(c(n, RMSE, MAE, R2, MEC), ~ mean(.x, na.rm = TRUE)),
                     MEC_q025 = q025(MEC), MEC_q975 = q975(MEC), .groups = "drop")
  readr::write_csv(cov_sum, out_file("P4_cover_metrics.csv"))
  cat("\nCover model (samples with seagrass present):\n")
  print_tbl(cov_sum)

  conf <- cov |>
    dplyr::filter(obs > 0) |>
    dplyr::mutate(obs_class = cover_class(obs), pred_class = cover_class(pred)) |>
    dplyr::count(design, obs_class, pred_class, name = "n")
  readr::write_csv(conf, out_file("P4_cover_confusion.csv"))

  conf_stats <- conf |>
    dplyr::group_by(design) |>
    dplyr::group_modify(function(d, key) {
      tab <- xtabs(n ~ obs_class + pred_class, data = d)
      N <- sum(tab); po <- sum(diag(tab)) / N
      pe <- sum(rowSums(tab) * colSums(tab)) / N^2
      tibble::tibble(overall_accuracy = po, kappa = (po - pe) / (1 - pe),
                     producer_acc = paste(sprintf("%s=%.2f", rownames(tab), diag(tab) / rowSums(tab)), collapse = "; "),
                     user_acc = paste(sprintf("%s=%.2f", colnames(tab), diag(tab) / colSums(tab)), collapse = "; "))
    }) |>
    dplyr::ungroup()
  cat("\nCover class accuracy (KepMenLH classes):\n")
  print(as.data.frame(conf_stats), row.names = FALSE)
  readr::write_csv(conf_stats, out_file("P4_cover_class_accuracy.csv"))

  cov_bin <- cov |>
    dplyr::filter(obs > 0, design %in% c("pixel", "cluster")) |>
    dplyr::mutate(bin = cut(obs, breaks = seq(0, 100, 10), right = FALSE, include.lowest = TRUE),
                  err = pred - obs) |>
    dplyr::group_by(design, bin) |>
    dplyr::summarise(n_per_repeat = dplyr::n() / dplyr::n_distinct(rep),
                     mean_err = mean(err), median_abs_err = median(abs(err)), .groups = "drop")
  readr::write_csv(cov_bin, out_file("P4_cover_error_by_bin.csv"))
  print_tbl(cov_bin, 2)
}

# =============================================================================
# STEP 7: Shrinkage toward the training mean (leak-free)
# Group mean of observed AGC vs group mean of out-of-fold pixel-based predictions
# =============================================================================
national_mean <- mean(samples$AGC_pred[samples$is_agc])
per_sample <- oof |>
  dplyr::group_by(design, gee_id, g_cluster, ecoregion, obs) |>
  dplyr::summarise(pred = mean(pred_pixel), .groups = "drop")

shr_groups <- dplyr::bind_rows(
  per_sample |> dplyr::filter(design == "cluster", grepl("^c", g_cluster)) |>
    dplyr::group_by(level = "DBSCAN cluster (cluster design)", group = g_cluster) |>
    dplyr::summarise(n = dplyr::n(), ecoregion = dplyr::first(ecoregion),
                     mean_obs = mean(obs), mean_pred = mean(pred), .groups = "drop") |>
    dplyr::filter(n >= 10),
  per_sample |> dplyr::filter(design == "ecoregion") |>
    dplyr::group_by(level = "Ecoregion (LOLO design)", group = ecoregion) |>
    dplyr::summarise(n = dplyr::n(), ecoregion = dplyr::first(ecoregion),
                     mean_obs = mean(obs), mean_pred = mean(pred), .groups = "drop") |>
    dplyr::filter(n >= CFG$min_n_eco)
)
readr::write_csv(shr_groups, out_file("P4_shrinkage_groups.csv"))

shr_fit <- shr_groups |>
  dplyr::group_by(level) |>
  dplyr::group_modify(function(d, key) {
    if (nrow(d) < 3) return(tibble::tibble(n_groups = nrow(d)))
    f <- lm(mean_pred ~ mean_obs, data = d)
    b <- coef(f)
    tibble::tibble(n_groups = nrow(d), slope = b[[2]], intercept = b[[1]],
                   crossover_gC_m2 = b[[1]] / (1 - b[[2]]),
                   r = cor(d$mean_obs, d$mean_pred),
                   spearman_obs_vs_pctdiff = cor(d$mean_obs, (d$mean_pred - d$mean_obs) / d$mean_obs,
                                                 method = "spearman"))
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(national_training_mean = national_mean)
cat("\nShrinkage fit (slope < 1 = predictions compressed toward the training mean):\n")
print_tbl(shr_fit)
readr::write_csv(shr_fit, out_file("P4_shrinkage_fit.csv"))

if (nrow(shr_groups) > 0) {
  print(
    ggplot(shr_groups, aes(x = mean_obs, y = mean_pred)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
      geom_hline(yintercept = national_mean, linetype = "dotted", colour = "firebrick") +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "steelblue") +
      geom_point(aes(size = n), alpha = 0.7) +
      facet_wrap(~ level) +
      labs(title = "Group mean AGC: field vs out-of-fold pixel-based prediction",
           subtitle = "Dashed = 1:1, dotted = national training mean",
           x = "Mean field AGC (gC m-2)", y = "Mean predicted AGC (gC m-2)", size = "n") +
      theme_minimal(base_size = 11)
  )
}

# =============================================================================
# STEP 8 (optional): ROI shrinkage check, leak-free (cluster design)
# =============================================================================
if (!is.na(CFG$roi_path) && file.exists(CFG$roi_path)) {
  sf::sf_use_s2(FALSE)
  roi <- sf::st_transform(sf::st_read(CFG$roi_path, quiet = TRUE), 4326)
  pts <- sf::st_as_sf(as.data.frame(samples[samples$is_agc, c("gee_id", "year", "xcoord", "ycoord")]),
                      coords = c("xcoord", "ycoord"), crs = 4326)
  j <- suppressMessages(sf::st_join(pts, roi[, CFG$roi_id_field], join = sf::st_within)) |>
    sf::st_drop_geometry() |>
    dplyr::filter(!is.na(.data[[CFG$roi_id_field]])) |>
    dplyr::rename(loc_id = dplyr::all_of(CFG$roi_id_field))
  roi_year <- j |> dplyr::count(loc_id, year) |> dplyr::group_by(loc_id) |>
    dplyr::slice_max(n, n = 1, with_ties = FALSE) |> dplyr::ungroup()
  roi_tbl <- j |>
    dplyr::semi_join(roi_year, by = c("loc_id", "year")) |>
    dplyr::inner_join(per_sample |> dplyr::filter(design == "cluster"), by = "gee_id") |>
    dplyr::group_by(loc_id, year) |>
    dplyr::summarise(n = dplyr::n(), mean_obs = mean(obs), mean_pred_oof = mean(pred),
                     pct_diff = 100 * (mean_pred_oof - mean_obs) / mean_obs, .groups = "drop")
  readr::write_csv(roi_tbl, out_file("P4_shrinkage_roi.csv"))
  print_tbl(roi_tbl, 2)
  if (nrow(roi_tbl) >= 3) {
    b <- coef(lm(mean_pred_oof ~ mean_obs, data = roi_tbl))
    cat(sprintf("ROI (out-of-fold): slope %.2f, intercept %.2f, crossover %.1f gC m-2\n",
                b[[2]], b[[1]], b[[1]] / (1 - b[[2]])))
  }
} else {
  cat("\nROI check skipped (set CFG$roi_path to the ROI polygons to run it).\n")
}

cat("\nP4_04_metrics.R COMPLETE\n")
