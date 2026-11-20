# =============================================================================
# P4_04_metrics.R
# Paper 4, design v3: accuracy at the sample level and at the inventory level.
#
# Two sample sets:
#   all  = every AGC sample; methods Tier 1, Tier 2, Local mean, Pixel-based GSE
#   s2   = AGC samples with a clear Sentinel-2 pixel; methods Tier 1, Tier 2,
#          Local mean, Pixel-based GSE, Pixel-based S2 and (check) Pixel-based S2
#          (median), every one trained on the same subset (the *_sub columns and
#          pred_s2*), with distances to the nearest training sample of the subset.
#          Only when P4_03 ran the Sentinel-2 models. Compare methods within a set;
#          "all" vs "s2" for GSE shows the effect of the smaller sample set.
#
# Inventory level (STEP 6): out-of-fold predictions are averaged per unit, in the
# design that holds the whole unit out, and compared with the unit's field mean:
#   site      100 m block, main survey year, n >= CFG$min_n_site     (b100 design)
#   meadow    DBSCAN cluster, n >= CFG$min_n_meadow                  (cluster design)
#   ecoregion MEOW ecoregion, n >= CFG$min_n_eco                     (ecoregion design)
# Metrics on the unit means (RMSE, MAE, R2, MEC, bias, slope) per repeat, then
# the mean and the 2.5-97.5 percentiles across repeats.
#
# Reads : P4_oof_agc_<design>.csv, P4_samples.csv
# Writes: P4_metrics_by_repeat.csv     sample level, design x repeat x method x set
#         P4_metrics_summary.csv       mean and 2.5-97.5 percentiles + median distance
#         P4_metrics_ecoregion.csv     per-ecoregion metrics (LOLO), n flag
#         P4_metrics_ecoregion_equal_weight.csv
#         P4_lolo_region_means.csv     regional mean AGC per method when the ecoregion is held out
#         P4_repeat_stability.csv      cumulative mean MEC by repeat
#         P4_error_by_agc_bin.csv      absolute error per 10 gC m-2 bin (pixel-based GSE)
#         P4_inventory_units.csv       every unit: field mean and mean prediction per method
#         P4_inventory_metrics_by_repeat.csv, P4_inventory_metrics.csv
#         P4_inventory_site_by_ecoregion.csv
#         P4_shrinkage_groups.csv, P4_shrinkage_fit.csv
#         P4_measured_metrics.csv      (only when P4_samples.csv has AGC_measured)
#         P4_fig_*.png
# =============================================================================

source("R/P4_00_config_functions.R")

oof <- read_oof()
samples <- read_samples()

design_levels <- CFG$designs
q025 <- function(x) quantile(x, 0.025, na.rm = TRUE)
q975 <- function(x) quantile(x, 0.975, na.rm = TRUE)

S2_ON <- "pred_s2" %in% names(oof) && any(is.finite(oof$pred_s2))
SETS <- list(all = c("Tier 1" = "pred_t1", "Tier 2" = "pred_t2", "Local mean" = "pred_local",
                     "Pixel-based GSE" = "pred_pixel"))
if (S2_ON) {
  SETS$s2 <- c("Tier 1" = "pred_t1_sub", "Tier 2" = "pred_t2_sub", "Local mean" = "pred_local_sub",
               "Pixel-based GSE" = "pred_gse_sub", "Pixel-based S2" = "pred_s2")
  if ("pred_s2med" %in% names(oof) && any(is.finite(oof$pred_s2med)))
    SETS$s2 <- c(SETS$s2, "Pixel-based S2 (median)" = "pred_s2med")
}
method_levels <- METHOD_LEVELS
cat(sprintf("Sample sets: %s\n", paste(names(SETS), collapse = ", ")))

# Long table: one row per test sample x method x set
to_long <- function(d, set) {
  cols <- SETS[[set]]
  if (set == "s2") d <- d |> dplyr::filter(has_s2, is.finite(pred_s2)) |> dplyr::mutate(dist_m = dist_s2_m)
  d |>
    tidyr::pivot_longer(dplyr::all_of(unname(cols)), names_to = "col", values_to = "pred") |>
    dplyr::mutate(sample_set = set,
                  method = factor(names(cols)[match(col, cols)], levels = method_levels),
                  design = factor(design, levels = design_levels))
}
long <- dplyr::bind_rows(lapply(names(SETS), function(s) to_long(oof, s)))

# =============================================================================
# STEP 1: Sample-level metrics per design x repeat x method (pooled out-of-fold)
# =============================================================================
m_rep <- long |>
  dplyr::group_by(sample_set, design, rep, method) |>
  dplyr::summarise(reg_metrics(obs, pred), .groups = "drop")

# median distance of the test samples of each set (Tier 1 rows: one per test sample)
dist_rep <- long |>
  dplyr::filter(method == "Tier 1") |>
  dplyr::group_by(sample_set, design, rep) |>
  dplyr::summarise(median_dist_m = median(dist_m), .groups = "drop")

m_rep <- dplyr::left_join(m_rep, dist_rep, by = c("sample_set", "design", "rep"))
readr::write_csv(m_rep, out_file("P4_metrics_by_repeat.csv"))

m_sum <- m_rep |>
  dplyr::group_by(sample_set, design, method) |>
  dplyr::summarise(
    n_repeats = dplyr::n(), n = mean(n),
    median_dist_m = mean(median_dist_m),
    dplyr::across(c(RMSE, MAE, R2, MEC, bias_pct, slope),
                  list(mean = ~ mean(.x, na.rm = TRUE), q025 = q025, q975 = q975)),
    .groups = "drop") |>
  dplyr::arrange(sample_set, design, method)
print_tbl(m_sum |> dplyr::select(sample_set, design, method, n, median_dist_m, RMSE_mean, MAE_mean,
                                 R2_mean, MEC_mean, MEC_q025, MEC_q975, bias_pct_mean))
readr::write_csv(m_sum, out_file("P4_metrics_summary.csv"))

t2_fb <- oof |>
  dplyr::group_by(design) |>
  dplyr::summarise(pct_tier2_fallback = 100 * mean(t2_fallback),
                   pct_local_fallback = 100 * mean(local_fallback), .groups = "drop")
cat("\nTest samples whose ecoregion had < CFG$min_n_eco training samples (Tier 2 = Tier 1),\n",
    "and whose meadow had < CFG$min_n_local training samples (Local mean = Tier 2):\n")
print(as.data.frame(t2_fb), row.names = FALSE)

if (S2_ON) {
  s2_share <- oof |>
    dplyr::group_by(design) |>
    dplyr::summarise(pct_test_with_s2 = 100 * mean(has_s2 & is.finite(pred_s2)), .groups = "drop")
  cat("\nTest samples in the Sentinel-2 set:\n")
  print(as.data.frame(s2_share), row.names = FALSE)
}

# =============================================================================
# STEP 2: Leave-one-ecoregion-out, per ecoregion (all samples)
# =============================================================================
eco_long <- long |> dplyr::filter(design == "ecoregion", sample_set == "all")
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

  # Regional mean AGC: what matters for an area total (Tier 2 = Tier 1 in this design)
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
  cat(sprintf("Pixel-based GSE closer to the field mean than Tier 1 in %d of %d ecoregions\n",
              sum(rr$pixel_closer_than_tier1), nrow(rr)))
}

# =============================================================================
# STEP 3: Stability of repeats (cumulative mean MEC)
# =============================================================================
stab <- m_rep |>
  dplyr::filter(design != "ecoregion") |>
  dplyr::arrange(sample_set, design, method, rep) |>
  dplyr::group_by(sample_set, design, method) |>
  dplyr::mutate(cum_mean_MEC = cumsum(MEC) / dplyr::row_number()) |>
  dplyr::ungroup()
readr::write_csv(stab, out_file("P4_repeat_stability.csv"))

if (max(stab$rep) > 1) {
  stab_last <- stab |>
    dplyr::group_by(sample_set, design, method) |>
    dplyr::summarise(last_change = abs(dplyr::last(cum_mean_MEC) - dplyr::nth(cum_mean_MEC, -2)),
                     stable = last_change < 0.01, .groups = "drop")
  cat("\nChange in cumulative mean MEC at the last repeat (stable if < 0.01):\n")
  print_tbl(stab_last, 4)
  if (any(!stab_last$stable))
    cat("Some designs are not yet stable: raise CFG$n_repeats and rerun P4_03 for those designs.\n")
}

p_stab <- ggplot(stab, aes(x = rep, y = cum_mean_MEC, colour = method, linetype = sample_set)) +
  geom_line() + geom_point(size = 1) +
  facet_wrap(~ design, nrow = 1) +
  labs(title = "Cumulative mean MEC by number of repeats", x = "Repeat", y = "Cumulative mean MEC",
       colour = NULL, linetype = "Sample set") +
  theme_minimal(base_size = 11)
print(p_stab)
ggsave(out_file("P4_fig_repeat_stability.png"), p_stab, width = 10, height = 4, dpi = 300)

# =============================================================================
# STEP 4: MEC against median distance to the nearest training sample
# =============================================================================
p_dist <- ggplot(m_sum, aes(x = median_dist_m, y = MEC_mean, colour = method, group = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line() +
  geom_errorbar(aes(ymin = MEC_q025, ymax = MEC_q975), width = 0.04) +
  geom_point(size = 2) +
  geom_text(data = dplyr::filter(m_sum, method == "Pixel-based GSE"),
            aes(label = design), vjust = -1.1, size = 3, colour = "grey20", show.legend = FALSE) +
  facet_wrap(~ sample_set, labeller = labeller(sample_set = c(all = "All AGC samples",
                                                               s2 = "Samples with a clear Sentinel-2 pixel"))) +
  scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  labs(title = "Sample-level accuracy against distance from field data",
       x = "Median distance from test to nearest training sample (m, log scale)",
       y = "MEC (mean, 2.5-97.5 percentiles across repeats)", colour = NULL) +
  theme_minimal(base_size = 11)
print(p_dist)
ggsave(out_file("P4_fig_mec_distance.png"), p_dist, width = 10, height = 5, dpi = 300)

# =============================================================================
# STEP 5: Pixel-based GSE absolute error per 10 gC m-2 bin (out-of-fold)
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

p_bin <- ggplot(err_bin, aes(x = bin, y = median_abs_err)) +
  geom_pointrange(aes(ymin = q25, ymax = q75)) +
  facet_wrap(~ design, ncol = 1) +
  labs(title = "Pixel-based GSE absolute error by observed AGC (out-of-fold)",
       x = "Observed AGC bin (gC m-2)", y = "Absolute error (gC m-2), median and IQR") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_bin)
ggsave(out_file("P4_fig_error_by_agc_bin.png"), p_bin, width = 8, height = 6, dpi = 300)

# =============================================================================
# STEP 6: Inventory-level validation (unit means of out-of-fold predictions)
# =============================================================================
agc <- samples |> dplyr::filter(is_agc)
UNITS <- list(
  site      = list(design = "b100",      label = "Site (100 m block)",     min_n = CFG$min_n_site),
  meadow    = list(design = "cluster",   label = "Meadow (DBSCAN cluster)", min_n = CFG$min_n_meadow),
  ecoregion = list(design = "ecoregion", label = "Ecoregion",               min_n = CFG$min_n_eco)
)

# Members of every unit (gee_id, year, unit); the threshold is applied per sample set
# site: the survey year with most AGC samples in the block (earliest year on ties)
site_year <- agc |>
  dplyr::count(g_b100, year, name = "n_year") |>
  dplyr::arrange(g_b100, dplyr::desc(n_year), year) |>
  dplyr::distinct(g_b100, .keep_all = TRUE)
members <- dplyr::bind_rows(
  agc |> dplyr::semi_join(site_year, by = c("g_b100", "year")) |>
    dplyr::transmute(level = "site", unit = paste(g_b100, year, sep = "|"), gee_id, year),
  agc |> dplyr::filter(grepl("^c", g_cluster)) |>
    dplyr::transmute(level = "meadow", unit = g_cluster, gee_id, year),
  agc |> dplyr::transmute(level = "ecoregion", unit = ecoregion, gee_id, year)
)

unit_means <- function(lv, set) {
  U <- UNITS[[lv]]
  cols <- SETS[[set]]
  d <- oof |>
    dplyr::filter(design == U$design) |>
    dplyr::inner_join(members |> dplyr::filter(level == lv) |> dplyr::select(gee_id, year, unit),
                      by = c("gee_id", "year"))
  if (set == "s2") d <- d |> dplyr::filter(has_s2, is.finite(pred_s2)) |> dplyr::mutate(dist_m = dist_s2_m)
  d |>
    dplyr::group_by(rep, unit) |>
    # n_years before year: summarise() evaluates in order
    dplyr::summarise(n = dplyr::n(), ecoregion = dplyr::first(ecoregion),
                     n_years = dplyr::n_distinct(year), year = min(year), obs = mean(obs),
                     dplyr::across(dplyr::all_of(unname(cols)), mean),
                     median_dist_m = median(dist_m), .groups = "drop") |>
    dplyr::filter(n >= U$min_n) |>
    tidyr::pivot_longer(dplyr::all_of(unname(cols)), names_to = "col", values_to = "pred") |>
    dplyr::mutate(level = lv, sample_set = set,
                  method = factor(names(cols)[match(col, cols)], levels = method_levels)) |>
    dplyr::select(-col)
}
inv <- purrr::map_dfr(names(SETS), function(s) purrr::map_dfr(names(UNITS), unit_means, set = s)) |>
  dplyr::mutate(level = factor(level, levels = names(UNITS)))

inv_units <- inv |>
  dplyr::group_by(sample_set, level, unit, ecoregion, year, n_years, method) |>
  dplyr::summarise(n = mean(n), obs = mean(obs), pred = mean(pred),
                   median_dist_m = mean(median_dist_m), .groups = "drop") |>
  tidyr::pivot_wider(names_from = method, values_from = pred) |>
  dplyr::arrange(sample_set, level, obs)
readr::write_csv(inv_units, out_file("P4_inventory_units.csv"))

inv_rep <- inv |>
  dplyr::group_by(sample_set, level, method, rep) |>
  dplyr::summarise(n_units = dplyr::n(), n_samples = sum(n),
                   reg_metrics(obs, pred), .groups = "drop")
readr::write_csv(inv_rep, out_file("P4_inventory_metrics_by_repeat.csv"))

inv_sum <- inv_rep |>
  dplyr::group_by(sample_set, level, method) |>
  dplyr::summarise(n_repeats = dplyr::n(), n_units = mean(n_units), n_samples = mean(n_samples),
                   dplyr::across(c(RMSE, MAE, R2, MEC, bias_pct, slope),
                                 list(mean = ~ mean(.x, na.rm = TRUE), q025 = q025, q975 = q975)),
                   .groups = "drop") |>
  dplyr::arrange(sample_set, level, method)
cat("\nInventory-level accuracy (metrics on unit means; mean over repeats):\n")
print_tbl(inv_sum |> dplyr::select(sample_set, level, method, n_units, RMSE_mean, MAE_mean, R2_mean,
                                   MEC_mean, MEC_q025, MEC_q975, bias_pct_mean, slope_mean), 3)
readr::write_csv(inv_sum, out_file("P4_inventory_metrics.csv"))

# Sites per ecoregion (sites are mostly in Lesser Sunda)
site_eco <- inv |>
  dplyr::filter(level == "site") |>
  dplyr::group_by(sample_set, ecoregion, method, rep) |>
  dplyr::summarise(n_units = dplyr::n(), reg_metrics(obs, pred), .groups = "drop") |>
  dplyr::group_by(sample_set, ecoregion, method) |>
  dplyr::summarise(n_sites = mean(n_units),
                   dplyr::across(c(RMSE, MAE, R2, MEC, bias_pct, slope), ~ mean(.x, na.rm = TRUE)),
                   .groups = "drop") |>
  dplyr::mutate(reported = n_sites >= CFG$min_sites_eco) |>
  dplyr::arrange(sample_set, dplyr::desc(n_sites), method)
cat("\nSite-level accuracy per ecoregion:\n")
print_tbl(site_eco |> dplyr::filter(reported), 3)
readr::write_csv(site_eco, out_file("P4_inventory_site_by_ecoregion.csv"))

nat_mean <- mean(agc$AGC_pred)
p_inv <- ggplot(inv_units |>
                  tidyr::pivot_longer(dplyr::any_of(method_levels), names_to = "method", values_to = "pred") |>
                  dplyr::filter(!is.na(pred)) |>
                  dplyr::mutate(method = factor(method, levels = method_levels)),
                aes(x = obs, y = pred)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = nat_mean, linetype = "dotted", colour = "firebrick") +
  geom_point(aes(size = n), alpha = 0.5, colour = "steelblue4") +
  facet_grid(sample_set + level ~ method) +
  scale_size_area(max_size = 3) +
  labs(title = "Unit mean AGC: field vs out-of-fold prediction (unit held out)",
       subtitle = "Dashed = 1:1, dotted = national mean of the AGC samples",
       x = "Field mean (gC m-2)", y = "Mean prediction (gC m-2)", size = "n") +
  theme_bw(base_size = 9)
print(p_inv)
ggsave(out_file("P4_fig_inventory_units.png"), p_inv,
       width = 1 + 2.3 * length(unique(unlist(lapply(SETS, names)))), height = 2.2 * length(SETS) * 3,
       dpi = 300, limitsize = FALSE)

# =============================================================================
# STEP 7: Shrinkage toward the training mean (leak-free, pixel-based GSE)
# Group mean of observed AGC vs group mean of out-of-fold pixel-based predictions
# =============================================================================
per_sample <- oof |>
  dplyr::group_by(design, gee_id, g_cluster, ecoregion, obs) |>
  dplyr::summarise(pred = mean(pred_pixel), .groups = "drop")

shr_groups <- dplyr::bind_rows(
  per_sample |> dplyr::filter(design == "cluster", grepl("^c", g_cluster)) |>
    dplyr::group_by(level = "DBSCAN cluster (cluster design)", group = g_cluster) |>
    dplyr::summarise(n = dplyr::n(), ecoregion = dplyr::first(ecoregion),
                     mean_obs = mean(obs), mean_pred = mean(pred), .groups = "drop") |>
    dplyr::filter(n >= CFG$min_n_meadow),
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
  dplyr::mutate(national_training_mean = nat_mean)
cat("\nShrinkage fit (slope < 1 = predictions compressed toward the training mean):\n")
print_tbl(shr_fit)
readr::write_csv(shr_fit, out_file("P4_shrinkage_fit.csv"))

# =============================================================================
# STEP 8 (optional): out-of-fold predictions against measured AGC
# =============================================================================
if ("AGC_measured" %in% names(samples) && any(is.finite(samples$AGC_measured))) {
  meas <- samples |> dplyr::filter(is_agc, is.finite(AGC_measured)) |>
    dplyr::select(gee_id, year, AGC_measured)
  mm <- long |>
    dplyr::inner_join(meas, by = c("gee_id", "year")) |>
    dplyr::group_by(sample_set, design, rep, method) |>
    dplyr::summarise(reg_metrics(AGC_measured, pred), .groups = "drop") |>
    dplyr::group_by(sample_set, design, method) |>
    dplyr::summarise(dplyr::across(c(n, RMSE, MAE, R2, MEC, bias_pct, slope), ~ mean(.x, na.rm = TRUE)),
                     .groups = "drop")
  ref <- reg_metrics(meas$AGC_measured, samples$AGC_pred[match(paste(meas$gee_id, meas$year),
                                                               paste(samples$gee_id, samples$year))])
  cat(sprintf("\nMeasured AGC: %d samples | AGC_pred (cover-based reference) vs measured: RMSE %.2f, MEC %.3f\n",
              nrow(meas), ref$RMSE, ref$MEC))
  print_tbl(mm, 3)
  readr::write_csv(dplyr::bind_rows(mm, ref |> dplyr::mutate(sample_set = "reference",
                                                              design = "none", method = "AGC_pred")),
                   out_file("P4_measured_metrics.csv"))
}

cat("\nP4_04_metrics.R COMPLETE\n")
