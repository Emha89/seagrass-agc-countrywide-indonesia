# =============================================================================
# P4_02_tuning.R
# Paper 4, design v3: Random Forest hyperparameters.
#
#   AGC          GSE, all AGC samples (the model behind the national map; as P4v2)
#   AGC_GSE_SUB  GSE, AGC samples with a clear Sentinel-2 pixel (has_s2)
#   AGC_S2       Sentinel-2 percentiles, the same samples and folds as AGC_GSE_SUB
#   AGC_S2MED    Sentinel-2 median of each band (check; CFG$s2_median_check)
# The subset models share samples, folds, seeds and grid, so the GSE and
# Sentinel-2 benchmark differ only in the predictors. mtry values above the number
# of predictors are dropped; a small predictor set (the median model) gets
# floor(sqrt(p)) and floor(p / 3) added, so it is not limited to bagging. The grid
# of the GSE models (64 predictors) is that of P4v2.
#
# Tuning uses pixel-grouped k-fold CV (quadrats of one 10-m pixel never split),
# so the grid search itself does not reward memorising shared pixels.
# The selected settings are then fixed in every validation design (P4_03);
# the AGC settings are the GEE smileRandomForest settings of the map.
#
# Reads : P4_samples.csv
# Writes: P4_tuning_AGC.csv, P4_tuning_AGC_GSE_SUB.csv, P4_tuning_AGC_S2.csv,
#         P4_tuning_AGC_S2MED.csv, P4_best_params.csv, P4_fig_tuning_<model>.png
# =============================================================================

source("R/P4_00_config_functions.R")

samples <- read_samples()

tune_grid <- function(df, response, label, vars) {
  gt <- df |>
    dplyr::group_by(group = g_pixel) |>
    dplyr::summarise(n_agc = dplyr::n(), mean_agc = mean(.data[[response]]), n_cov = 0,
                     .groups = "drop")
  folds <- assign_folds(as.data.frame(gt), CFG$tune_folds, CFG$seed)
  f <- folds[df$g_pixel]

  p_n <- length(vars)
  mtry_set <- CFG$grid_mtry[CFG$grid_mtry <= p_n]
  if (p_n < max(CFG$grid_mtry))
    mtry_set <- sort(unique(c(mtry_set, max(1, floor(sqrt(p_n))), max(1, floor(p_n / 3)))))
  grid <- expand.grid(mtry = mtry_set, min_node_size = CFG$grid_node,
                      sample_fraction = CFG$grid_frac)
  cat(sprintf("\n%s: %d samples, %d predictors, %d grid cells x %d folds\n",
              label, nrow(df), length(vars), nrow(grid), CFG$tune_folds))

  res <- vector("list", nrow(grid))
  t0 <- Sys.time()
  for (i in seq_len(nrow(grid))) {
    p <- as.list(grid[i, ])
    pred <- rep(NA_real_, nrow(df))
    for (k in seq_len(CFG$tune_folds)) {
      tr <- f != k
      m <- fit_rf(df[tr, ], response, p, CFG$seed + k, vars)
      pred[!tr] <- predict_rf(m, df[!tr, ])
    }
    res[[i]] <- dplyr::bind_cols(as.data.frame(p), reg_metrics(df[[response]], pred))
    cat(sprintf("  [%2d/%d] mtry=%-2d node=%-2d frac=%.1f  RMSE=%.3f  MEC=%.3f  (%.1f min)\n",
                i, nrow(grid), p$mtry, p$min_node_size, p$sample_fraction,
                res[[i]]$RMSE, res[[i]]$MEC, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  dplyr::bind_rows(res) |> dplyr::mutate(model = label, n_predictors = length(vars), .before = 1)
}

# =============================================================================
# STEP 1: GSE model, all AGC samples
# =============================================================================
agc <- samples |> dplyr::filter(is_agc)
tunes <- list(AGC = tune_grid(agc, CFG$agc_response, "AGC", CFG$predictors))

# =============================================================================
# STEP 2: GSE and Sentinel-2 on the Sentinel-2 subset (same samples and folds)
# =============================================================================
S2_ON <- isTRUE(CFG$run_s2) && "has_s2" %in% names(samples) && sum(samples$has_s2) > 0
if (S2_ON) {
  sub <- samples |> dplyr::filter(has_s2)
  cat(sprintf("\nSentinel-2 subset: %d of %d AGC samples (%.1f%%)\n",
              nrow(sub), nrow(agc), 100 * nrow(sub) / nrow(agc)))
  tunes$AGC_GSE_SUB <- tune_grid(sub, CFG$agc_response, "AGC_GSE_SUB", CFG$predictors)
  tunes$AGC_S2      <- tune_grid(sub, CFG$agc_response, "AGC_S2", CFG$s2_vars)
  if (isTRUE(CFG$s2_median_check))
    tunes$AGC_S2MED <- tune_grid(sub, CFG$agc_response, "AGC_S2MED", CFG$s2_med_vars)
} else {
  cat("\nSentinel-2 models skipped (no has_s2 samples, or CFG$run_s2 = FALSE).\n")
}
for (nm in names(tunes)) readr::write_csv(tunes[[nm]], out_file(sprintf("P4_tuning_%s.csv", nm)))

# =============================================================================
# STEP 3: Best settings (lowest RMSE) and the GEE equivalents
# =============================================================================
best <- dplyr::bind_rows(tunes) |>
  dplyr::group_by(model) |>
  dplyr::slice_min(RMSE, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::mutate(num_trees = CFG$num_trees,
                gee_numberOfTrees = num_trees,
                gee_variablesPerSplit = mtry,
                gee_minLeafPopulation = min_node_size,
                gee_bagFraction = sample_fraction) |>
  dplyr::arrange(match(model, c("AGC", "AGC_GSE_SUB", "AGC_S2", "AGC_S2MED")))
print(as.data.frame(best), row.names = FALSE)
readr::write_csv(best, out_file("P4_best_params.csv"))
cat("Saved: P4_best_params.csv\n")
cat("Note: ranger min.node.size and GEE minLeafPopulation are close but not identical rules.\n")

# =============================================================================
# STEP 4: Grid heatmaps
# =============================================================================
for (tb in tunes) {
  p_grid <- ggplot(tb, aes(x = factor(mtry), y = factor(min_node_size), fill = RMSE)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", RMSE)), size = 3, colour = "white") +
    facet_wrap(~ sample_fraction, labeller = label_both) +
    scale_fill_gradient(low = "steelblue", high = "firebrick") +
    labs(title = sprintf("%s grid search, pixel-grouped %d-fold CV (n = %d)", tb$model[1],
                         CFG$tune_folds, tb$n[1]),
         x = "mtry", y = "min.node.size") +
    theme_minimal(base_size = 11)
  print(p_grid)
  ggsave(out_file(sprintf("P4_fig_tuning_%s.png", tb$model[1])), p_grid, width = 8, height = 4, dpi = 300)
}

cat("\nP4_02_tuning.R COMPLETE\n")
