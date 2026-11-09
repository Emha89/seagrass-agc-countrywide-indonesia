# =============================================================================
# P4_02_tuning.R
# Paper 4 revised design: Random Forest hyperparameters for the GSE-only AGC
# model and the GSE-only cover model.
#
# Tuning uses pixel-grouped k-fold CV (quadrats of one 10-m pixel never split),
# so the grid search itself does not reward memorising shared pixels.
# The selected settings are then fixed in every validation design (P4_03)
# and copied to GEE (smileRandomForest).
#
# Reads : P4_samples.csv
# Writes: P4_tuning_AGC.csv, P4_tuning_COVER.csv, P4_best_params.csv
# =============================================================================

# setwd("path/to/seagrass-agc-countrywide-indonesia/R")   # folder that holds these scripts
source("P4_00_config_functions.R")

samples <- read_samples()

tune_grid <- function(df, response, label) {
  gt <- df |>
    dplyr::group_by(group = g_pixel) |>
    dplyr::summarise(n_agc = dplyr::n(), mean_agc = mean(.data[[response]]), n_cov = 0,
                     .groups = "drop")
  folds <- assign_folds(as.data.frame(gt), CFG$tune_folds, CFG$seed)
  f <- folds[df$g_pixel]

  grid <- expand.grid(mtry = CFG$grid_mtry, min_node_size = CFG$grid_node,
                      sample_fraction = CFG$grid_frac)
  grid <- grid[grid$mtry <= length(CFG$predictors), ]
  cat(sprintf("\n%s: %d samples, %d grid cells x %d folds\n",
              label, nrow(df), nrow(grid), CFG$tune_folds))

  res <- vector("list", nrow(grid))
  t0 <- Sys.time()
  for (i in seq_len(nrow(grid))) {
    p <- as.list(grid[i, ])
    pred <- rep(NA_real_, nrow(df))
    for (k in seq_len(CFG$tune_folds)) {
      tr <- f != k
      m <- fit_rf(df[tr, ], response, p, CFG$seed + k)
      pred[!tr] <- predict_rf(m, df[!tr, ])
    }
    res[[i]] <- dplyr::bind_cols(as.data.frame(p), reg_metrics(df[[response]], pred))
    cat(sprintf("  [%2d/%d] mtry=%-2d node=%-2d frac=%.1f  RMSE=%.3f  MEC=%.3f  (%.1f min)\n",
                i, nrow(grid), p$mtry, p$min_node_size, p$sample_fraction,
                res[[i]]$RMSE, res[[i]]$MEC, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  dplyr::bind_rows(res) |> dplyr::mutate(model = label, .before = 1)
}

# =============================================================================
# STEP 1: AGC model
# =============================================================================
agc <- samples |> dplyr::filter(is_agc)
tune_agc <- tune_grid(agc, CFG$agc_response, "AGC")
readr::write_csv(tune_agc, out_file("P4_tuning_AGC.csv"))

# =============================================================================
# STEP 2: Cover model
# =============================================================================
cov <- samples |> dplyr::filter(is_cov)
tune_cov <- tune_grid(cov, CFG$cover_response, "COVER")
readr::write_csv(tune_cov, out_file("P4_tuning_COVER.csv"))

# =============================================================================
# STEP 3: Best settings (lowest RMSE) and the GEE equivalents
# =============================================================================
best <- dplyr::bind_rows(tune_agc, tune_cov) |>
  dplyr::group_by(model) |>
  dplyr::slice_min(RMSE, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::mutate(num_trees = CFG$num_trees,
                gee_numberOfTrees = num_trees,
                gee_variablesPerSplit = mtry,
                gee_minLeafPopulation = min_node_size,
                gee_bagFraction = sample_fraction,
                cover_train_set = CFG$cover_train_set)
print(as.data.frame(best), row.names = FALSE)
readr::write_csv(best, out_file("P4_best_params.csv"))
cat("Saved: P4_best_params.csv\n")
cat("Note: ranger min.node.size and GEE minLeafPopulation are close but not identical rules.\n")

# =============================================================================
# STEP 4: Grid heatmaps (printed only)
# =============================================================================
for (tb in list(tune_agc, tune_cov)) {
  print(
    ggplot(tb, aes(x = factor(mtry), y = factor(min_node_size), fill = RMSE)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = sprintf("%.2f", RMSE)), size = 3, colour = "white") +
      facet_wrap(~ sample_fraction, labeller = label_both) +
      scale_fill_gradient(low = "steelblue", high = "firebrick") +
      labs(title = sprintf("%s grid search, pixel-grouped %d-fold CV", tb$model[1], CFG$tune_folds),
           x = "mtry", y = "min.node.size") +
      theme_minimal(base_size = 11)
  )
}

cat("\nP4_02_tuning.R COMPLETE\n")
