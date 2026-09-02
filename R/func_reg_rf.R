# =============================================================================
# func_reg_rf.R
# RF utility functions -- regression stages (SPC, AGB, carbon index, AGC).
# Grid tuning (ranger), Monte Carlo robustness (randomForest), apply-model
# with predictor-type coercion, and CI summary.
#
# NOTE: build_rf_formula() is also defined, identically, in func_class_rf.R.
# If both files are sourced into the same session, the one sourced last
# simply overwrites the other -- harmless since the implementations are
# identical, just a minor redundancy worth knowing about.
# =============================================================================

library(dplyr)
library(randomForest)
library(caret)
library(tibble)
library(purrr)

`%||%` <- function(x, y) if (!is.null(x)) x else y

# -----------------------------------------------------------------------------
# build_rf_formula
# Builds a formula object (response ~ covar1 + covar2 + ...) from a
# character vector of covariate names and a response name.
# -----------------------------------------------------------------------------
build_rf_formula <- function(covars, response) {
  if (length(covars) == 0 || is.null(covars))
    stop("Covariate list is empty or NULL.")
  as.formula(paste(response, "~", paste(covars, collapse = " + ")))
}

# -----------------------------------------------------------------------------
# evaluate_rf_grid_ranger_regression
# Grid search over ntree / mtry / nodesize / sample_fraction using ranger,
# with n_iter random 70/30 (or train_frac) splits per grid cell. Selects the
# best combination by lowest median RMSE across iterations (median, not
# mean, for robustness to any single unusually bad split).
# Returns list(grid_results, best_params).
# -----------------------------------------------------------------------------
evaluate_rf_grid_ranger_regression <- function(
    data, covars, response,
    ntree_values           = c(300, 500, 700),
    mtry_values            = NULL,
    nodesize_values        = c(3, 5, 7),
    sample_fraction_values = c(0.6, 0.8, 1.0),
    n_iter = 5, train_frac = 0.7, seed = 42
) {
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("Package 'ranger' is required.")
  if (is.null(mtry_values)) {
    p           <- length(covars)
    mtry_values <- unique(c(floor(sqrt(p)), floor(p / 3)))
  }

  results <- expand.grid(
    ntree = ntree_values, mtry = mtry_values,
    nodesize = nodesize_values, sample_fraction = sample_fraction_values,
    stringsAsFactors = FALSE
  )
  results$rmse_median <- NA_real_
  results$rmse_mean   <- NA_real_
  results$rmse_sd     <- NA_real_
  results$ci_lower    <- NA_real_
  results$ci_upper    <- NA_real_

  for (i in seq_len(nrow(results))) {
    row      <- results[i, ]
    rmse_vec <- numeric(n_iter)
    for (j in seq_len(n_iter)) {
      set.seed(seed + j)
      idx        <- sample(seq_len(nrow(data)), floor(train_frac * nrow(data)))
      train_data <- data[idx, ]
      test_data  <- data[-idx, ]
      model <- ranger::ranger(
        formula         = build_rf_formula(covars, response),
        data            = train_data,
        num.trees       = row$ntree,
        mtry            = row$mtry,
        min.node.size   = row$nodesize,
        sample.fraction = row$sample_fraction,
        replace         = FALSE,
        importance      = "impurity",
        seed            = seed + j
      )
      preds        <- predict(model, data = test_data)$predictions
      rmse_vec[j]  <- sqrt(mean((preds - test_data[[response]])^2, na.rm = TRUE))
    }
    results$rmse_median[i] <- median(rmse_vec)
    results$rmse_mean[i]   <- mean(rmse_vec)
    results$rmse_sd[i]     <- sd(rmse_vec)
    results$ci_lower[i]    <- quantile(rmse_vec, 0.025)
    results$ci_upper[i]    <- quantile(rmse_vec, 0.975)
  }
  best <- results[which.min(results$rmse_median), ]
  list(grid_results = results, best_params = best)
}

# -----------------------------------------------------------------------------
# run_montecarlo_rf_regression
# n_iter random train/test splits (train_frac each) -> RMSE, MAE, R2 per
# iteration, using randomForest (the final-evaluation engine throughout
# this pipeline, distinct from ranger used for grid tuning above).
# -----------------------------------------------------------------------------
run_montecarlo_rf_regression <- function(
    data, covars, response,
    n_iter = 100, ntree = 500, mtry = NULL,
    nodesize = 5, sample_fraction = 0.8,
    train_frac = 0.7, seed = 42
) {
  set.seed(seed)
  seeds <- sample(1:10000, n_iter)

  purrr::map_dfr(seeds, function(s) {
    set.seed(s)
    idx   <- sample(seq_len(nrow(data)), floor(train_frac * nrow(data)))
    train <- data[idx, ]
    test  <- data[-idx, ] %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(covars), ~ !is.na(.)))

    model <- randomForest(
      x        = train[, covars],
      y        = train[[response]],
      ntree    = ntree,
      mtry     = mtry %||% floor(sqrt(length(covars))),
      nodesize = nodesize,
      replace  = FALSE,
      sampsize = floor(sample_fraction * nrow(train))
    )
    preds  <- predict(model, newdata = test[, covars])
    actual <- test[[response]]
    tibble::tibble(
      iteration = s,
      RMSE = caret::RMSE(preds, actual),
      MAE  = caret::MAE(preds, actual),
      R2   = caret::R2(preds, actual)
    )
  })
}

# -----------------------------------------------------------------------------
# summarize_ci_regression
# Mean, SD, and 95% CI for RMSE / MAE / R2 across Monte Carlo iterations.
# -----------------------------------------------------------------------------
summarize_ci_regression <- function(df) {
  df %>% dplyr::summarise(
    RMSE_mean  = mean(RMSE,  na.rm = TRUE),
    RMSE_sd    = sd(RMSE,    na.rm = TRUE),
    RMSE_lower = quantile(RMSE, 0.025, na.rm = TRUE),
    RMSE_upper = quantile(RMSE, 0.975, na.rm = TRUE),
    MAE_mean   = mean(MAE,   na.rm = TRUE),
    MAE_sd     = sd(MAE,     na.rm = TRUE),
    MAE_lower  = quantile(MAE,  0.025, na.rm = TRUE),
    MAE_upper  = quantile(MAE,  0.975, na.rm = TRUE),
    R2_mean    = mean(R2,    na.rm = TRUE),
    R2_sd      = sd(R2,      na.rm = TRUE),
    R2_lower   = quantile(R2,   0.025, na.rm = TRUE),
    R2_upper   = quantile(R2,   0.975, na.rm = TRUE)
  )
}

# -----------------------------------------------------------------------------
# apply_model_regression
# Applies a saved RF model to new data, coercing each predictor column to
# the type recorded in predictor_structure (as saved alongside the model
# at training time) before prediction.
# -----------------------------------------------------------------------------
apply_model_regression <- function(df, predictors, model, predictor_structure) {
  df_pred <- df %>% dplyr::select(dplyr::all_of(predictors))
  for (var in predictors) {
    cls <- predictor_structure$classes[[var]][1]
    if      (cls == "numeric")   df_pred[[var]] <- as.numeric(df_pred[[var]])
    else if (cls == "integer")   df_pred[[var]] <- as.integer(df_pred[[var]])
    else if (cls == "factor") {
      lvls <- predictor_structure$levels[[var]]
      df_pred[[var]] <- factor(df_pred[[var]], levels = lvls)
    } else if (cls == "character") df_pred[[var]] <- as.character(df_pred[[var]])
  }
  predict(model, newdata = df_pred)
}