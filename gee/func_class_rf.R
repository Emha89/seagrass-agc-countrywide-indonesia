# =============================================================================
# func_class_rf.R
# RF utility functions -- classification stages (PA, leaf morphology).
# Grid tuning binary/multiclass (ranger), Monte Carlo robustness
# (randomForest), threshold sweep, and CI summary.
#
# NOTE: build_rf_formula() is also defined, identically, in func_reg_rf.R.
# If both files are sourced into the same session, the one sourced last
# simply overwrites the other -- harmless since the implementations are
# identical, just a minor redundancy worth knowing about.
#
# THRESHOLD DEFAULT NOTE: this file's functions don't share one consistent
# default classification threshold -- evaluate_rf_grid_ranger_classification()
# defaults to 0.7, while run_montecarlo_rf_classification() defaults to 0.5.
# evaluate_thresholds_classification() takes no single default since it
# sweeps a list of candidate thresholds. Same pattern noted in the Study 2
# repository's equivalent file -- worth confirming which threshold each
# calling script actually passes in, since the two defaults differ.
# =============================================================================

library(dplyr)
library(randomForest)
library(caret)
library(tibble)
library(purrr)

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
# evaluate_rf_grid_ranger_classification
# Grid search over ntree / mtry / nodesize / sample_fraction for the binary
# PA classifier, using ranger. Splits are stratified via
# caret::createDataPartition() (keeps the PA=1/PA=0 ratio consistent
# between train and test, unlike the regression grid search's plain random
# sampling -- deliberate, since regression has no class balance to
# preserve). Selects the best combination by highest mean accuracy.
# Returns list(grid_results, best_params).
# -----------------------------------------------------------------------------
evaluate_rf_grid_ranger_classification <- function(
    data, covars, response,
    ntree_values, mtry_values, nodesize_values, sample_fraction_values,
    threshold = 0.7, n_iter = 5, train_frac = 0.7, seed = 42
) {
  if (!requireNamespace("ranger", quietly = TRUE)) stop("Package 'ranger' required.")
  set.seed(seed)
  grid <- expand.grid(
    ntree = ntree_values, mtry = mtry_values,
    nodesize = nodesize_values, sample_fraction = sample_fraction_values,
    stringsAsFactors = FALSE
  )
  grid$accuracy_mean <- NA_real_; grid$accuracy_sd <- NA_real_
  grid$ci_lower <- NA_real_;      grid$ci_upper    <- NA_real_

  for (i in seq_len(nrow(grid))) {
    acc_iter <- numeric(n_iter)
    for (iter in seq_len(n_iter)) {
      idx        <- caret::createDataPartition(data[[response]], p = train_frac, list = FALSE)
      train_data <- data[idx, ]; test_data <- data[-idx, ]
      model <- ranger::ranger(
        formula         = build_rf_formula(covars, response),
        data            = train_data,
        num.trees       = grid$ntree[i],  mtry = grid$mtry[i],
        min.node.size   = grid$nodesize[i],
        sample.fraction = grid$sample_fraction[i],
        probability = TRUE, classification = TRUE, seed = seed + iter
      )
      probs  <- predict(model, data = test_data)$predictions[, "1"]
      preds  <- ifelse(probs > threshold, 1, 0)
      actual <- as.numeric(as.character(test_data[[response]]))
      acc_iter[iter] <- mean(preds == actual, na.rm = TRUE)
    }
    grid$accuracy_mean[i] <- mean(acc_iter); grid$accuracy_sd[i] <- sd(acc_iter)
    grid$ci_lower[i] <- quantile(acc_iter, 0.025); grid$ci_upper[i] <- quantile(acc_iter, 0.975)
  }
  best_idx <- which.max(grid$accuracy_mean)
  list(grid_results = grid, best_params = grid[best_idx, ])
}

# -----------------------------------------------------------------------------
# evaluate_rf_grid_ranger_multiclass
# Same grid-search pattern as the binary version above, for the 3-class
# leaf morphology (morph3) response. Predicted class is taken as the
# column with the highest predicted probability (max.col).
# Returns list(grid_results, best_params).
# -----------------------------------------------------------------------------
evaluate_rf_grid_ranger_multiclass <- function(
    data, covars, response,
    ntree_values, mtry_values, nodesize_values, sample_fraction_values,
    n_iter = 5, train_frac = 0.7, seed = 42
) {
  if (!requireNamespace("ranger", quietly = TRUE)) stop("Package 'ranger' required.")
  set.seed(seed)
  grid <- expand.grid(
    ntree = ntree_values, mtry = mtry_values,
    nodesize = nodesize_values, sample_fraction = sample_fraction_values,
    stringsAsFactors = FALSE
  )
  grid$accuracy_mean <- NA_real_; grid$accuracy_sd <- NA_real_
  grid$ci_lower <- NA_real_;      grid$ci_upper    <- NA_real_

  for (i in seq_len(nrow(grid))) {
    acc_iter <- numeric(n_iter)
    for (iter in seq_len(n_iter)) {
      idx        <- caret::createDataPartition(data[[response]], p = train_frac, list = FALSE)
      train_data <- data[idx, ]; test_data <- data[-idx, ]
      model <- ranger::ranger(
        formula         = build_rf_formula(covars, response),
        data            = train_data,
        num.trees       = grid$ntree[i],  mtry = grid$mtry[i],
        min.node.size   = grid$nodesize[i],
        sample.fraction = grid$sample_fraction[i],
        probability = TRUE, classification = TRUE, seed = seed + iter
      )
      probs       <- predict(model, data = test_data)$predictions
      pred_class  <- colnames(probs)[max.col(probs)]
      acc_iter[iter] <- mean(pred_class == as.character(test_data[[response]]), na.rm = TRUE)
    }
    grid$accuracy_mean[i] <- mean(acc_iter); grid$accuracy_sd[i] <- sd(acc_iter)
    grid$ci_lower[i] <- quantile(acc_iter, 0.025); grid$ci_upper[i] <- quantile(acc_iter, 0.975)
  }
  best_idx <- which.max(grid$accuracy_mean)
  list(grid_results = grid, best_params = grid[best_idx, ])
}

# -----------------------------------------------------------------------------
# evaluate_thresholds_classification
# Sweeps a list of candidate PA probability thresholds, repeating n_iter
# random 70/30 splits at each threshold, and reports mean/SD accuracy per
# threshold -- used to choose an operating threshold for the PA model.
# Uses base randomForest (not ranger), unlike the grid-search functions
# above.
# -----------------------------------------------------------------------------
evaluate_thresholds_classification <- function(
    data, covars, response,
    thresholds, n_iter = 10, ntree = 500, mtry = NULL
) {
  purrr::map_dfr(thresholds, function(thresh) {
    acc_vec <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      set.seed(i + 100)
      idx        <- caret::createDataPartition(data[[response]], p = 0.7, list = FALSE)
      train_data <- data[idx, ]; test_data <- data[-idx, ]
      model <- randomForest(
        formula   = build_rf_formula(covars, response),
        data      = train_data, ntree = ntree, mtry = mtry, importance = FALSE
      )
      probs      <- predict(model, newdata = test_data, type = "prob")[, "1"]
      preds      <- ifelse(probs > thresh, 1, 0)
      actuals    <- as.numeric(as.character(test_data[[response]]))
      acc_vec[i] <- mean(preds == actuals, na.rm = TRUE)
    }
    tibble::tibble(threshold = thresh, accuracy = mean(acc_vec), accuracy_sd = sd(acc_vec))
  })
}

# -----------------------------------------------------------------------------
# run_montecarlo_rf_classification
# n_iter random 70/30 (stratified) splits for the binary PA model ->
# accuracy, precision, recall, F1 per iteration, via caret's
# confusionMatrix(). Uses randomForest (the final-evaluation engine
# throughout this pipeline, distinct from ranger used for grid tuning).
# -----------------------------------------------------------------------------
run_montecarlo_rf_classification <- function(
    data, covars, response,
    n_iter = 100, train_frac = 0.7,
    ntree = 500, mtry = NULL, threshold = 0.5, seed = 42
) {
  set.seed(seed)
  seeds <- sample(1:10000, n_iter)
  purrr::map_dfr(seeds, function(s) {
    set.seed(s)
    idx        <- caret::createDataPartition(data[[response]], p = train_frac, list = FALSE)
    train_data <- data[idx, ]
    test_data  <- data[-idx, ] %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(covars), ~ !is.na(.)))
    model  <- randomForest(
      x = train_data[, covars], y = as.factor(train_data[[response]]),
      ntree = ntree, mtry = mtry, replace = FALSE
    )
    actual <- factor(test_data[[response]], levels = c(0, 1))
    probs  <- predict(model, newdata = test_data[, covars], type = "prob")[, "1"]
    preds  <- factor(ifelse(probs >= threshold, 1, 0), levels = c(0, 1))
    cm     <- caret::confusionMatrix(preds, actual, positive = "1")
    tibble::tibble(
      iteration = s,
      accuracy  = as.numeric(cm$overall["Accuracy"]),
      precision = as.numeric(cm$byClass["Precision"]),
      recall    = as.numeric(cm$byClass["Recall"]),
      f1        = as.numeric(cm$byClass["F1"])
    )
  })
}

# -----------------------------------------------------------------------------
# run_montecarlo_rf_multiclass
# Same pattern as run_montecarlo_rf_classification(), for the 3-class leaf
# morphology (morph3) response -- reports overall accuracy and kappa per
# iteration.
# -----------------------------------------------------------------------------
run_montecarlo_rf_multiclass <- function(
    data, covars, response,
    n_iter = 100, train_frac = 0.7,
    ntree = 500, mtry = NULL, seed = 42
) {
  set.seed(seed)
  seeds <- sample(1:10000, n_iter)
  purrr::map_dfr(seeds, function(s) {
    set.seed(s)
    idx        <- caret::createDataPartition(data[[response]], p = train_frac, list = FALSE)
    train_data <- data[idx, ]
    test_data  <- data[-idx, ] %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(covars), ~ !is.na(.)))
    model  <- randomForest(
      x = train_data[, covars], y = as.factor(train_data[[response]]),
      ntree = ntree, mtry = mtry, replace = FALSE
    )
    actual <- factor(test_data[[response]])
    preds  <- factor(predict(model, newdata = test_data[, covars]), levels = levels(actual))
    cm     <- caret::confusionMatrix(preds, actual)
    tibble::tibble(
      iteration = s,
      accuracy  = as.numeric(cm$overall["Accuracy"]),
      kappa     = as.numeric(cm$overall["Kappa"])
    )
  })
}

# -----------------------------------------------------------------------------
# summarize_ci_classification
# Mean, SD, and 95% CI for accuracy / precision / recall / F1 across Monte
# Carlo iterations (binary PA model).
# -----------------------------------------------------------------------------
summarize_ci_classification <- function(df) {
  df %>% dplyr::summarise(
    accuracy_mean   = mean(accuracy,  na.rm = TRUE),
    accuracy_sd     = sd(accuracy,    na.rm = TRUE),
    accuracy_lower  = quantile(accuracy,  0.025, na.rm = TRUE),
    accuracy_upper  = quantile(accuracy,  0.975, na.rm = TRUE),
    precision_mean  = mean(precision, na.rm = TRUE),
    precision_sd    = sd(precision,   na.rm = TRUE),
    precision_lower = quantile(precision, 0.025, na.rm = TRUE),
    precision_upper = quantile(precision, 0.975, na.rm = TRUE),
    recall_mean     = mean(recall,    na.rm = TRUE),
    recall_sd       = sd(recall,      na.rm = TRUE),
    recall_lower    = quantile(recall,    0.025, na.rm = TRUE),
    recall_upper    = quantile(recall,    0.975, na.rm = TRUE),
    f1_mean         = mean(f1,        na.rm = TRUE),
    f1_sd           = sd(f1,          na.rm = TRUE),
    f1_lower        = quantile(f1,        0.025, na.rm = TRUE),
    f1_upper        = quantile(f1,        0.975, na.rm = TRUE)
  )
}
