# =============================================================================
# A1_rf_model_PA.R
# RF classification -- presence/absence (PA).
#
# Reads: master_raw_YYYY.csv  (from 01_build_master_raw.R)
# Saves:
#   final_rf_model_PA.rds
#   rf_predictor_structure_PA.rds
#   training_PA.csv              -- filtered training subset
#   var_importance_rf_PA.csv
#   threshold_eval_PA.csv
#   grid_rf_tuning_PA.csv / grid_rf_best_PA.csv
#   rf_pa_mc_results.csv / rf_pa_mc_summary_ci.csv
#
# PREDICTOR NOTE: env_vars is depth only here -- distToLand is present as a
# commented-out alternative (kept below) but not active. This differs from
# Study 2's PA model, which uses both depth and distToLand. Kept exactly
# as in the original script; worth confirming this is the intended
# national-scale predictor set before treating it as final.
#
# Threshold: auto-tuned threshold is computed and printed for comparison,
# but the actual model uses a manually fixed threshold of 0.6 (same
# pattern as Study 2's PA model). This value is explicitly passed to the
# Monte Carlo evaluation call, so it's unaffected by the differing
# default thresholds across functions in func_class_rf.R (0.5 vs 0.7).
#
# NOTE: normalise_master_types() is duplicated here from
# 01_build_master_raw.R (not sourced from a shared location) -- harmless
# redundancy, not a bug.
#
# Data availability: this script's inputs are not included in this
# repository.
# =============================================================================

library(dplyr); library(readr); library(randomForest)
library(caret);  library(ranger); library(purrr)
library(ggplot2); library(tibble); library(tidyr); library(maps)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "func_class_rf.R"))


# Helper: normalise column types -- called after column selection
# so only known columns remain; no ambiguous metadata columns
normalise_master_types <- function(df) {
  num_cols <- intersect(
    c("PA", "tSPC",
      "AGB_pred", "AGB_low", "AGB_up", "AGB_CIwidt",
      "AGC_pred", "AGC_low", "AGC_up", "AGC_CIwidt",
      "carbon_index", "xcoord", "ycoord",
      "depth", "distToLand",
      "Ea_SPC","Th_SPC","Cr_SPC","Cs_SPC","Si_SPC","Hu_SPC",
      "Ho_SPC","Hp_SPC","Tc_SPC","Hm_SPC","Hs_SPC","Hd_SPC",
      paste0("GSE_A", sprintf("%02d", 0:63))),
    names(df))
  char_cols <- intersect(
    c("gee_id", "compositio", "sg_morpho", "morph3"),
    names(df))
  int_cols  <- intersect(c("year_gt"), names(df))
  
  for (col in num_cols)  df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  for (col in char_cols) df[[col]] <- as.character(df[[col]])
  for (col in int_cols)  df[[col]] <- suppressWarnings(as.integer(df[[col]]))
  df
}


result_dir <- here("result")
years      <- as.character(2017:2024)
dir.create(result_dir, showWarnings = FALSE)

# =============================================================================
# STEP 1: Predictors
# =============================================================================
env_vars      <- "depth" # c("depth", "distToLand")
gse_vars      <- paste0("GSE_A", sprintf("%02d", 0:63))
pa_predictors <- c(env_vars, gse_vars)
response_var  <- "PA"

# =============================================================================
# STEP 2: Load training data from master_raw_YYYY
# =============================================================================
cat("Loading PA training data from master_raw_YYYY...\n")

training_df <- purrr::map_dfr(years, function(yr) {
  path <- file.path(result_dir, paste0("master_raw_", yr, ".csv"))
  if (!file.exists(path)) { message(sprintf("master_raw_%s not found", yr)); return(NULL) }
  df <- readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
  df <- normalise_master_types(df)
  df %>%
    dplyr::filter(PA %in% c(0, 1)) %>%
    dplyr::mutate(year = as.factor(yr), PA = as.factor(PA))
})

cat(sprintf("\nTotal rows before NA filter: %d\n", nrow(training_df)))
cat("PA distribution:\n"); print(table(training_df$PA))
cat("Samples per year:\n")
print(training_df %>% count(year, name = "n_samples"))

# Optional class balancing
use_balanced <- TRUE
if (use_balanced) {
  n_pa1 <- sum(training_df$PA == "1")
  n_pa0 <- sum(training_df$PA == "0")
  if (n_pa0 > n_pa1) {
    df_0 <- training_df %>% filter(PA == "0") %>% sample_n(n_pa1)
    df_1 <- training_df %>% filter(PA == "1")
    training_df <- bind_rows(df_0, df_1)
    cat(sprintf("Balanced: PA=0 downsampled to %d (matched PA=1=%d)\n", n_pa1, n_pa1))
  }
}

# NA filter on predictors + response
n_before    <- nrow(training_df)
training_df <- training_df %>%
  filter(if_all(all_of(c(response_var, pa_predictors)), ~ !is.na(.)))
cat(sprintf("After NA filter: %d rows (removed %d)\n\n",
            nrow(training_df), n_before - nrow(training_df)))

# Save filtered training subset
write_csv(training_df, file.path(result_dir, "training_PA.csv"))
cat("Saved: training_PA.csv\n\n")

# =============================================================================
# STEP 3: Spatial map
# =============================================================================
x_col <- ifelse("xcoord" %in% names(training_df), "xcoord", "lon")
y_col <- ifelse("ycoord" %in% names(training_df), "ycoord", "lat")
if (all(c(x_col, y_col) %in% names(training_df))) {
  map_df <- training_df %>%
    mutate(x = as.numeric(.data[[x_col]]), y = as.numeric(.data[[y_col]])) %>%
    filter(!is.na(x), !is.na(y))
  print(
    ggplot() +
      geom_polygon(data = ggplot2::map_data("world", region = "Indonesia"),
                   aes(x = long, y = lat, group = group),
                   fill = "#d4e6c3", colour = "grey50", linewidth = 0.3) +
      geom_point(data = map_df,
                 aes(x = x, y = y,
                     colour = factor(PA, levels = c("0","1"),
                                     labels = c("Absent","Present"))),
                 size = 1, alpha = 0.4) +
      scale_colour_manual(values = c("Absent"="#e74c3c","Present"="#2980b9"),
                          name = "Seagrass PA") +
      coord_quickmap(xlim = range(map_df$x) + c(-1,1), ylim = range(map_df$y) + c(-1,1)) +
      labs(title = sprintf("PA Training Samples -- All Years (n = %d)", nrow(map_df)),
           x = "Longitude", y = "Latitude") +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
  )
}

# =============================================================================
# STEP 4: Hyperparameter tuning (cached)
# =============================================================================
grid_rds <- file.path(result_dir, "grid_rf_tuning_PA.rds")
tune_result <- if (file.exists(grid_rds)) {
  cat("Loading cached tuning...\n"); readRDS(grid_rds)
} else {
  cat("Running grid tuning...\n")
  res <- evaluate_rf_grid_ranger_classification(
    data = training_df, covars = pa_predictors, response = response_var,
    ntree_values = c(300, 500, 700), mtry_values = 4:10,
    nodesize_values = c(3, 5, 7), sample_fraction_values = c(0.6, 0.7, 0.8),
    n_iter = 5, train_frac = 0.7, seed = 42)
  saveRDS(res, grid_rds); res
}
write_csv(tune_result$grid_results, file.path(result_dir, "grid_rf_tuning_PA.csv"))
write_csv(tune_result$best_params,  file.path(result_dir, "grid_rf_best_PA.csv"))
best         <- tune_result$best_params
model_params <- list(ntree = best$ntree[1], mtry = best$mtry[1],
                     nodesize = best$nodesize[1], sample_fraction = best$sample_fraction[1])
cat("Best hyperparameters:\n"); print(model_params)

# =============================================================================
# STEP 5: Threshold evaluation
# =============================================================================
threshold_result <- evaluate_thresholds_classification(
  data = training_df, covars = pa_predictors, response = response_var,
  thresholds = seq(0.1, 0.9, 0.1), n_iter = 10,
  ntree = model_params$ntree, mtry = model_params$mtry)
write_csv(threshold_result, file.path(result_dir, "threshold_eval_PA.csv"))
auto_thresh           <- threshold_result$threshold[which.max(threshold_result$accuracy)]
model_params$threshold <- 0.6
cat(sprintf("Auto threshold: %.2f | Used: %.2f\n", auto_thresh, model_params$threshold))

# Threshold plot
print(
  ggplot(threshold_result, aes(x = threshold, y = accuracy)) +
    geom_line(colour = "grey40", linewidth = 1) + geom_point(size = 2) +
    geom_text(aes(label = round(accuracy, 3)), vjust = -0.8, size = 3.5) +
    geom_vline(xintercept = model_params$threshold, linetype = "dashed", colour = "red") +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(title = "PA Threshold vs Accuracy", x = "Threshold", y = "Accuracy") +
    theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold"))
)

# =============================================================================
# STEP 6: Train final model
# =============================================================================
final_rf <- randomForest(
  formula  = build_rf_formula(pa_predictors, response_var),
  data     = training_df[, c(response_var, pa_predictors)],
  ntree    = model_params$ntree,   mtry     = model_params$mtry,
  nodesize = model_params$nodesize,
  sampsize = floor(model_params$sample_fraction * nrow(training_df)),
  replace  = FALSE, importance = TRUE)
saveRDS(final_rf, file.path(result_dir, "final_rf_model_PA.rds"))
cat("Saved: final_rf_model_PA.rds\n")

# =============================================================================
# STEP 7: Predict on training data
# =============================================================================
training_df$PA_prob <- predict(final_rf, newdata = training_df, type = "prob")[, "1"]
training_df$PA_pred <- as.integer(training_df$PA_prob > model_params$threshold)

# =============================================================================
# STEP 8: Internal evaluation
# =============================================================================
cm <- caret::confusionMatrix(
  factor(training_df$PA_pred, levels = c(0,1)),
  factor(as.numeric(as.character(training_df$PA)), levels = c(0,1)), positive = "1")
cat("\n=== PA Internal Evaluation (Training) ===\n"); print(cm)

# =============================================================================
# STEP 9: Variable importance + predictor structure
# =============================================================================
var_imp_df <- as.data.frame(randomForest::importance(final_rf)) %>%
  tibble::rownames_to_column("variable") %>%
  select(variable, importance = MeanDecreaseGini) %>% arrange(desc(importance))
write_csv(var_imp_df, file.path(result_dir, "var_importance_rf_PA.csv"))

predictor_structure <- list(
  classes = sapply(training_df[, pa_predictors], class),
  levels  = lapply(training_df[, pa_predictors], function(x) if(is.factor(x)) levels(x) else NULL))
saveRDS(predictor_structure, file.path(result_dir, "rf_predictor_structure_PA.rds"))

# =============================================================================
# STEP 10: Monte Carlo
# =============================================================================
cat("\nRunning Monte Carlo (n=100)...\n")
mc_result  <- run_montecarlo_rf_classification(
  data = training_df, covars = pa_predictors, response = response_var,
  n_iter = 100, ntree = model_params$ntree, mtry = model_params$mtry,
  threshold = model_params$threshold, train_frac = 0.7, seed = 42)
write_csv(mc_result, file.path(result_dir, "rf_pa_mc_results.csv"))
summary_ci <- summarize_ci_classification(mc_result)
write_csv(summary_ci, file.path(result_dir, "rf_pa_mc_summary_ci.csv"))
cat("=== PA Monte Carlo CI ===\n"); print(summary_ci)

# =============================================================================
# STEP 11: Supplementary plots
# =============================================================================
print(
  ggplot(tune_result$grid_results,
         aes(x = factor(mtry), y = factor(nodesize), fill = accuracy_mean)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = round(accuracy_mean, 3)), size = 2.5, colour = "white") +
    scale_fill_gradient(low = "steelblue", high = "red", name = "Accuracy") +
    facet_wrap(~ ntree, labeller = label_both) +
    labs(title = "PA Grid Tuning Heatmap", x = "mtry", y = "nodesize") +
    theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold")))

print(var_imp_df %>% slice_head(n = 20) %>%
        ggplot(aes(x = reorder(variable, importance), y = importance)) +
        geom_col(fill = "grey40") + coord_flip() +
        labs(title = "Top 20 Variable Importance -- PA", x = NULL, y = "MeanDecreaseGini") +
        theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

print(mc_result %>%
        pivot_longer(c(accuracy, precision, recall, f1), names_to = "metric", values_to = "value") %>%
        ggplot(aes(x = value)) + geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
        facet_wrap(~ metric, scales = "free") +
        labs(title = "PA Monte Carlo Distributions (n=100)", x = "Value", y = "Count") +
        theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold")))

cat("\nA1_rf_model_PA.R COMPLETE\n")
