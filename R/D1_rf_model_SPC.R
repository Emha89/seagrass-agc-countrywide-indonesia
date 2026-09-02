# =============================================================================
# D1_rf_model_SPC.R
# RF regression -- seagrass percent cover (SPC / tSPC).
#
# Reads: master_raw_YYYY.csv  (from 01_build_master_raw.R)
# Saves:
#   final_rf_model_PCT.rds
#   rf_predictor_structure_PCT.rds
#   training_SPC.csv
#   var_importance_rf_PCT.csv
#   grid_rf_tuning_PCT.csv / grid_rf_best_PCT.csv
#   rf_pct_mc_results.csv / rf_pct_mc_summary_ci.csv
#
# DEPENDENCY: morph_pred (the P_* predictors) are NOT in master_raw --
# they come from B2_rf_applyModel_leafMorpho.R's output
# (predicted_MORPH3_probs_<year>.csv). This script must therefore run
# AFTER B2 has produced that file for each year; if it doesn't exist yet,
# P_* stays NA and those rows get dropped in the NA filter below,
# leaving little or no usable training data.
#
# Grid tuning search here (ntree 700/900, mtry 8/10) is narrower than the
# other stages (which search ntree 300-900, mtry 4:10) -- looks like a
# deliberately refined range rather than the full default grid.
#
# Data availability: this script's inputs are not included in this
# repository.
# =============================================================================

library(dplyr); library(readr); library(randomForest)
library(caret);  library(ranger); library(purrr)
library(ggplot2); library(tibble); library(tidyr); library(maps)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "func_reg_rf.R"))


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
env_vars       <- "depth" # c("depth", "distToLand")
gse_vars       <- paste0("GSE_A", sprintf("%02d", 0:63))
morph_pred     <- c("P_mixed_short_plus_mono_short", "P_mixed_long", "P_mono_Ea")
pct_predictors <- c(env_vars, gse_vars, morph_pred)
response_var   <- "tSPC"

# =============================================================================
# STEP 2: Load training data from master_raw_YYYY
# NOTE: morph_pred (P_*) are NOT in master_raw -- they come from
# B2_rf_applyModel_leafMorpho.R. If predicted_MORPH3_probs_YYYY exists,
# join it; otherwise P_* remain NA and rows with NA predictors are
# removed in the NA filter step below.
# =============================================================================
cat("Loading SPC training data from master_raw_YYYY...\n")

training_df <- purrr::map_dfr(years, function(yr) {
  path <- file.path(result_dir, paste0("master_raw_", yr, ".csv"))
  if (!file.exists(path)) { message(sprintf("master_raw_%s not found", yr)); return(NULL) }
  df_raw <- readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
  df_raw <- normalise_master_types(df_raw)
  df <- df_raw %>%
    dplyr::filter(!is.na(tSPC), tSPC >= 0, tSPC <= 100) %>%
    dplyr::mutate(year = as.factor(yr))

  # Join Morph3 probabilities if available
  morph_path <- file.path(result_dir, paste0("predicted_MORPH3_probs_", yr, ".csv"))
  if (file.exists(morph_path)) {
    df_morph <- readr::read_csv(morph_path, show_col_types = FALSE) %>%
      dplyr::select(gee_id, dplyr::any_of(morph_pred))
    df <- dplyr::left_join(df, df_morph, by = "gee_id")
  } else {
    for (mc in morph_pred) if (!mc %in% names(df)) df[[mc]] <- NA_real_
  }
  df
})

cat(sprintf("\nRows with tSPC: %d\n", nrow(training_df)))
cat("Samples per year:\n")
print(training_df %>% count(year, name = "n_samples"))

# NA filter on predictors + response
n_before    <- nrow(training_df)
req_cols    <- intersect(c(response_var, pct_predictors), names(training_df))
training_df <- training_df %>%
  filter(if_all(all_of(req_cols), ~ !is.na(.)))
cat(sprintf("After NA filter: %d rows (removed %d)\n\n",
            nrow(training_df), n_before - nrow(training_df)))

write_csv(training_df, file.path(result_dir, "training_SPC.csv"))
cat("Saved: training_SPC.csv\n\n")

# =============================================================================
# STEP 3: Spatial map
# =============================================================================
x_col <- ifelse("xcoord" %in% names(training_df), "xcoord", "lon")
y_col <- ifelse("ycoord" %in% names(training_df), "ycoord", "lat")
if (all(c(x_col, y_col) %in% names(training_df))) {
  map_df <- training_df %>%
    mutate(x = as.numeric(.data[[x_col]]), y = as.numeric(.data[[y_col]])) %>% filter(!is.na(x))
  print(ggplot() +
    geom_polygon(data = ggplot2::map_data("world", region = "Indonesia"),
                 aes(x = long, y = lat, group = group),
                 fill = "#d4e6c3", colour = "grey50", linewidth = 0.3) +
    geom_point(data = map_df, aes(x = x, y = y, colour = tSPC), size = 1, alpha = 0.5) +
    scale_colour_gradient(name = "tSPC (%)", low = "lightyellow", high = "darkgreen") +
    coord_quickmap(xlim = range(map_df$x) + c(-1,1), ylim = range(map_df$y) + c(-1,1)) +
    labs(title = sprintf("SPC Training Samples -- All Years (n = %d)", nrow(map_df)),
         x = "Longitude", y = "Latitude") +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom"))
}

# =============================================================================
# STEP 4: Hyperparameter tuning (cached)
# =============================================================================
grid_rds <- file.path(result_dir, "grid_rf_tuning_PCT.rds")
tune_result <- if (file.exists(grid_rds)) {
  cat("Loading cached tuning...\n"); readRDS(grid_rds)
} else {
  cat("Running grid tuning...\n")
  res <- evaluate_rf_grid_ranger_regression(
    data = training_df, covars = pct_predictors, response = response_var,
    ntree_values = c(700, 900), mtry_values = c(8, 10),
    nodesize_values = c(7, 9), sample_fraction_values = c(0.7, 0.8, 0.9),
    n_iter = 5, train_frac = 0.7, seed = 42)
  saveRDS(res, grid_rds); res
}
write_csv(tune_result$grid_results, file.path(result_dir, "grid_rf_tuning_PCT.csv"))
write_csv(tune_result$best_params,  file.path(result_dir, "grid_rf_best_PCT.csv"))
best         <- tune_result$best_params
model_params <- list(ntree = best$ntree[1], mtry = best$mtry[1],
                     nodesize = best$nodesize[1], sample_fraction = best$sample_fraction[1])
cat("Best hyperparameters:\n"); print(model_params)

# =============================================================================
# STEP 5: Train final model
# =============================================================================
final_rf <- randomForest(
  formula  = build_rf_formula(pct_predictors, response_var),
  data     = training_df[, c(response_var, pct_predictors)],
  ntree = model_params$ntree, mtry = model_params$mtry, nodesize = model_params$nodesize,
  sampsize = floor(model_params$sample_fraction * nrow(training_df)),
  replace = FALSE, importance = TRUE)
saveRDS(final_rf, file.path(result_dir, "final_rf_model_PCT.rds"))
cat("Saved: final_rf_model_PCT.rds\n")

# =============================================================================
# STEP 6: Predict + internal evaluation
# =============================================================================
training_df$tSPC_pred <- predict(final_rf, newdata = training_df)
cat("\n=== SPC Internal Evaluation (Training) ===\n")
cat(sprintf("RMSE: %.3f\n", caret::RMSE(training_df$tSPC_pred, training_df$tSPC)))
cat(sprintf("MAE : %.3f\n", caret::MAE(training_df$tSPC_pred,  training_df$tSPC)))
cat(sprintf("R2  : %.3f\n", caret::R2(training_df$tSPC_pred,   training_df$tSPC)))

# =============================================================================
# STEP 7: Variable importance + predictor structure
# =============================================================================
var_imp_df <- as.data.frame(final_rf$importance) %>%
  tibble::rownames_to_column("variable") %>%
  select(variable, importance = IncNodePurity) %>% arrange(desc(importance))
write_csv(var_imp_df, file.path(result_dir, "var_importance_rf_PCT.csv"))
predictor_structure <- list(
  classes = sapply(training_df[, pct_predictors], class),
  levels  = lapply(training_df[, pct_predictors], function(x) if(is.factor(x)) levels(x) else NULL))
saveRDS(predictor_structure, file.path(result_dir, "rf_predictor_structure_PCT.rds"))

# =============================================================================
# STEP 8: Monte Carlo
# =============================================================================
cat("\nRunning Monte Carlo (n=100)...\n")
mc_result  <- run_montecarlo_rf_regression(
  data = training_df, covars = pct_predictors, response = response_var,
  n_iter = 100, ntree = model_params$ntree, mtry = model_params$mtry,
  nodesize = model_params$nodesize, sample_fraction = model_params$sample_fraction,
  train_frac = 0.7, seed = 42)
write_csv(mc_result, file.path(result_dir, "rf_pct_mc_results.csv"))
summary_ci <- summarize_ci_regression(mc_result)
write_csv(summary_ci, file.path(result_dir, "rf_pct_mc_summary_ci.csv"))
cat("=== SPC Monte Carlo CI ===\n"); print(summary_ci)

# =============================================================================
# STEP 9: Supplementary plots
# =============================================================================
print(ggplot(tune_result$grid_results,
    aes(x = factor(mtry), y = factor(nodesize), fill = rmse_median)) +
  geom_tile(colour = "white") + geom_text(aes(label = round(rmse_median, 2)), size = 2.5, colour = "white") +
  scale_fill_gradient(low = "steelblue", high = "red", name = "RMSE") +
  facet_wrap(~ ntree, labeller = label_both) +
  labs(title = "SPC Grid Tuning Heatmap", x = "mtry", y = "nodesize") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold")))

print(var_imp_df %>% slice_head(n = 20) %>%
  ggplot(aes(x = reorder(variable, importance), y = importance)) +
  geom_col(fill = "grey40") + coord_flip() +
  labs(title = "Top 20 Variable Importance -- SPC", x = NULL, y = "IncNodePurity") +
  theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

print(ggplot(training_df, aes(x = tSPC, y = tSPC_pred)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  labs(title = "SPC: Observed vs Predicted (Training)",
       x = "Observed tSPC (%)", y = "Predicted tSPC (%)") +
  theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

print(mc_result %>%
  pivot_longer(c(RMSE, MAE, R2), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = value)) + geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~ metric, scales = "free") +
  labs(title = "SPC Monte Carlo Distributions (n=100)", x = "Value", y = "Count") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold")))

cat("\nD1_rf_model_SPC.R COMPLETE\n")
