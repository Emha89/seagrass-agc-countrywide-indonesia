# =============================================================================
# C1_rf_model_cIndex.R
# RF regression -- carbon index (0-1).
#
# Reads: master_raw_YYYY.csv (carbon_index already in master from
# 01_build_master_raw.R)
# Saves:
#   final_rf_model_CINDEX.rds
#   rf_predictor_structure_CINDEX.rds
#   training_CINDEX.csv
#   var_importance_rf_CINDEX.csv
#   grid_rf_tuning_CINDEX.csv / grid_rf_best_CINDEX.csv
#   rf_cindex_mc_results.csv / rf_cindex_mc_summary_ci.csv
#
# PREDICTOR NOTE: env_vars is depth only here too -- third stage now
# (after A1_rf_model_PA.R and B1_rf_model_leafMorpho.R) confirming this
# as a consistent national-scale predictor choice, not stage-specific.
#
# TRAINING FILTER: rows with carbon_index == 0 are excluded from
# training ("structural zeros excluded -- modelling carbon presence",
# per the original script's own comment) -- this model is trained only
# on locations with positive carbon_index.
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
env_vars          <- "depth" # c("depth", "distToLand")
gse_vars          <- paste0("GSE_A", sprintf("%02d", 0:63))
cindex_predictors <- c(env_vars, gse_vars)
response_var      <- "carbon_index"

# =============================================================================
# STEP 2: Load training data -- rows where carbon_index > 0
# carbon_index is already in master_raw (derived in 01_build_master_raw.R)
# =============================================================================
cat("Loading CINDEX training data from master_raw_YYYY...\n")

training_df <- purrr::map_dfr(years, function(yr) {
  path <- file.path(result_dir, paste0("master_raw_", yr, ".csv"))
  if (!file.exists(path)) { message(sprintf("master_raw_%s not found", yr)); return(NULL) }
  df <- readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
  df <- normalise_master_types(df)
  df %>%
    dplyr::mutate(year = as.factor(yr))
})

# Distribution report before filtering
n_total   <- nrow(training_df)
n_zero    <- sum(training_df$carbon_index == 0, na.rm = TRUE)
n_nonzero <- sum(training_df$carbon_index  > 0, na.rm = TRUE)
n_na      <- sum(is.na(training_df$carbon_index))
cat(sprintf("\ncarbon_index distribution:\n  Total: %d | NA: %d | = 0: %d (%.1f%%) | > 0: %d (%.1f%%)\n",
            n_total, n_na, n_zero, 100*n_zero/n_total, n_nonzero, 100*n_nonzero/n_total))
cat(sprintf("  Min: %.4f | Max: %.4f | Mean: %.4f | Median: %.4f\n",
            min(training_df$carbon_index, na.rm = TRUE),
            max(training_df$carbon_index, na.rm = TRUE),
            mean(training_df$carbon_index, na.rm = TRUE),
            median(training_df$carbon_index, na.rm = TRUE)))

cat("\nPer-year breakdown:\n")
print(training_df %>%
  group_by(year) %>%
  summarise(n = n(),
            n_zero    = sum(carbon_index == 0, na.rm = TRUE),
            pct_zero  = round(n_zero/n*100, 1),
            mean_ci   = round(mean(carbon_index, na.rm = TRUE), 4),
            .groups   = "drop"))

# Filter: retain carbon_index > 0 (structural zeros excluded -- modelling carbon presence)
training_df <- training_df %>% filter(carbon_index > 0, carbon_index <= 1)
cat(sprintf("\nAfter zero-filter: %d rows retained\n", nrow(training_df)))
cat("Samples per year:\n")
print(training_df %>% count(year, name = "n_samples"))

# NA filter on predictors + response
n_before    <- nrow(training_df)
training_df <- training_df %>%
  filter(if_all(all_of(c(response_var, cindex_predictors)), ~ !is.na(.)))
cat(sprintf("After NA filter: %d rows (removed %d)\n\n",
            nrow(training_df), n_before - nrow(training_df)))

write_csv(training_df, file.path(result_dir, "training_CINDEX.csv"))
cat("Saved: training_CINDEX.csv\n\n")

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
                 fill = "grey90", colour = "grey50", linewidth = 0.3) +
    geom_point(data = map_df, aes(x = x, y = y, colour = carbon_index), size = 1, alpha = 0.6) +
    scale_colour_gradient(name = "Carbon Index", low = "lightgreen", high = "darkgreen") +
    coord_quickmap(xlim = range(map_df$x) + c(-1,1), ylim = range(map_df$y) + c(-1,1)) +
    labs(title = sprintf("CINDEX Training Samples -- All Years (n = %d)", nrow(map_df)),
         x = "Longitude", y = "Latitude") +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom"))
}

# =============================================================================
# STEP 4: Hyperparameter tuning (cached)
# =============================================================================
grid_rds <- file.path(result_dir, "grid_rf_tuning_CINDEX.rds")
tune_result <- if (file.exists(grid_rds)) {
  cat("Loading cached tuning...\n"); readRDS(grid_rds)
} else {
  cat("Running grid tuning...\n")
  res <- evaluate_rf_grid_ranger_regression(
    data = training_df, covars = cindex_predictors, response = response_var,
    ntree_values = c(300,500,700,900), mtry_values = 4:10,
    nodesize_values = c(3,5,7,9), sample_fraction_values = c(0.6,0.7,0.8,0.9),
    n_iter = 5, train_frac = 0.7, seed = 42)
  saveRDS(res, grid_rds); res
}
write_csv(tune_result$grid_results, file.path(result_dir, "grid_rf_tuning_CINDEX.csv"))
write_csv(tune_result$best_params,  file.path(result_dir, "grid_rf_best_CINDEX.csv"))
best         <- tune_result$best_params
model_params <- list(ntree = best$ntree[1], mtry = best$mtry[1],
                     nodesize = best$nodesize[1], sample_fraction = best$sample_fraction[1])
cat("Best hyperparameters:\n"); print(model_params)

# =============================================================================
# STEP 5: Train final model
# =============================================================================
final_rf <- randomForest(
  formula = build_rf_formula(cindex_predictors, response_var),
  data    = training_df[, c(response_var, cindex_predictors)],
  ntree = model_params$ntree, mtry = model_params$mtry, nodesize = model_params$nodesize,
  sampsize = floor(model_params$sample_fraction * nrow(training_df)),
  replace = FALSE, importance = TRUE)
saveRDS(final_rf, file.path(result_dir, "final_rf_model_CINDEX.rds"))
cat("Saved: final_rf_model_CINDEX.rds\n")

# =============================================================================
# STEP 6: Predict + internal evaluation
# =============================================================================
training_df$carbon_index_pred <- pmin(pmax(predict(final_rf, newdata = training_df), 0), 1)
cat("\n=== Carbon Index Internal Evaluation (Training) ===\n")
cat(sprintf("RMSE: %.4f\n", caret::RMSE(training_df$carbon_index_pred, training_df$carbon_index)))
cat(sprintf("MAE : %.4f\n", caret::MAE(training_df$carbon_index_pred,  training_df$carbon_index)))
cat(sprintf("R2  : %.4f\n", caret::R2(training_df$carbon_index_pred,   training_df$carbon_index)))

# =============================================================================
# STEP 7: Variable importance + predictor structure
# =============================================================================
var_imp_df <- as.data.frame(final_rf$importance) %>%
  tibble::rownames_to_column("variable") %>%
  select(variable, importance = IncNodePurity) %>% arrange(desc(importance))
write_csv(var_imp_df, file.path(result_dir, "var_importance_rf_CINDEX.csv"))
predictor_structure <- list(
  classes = sapply(training_df[, cindex_predictors], class),
  levels  = lapply(training_df[, cindex_predictors], function(x) if(is.factor(x)) levels(x) else NULL))
saveRDS(predictor_structure, file.path(result_dir, "rf_predictor_structure_CINDEX.rds"))

# =============================================================================
# STEP 8: Monte Carlo
# =============================================================================
cat("\nRunning Monte Carlo (n=100)...\n")
mc_result  <- run_montecarlo_rf_regression(
  data = training_df, covars = cindex_predictors, response = response_var,
  n_iter = 100, ntree = model_params$ntree, mtry = model_params$mtry,
  nodesize = model_params$nodesize, sample_fraction = model_params$sample_fraction,
  train_frac = 0.7, seed = 42)
write_csv(mc_result, file.path(result_dir, "rf_cindex_mc_results.csv"))
summary_ci <- summarize_ci_regression(mc_result)
write_csv(summary_ci, file.path(result_dir, "rf_cindex_mc_summary_ci.csv"))
cat("=== CINDEX Monte Carlo CI ===\n"); print(summary_ci)

# =============================================================================
# STEP 9: Supplementary plots
# =============================================================================
print(ggplot(tune_result$grid_results,
    aes(x = factor(mtry), y = factor(nodesize), fill = rmse_median)) +
  geom_tile(colour = "white") + geom_text(aes(label = round(rmse_median, 4)), size = 2.3, colour = "white") +
  scale_fill_gradient(low = "steelblue", high = "red") +
  facet_wrap(~ ntree, labeller = label_both) +
  labs(title = "CINDEX Grid Tuning Heatmap", x = "mtry", y = "nodesize") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold")))

print(var_imp_df %>% slice_head(n = 20) %>%
  ggplot(aes(x = reorder(variable, importance), y = importance)) +
  geom_col(fill = "grey40") + coord_flip() +
  labs(title = "Top 20 Variable Importance -- CINDEX", x = NULL, y = "IncNodePurity") +
  theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

print(ggplot(training_df, aes(x = carbon_index, y = carbon_index_pred)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  labs(title = "Carbon Index: Observed vs Predicted (Training)",
       x = "Observed", y = "Predicted") +
  theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

print(mc_result %>%
  pivot_longer(c(RMSE, MAE, R2), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = value)) + geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~ metric, scales = "free") +
  labs(title = "CINDEX Monte Carlo Distributions (n=100)", x = "Value", y = "Count") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold")))

cat("\nC1_rf_model_cIndex.R COMPLETE\n")
