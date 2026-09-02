# =============================================================================
# B1_rf_model_leafMorpho.R
# RF classification -- leaf morphology (3-class: morph3).
#
# Reads: master_raw_YYYY.csv (morph3 column already derived in
# 01_build_master_raw.R)
# Saves:
#   final_rf_model_MORPH3.rds
#   rf_predictor_structure_MORPH3.rds
#   training_MORPH3.csv
#   training_MORPH3_points_ALL.csv  -- GEE upload reference points
#   var_importance_rf_MORPH3.csv
#   confusion_matrix_MORPH3.csv
#   grid_rf_tuning_MORPH3.csv / grid_rf_best_MORPH3.csv
#   rf_morph3_mc_results.csv / rf_morph3_mc_summary_ci.csv
#
# PREDICTOR NOTE: env_vars is depth only here too (same as
# A1_rf_model_PA.R) -- distToLand present as a commented-out alternative
# but not active. Two stages now showing this same choice reinforces that
# it's a deliberate national-scale predictor decision, not specific to
# one stage.
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
set.seed(42)

# =============================================================================
# STEP 1: Predictors
# =============================================================================
env_vars         <- "depth" #  c("depth", "distToLand")
gse_vars         <- paste0("GSE_A", sprintf("%02d", 0:63))
morph_predictors <- c(env_vars, gse_vars)
response_var     <- "morph3"

# =============================================================================
# STEP 2: Load training data -- rows where morph3 is not NA
# =============================================================================
cat("Loading MORPH3 training data from master_raw_YYYY...\n")

training_df <- purrr::map_dfr(years, function(yr) {
  path <- file.path(result_dir, paste0("master_raw_", yr, ".csv"))
  if (!file.exists(path)) { message(sprintf("master_raw_%s not found", yr)); return(NULL) }
  df <- readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
  df <- normalise_master_types(df)
  df %>%
    dplyr::filter(!is.na(morph3)) %>%
    dplyr::mutate(year = as.factor(yr), morph3 = as.factor(morph3))
})

cat(sprintf("\nMORPH3 training rows: %d\n", nrow(training_df)))
cat("Samples per year:\n")
print(training_df %>% count(year, name = "n_samples"))
cat("Class distribution:\n")
print(training_df %>% count(morph3) %>% mutate(prop = round(n / sum(n), 3)))

# NA filter
n_before    <- nrow(training_df)
training_df <- training_df %>%
  filter(if_all(all_of(c(response_var, morph_predictors)), ~ !is.na(.)))
cat(sprintf("After NA filter: %d rows (removed %d)\n\n",
            nrow(training_df), n_before - nrow(training_df)))

write_csv(training_df, file.path(result_dir, "training_MORPH3.csv"))
cat("Saved: training_MORPH3.csv\n\n")

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
    geom_point(data = map_df, aes(x = x, y = y, colour = morph3), size = 1, alpha = 0.5) +
    coord_quickmap(xlim = range(map_df$x) + c(-1,1), ylim = range(map_df$y) + c(-1,1)) +
    labs(title = sprintf("MORPH3 Training Samples -- All Years (n = %d)", nrow(map_df)),
         x = "Longitude", y = "Latitude", colour = "Morphology") +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom"))
}

# =============================================================================
# STEP 4: Export GEE reference points
# =============================================================================
coord_cols <- intersect(c("xcoord","ycoord","lon","lat"), names(training_df))
gee_ref <- training_df %>%
  select(gee_id, year, morph3, any_of(coord_cols)) %>%
  mutate(year = as.integer(as.character(year)), morph3 = as.character(morph3)) %>%
  distinct(gee_id, .keep_all = TRUE)
write_csv(gee_ref, file.path(result_dir, "training_MORPH3_points_ALL.csv"))
cat(sprintf("Saved: training_MORPH3_points_ALL.csv (%d rows)\n\n", nrow(gee_ref)))

# =============================================================================
# STEP 5: Hyperparameter tuning (cached)
# =============================================================================
grid_rds <- file.path(result_dir, "grid_rf_tuning_MORPH3.rds")
tune_result <- if (file.exists(grid_rds)) {
  cat("Loading cached tuning...\n"); readRDS(grid_rds)
} else {
  cat("Running grid tuning...\n")
  res <- evaluate_rf_grid_ranger_multiclass(
    data = training_df, covars = morph_predictors, response = response_var,
    ntree_values = c(300,500,700,900), mtry_values = 4:10,
    nodesize_values = c(3,5,7,9), sample_fraction_values = c(0.6,0.7,0.8,0.9),
    n_iter = 5, train_frac = 0.7, seed = 42)
  saveRDS(res, grid_rds); res
}
write_csv(tune_result$grid_results, file.path(result_dir, "grid_rf_tuning_MORPH3.csv"))
write_csv(tune_result$best_params,  file.path(result_dir, "grid_rf_best_MORPH3.csv"))
best         <- tune_result$best_params
model_params <- list(ntree = best$ntree[1], mtry = best$mtry[1],
                     nodesize = best$nodesize[1], sample_fraction = best$sample_fraction[1])
cat("Best hyperparameters:\n"); print(model_params)

# =============================================================================
# STEP 6: Train final model
# =============================================================================
final_rf <- randomForest(
  formula = build_rf_formula(morph_predictors, response_var),
  data    = training_df[, c(response_var, morph_predictors)],
  ntree = model_params$ntree, mtry = model_params$mtry, nodesize = model_params$nodesize,
  sampsize = floor(model_params$sample_fraction * nrow(training_df)),
  replace = FALSE, importance = TRUE)
saveRDS(final_rf, file.path(result_dir, "final_rf_model_MORPH3.rds"))
cat("Saved: final_rf_model_MORPH3.rds\n")

# =============================================================================
# STEP 7: Predict class probabilities + internal evaluation
# =============================================================================
prob_mat          <- predict(final_rf, newdata = training_df, type = "prob")
prob_df           <- as.data.frame(prob_mat)
colnames(prob_df) <- paste0("P_", colnames(prob_df))
training_df_probs <- bind_cols(training_df, prob_df)
write_csv(training_df_probs, file.path(result_dir, "training_data_for_GEE_MORPH3_probs.csv"))

pred_class <- predict(final_rf, newdata = training_df)
cm         <- caret::confusionMatrix(pred_class, training_df$morph3)
cat("\n=== MORPH3 Internal Evaluation (Training) ===\n")
cat(sprintf("Accuracy: %.4f | Kappa: %.4f\n", cm$overall["Accuracy"], cm$overall["Kappa"]))
cat("Class-wise accuracy:\n"); print(round(prop.table(cm$table, margin = 2), 3))
write_csv(as.data.frame(cm$table), file.path(result_dir, "confusion_matrix_MORPH3.csv"))

# =============================================================================
# STEP 8: Variable importance + predictor structure
# =============================================================================
var_imp_df <- as.data.frame(final_rf$importance) %>%
  tibble::rownames_to_column("variable") %>% arrange(desc(MeanDecreaseGini))
write_csv(var_imp_df, file.path(result_dir, "var_importance_rf_MORPH3.csv"))
predictor_structure <- list(
  classes = sapply(training_df[, morph_predictors], class),
  levels  = lapply(training_df[, morph_predictors], function(x) if(is.factor(x)) levels(x) else NULL))
saveRDS(predictor_structure, file.path(result_dir, "rf_predictor_structure_MORPH3.rds"))

# =============================================================================
# STEP 9: Monte Carlo
# =============================================================================
cat("\nRunning Monte Carlo (n=100)...\n")
mc_result  <- run_montecarlo_rf_multiclass(
  data = training_df, covars = morph_predictors, response = response_var,
  n_iter = 100, ntree = model_params$ntree, mtry = model_params$mtry, seed = 42)
write_csv(mc_result, file.path(result_dir, "rf_morph3_mc_results.csv"))
summary_ci <- mc_result %>% summarise(
  accuracy_mean = mean(accuracy), accuracy_sd = sd(accuracy),
  accuracy_lower = quantile(accuracy, 0.025), accuracy_upper = quantile(accuracy, 0.975),
  kappa_mean = mean(kappa), kappa_sd = sd(kappa),
  kappa_lower = quantile(kappa, 0.025), kappa_upper = quantile(kappa, 0.975))
write_csv(summary_ci, file.path(result_dir, "rf_morph3_mc_summary_ci.csv"))
cat("=== MORPH3 Monte Carlo CI ===\n"); print(summary_ci)

# =============================================================================
# STEP 10: Supplementary plots
# =============================================================================
print(ggplot(tune_result$grid_results,
    aes(x = factor(mtry), y = factor(nodesize), fill = accuracy_mean)) +
  geom_tile(colour = "white") + geom_text(aes(label = round(accuracy_mean, 3)), size = 2.5, colour = "white") +
  scale_fill_gradient(low = "steelblue", high = "red", name = "Accuracy") +
  facet_wrap(~ ntree, labeller = label_both) +
  labs(title = "MORPH3 Grid Tuning Heatmap", x = "mtry", y = "nodesize") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold")))

print(var_imp_df %>% slice_head(n = 20) %>%
  ggplot(aes(x = reorder(variable, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_col(fill = "grey40") + coord_flip() +
  labs(title = "Top 20 Variable Importance -- MORPH3", x = NULL, y = "MeanDecreaseGini") +
  theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

print(ggplot(as.data.frame(cm$table), aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(colour = "white") + geom_text(aes(label = Freq), size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "MORPH3 Confusion Matrix (Training)", x = "Actual", y = "Predicted") +
  theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

print(mc_result %>%
  pivot_longer(c(accuracy, kappa), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = value)) + geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~ metric, scales = "free") +
  labs(title = "MORPH3 Monte Carlo Distributions (n=100)", x = "Value", y = "Count") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold")))

cat("\nB1_rf_model_leafMorpho.R COMPLETE\n")
