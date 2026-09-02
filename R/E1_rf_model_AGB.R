# =============================================================================
# E1_rf_model_AGB.R
# RF regression -- above-ground biomass (AGB).
#
# Reads: master_AGB_YYYY.csv (from 04_build_master_AGB_AGC.R)
#        Contains master_raw columns + tSPC_pred
# Saves:
#   final_rf_model_AGB.rds
#   rf_predictor_structure_AGB.rds
#   training_AGB.csv
#   var_importance_rf_AGB.csv
#   grid_rf_tuning_AGB.csv / grid_rf_best_AGB.csv
#   rf_agb_mc_results.csv / rf_agb_mc_summary_ci.csv
#
# BOOTSTRAP UNCERTAINTY: AGB_simulated is drawn uniformly between AGB_low
# and AGB_up per row -- training_df_final uses one draw (for the final
# saved model), while training_df_boot expands each row into n_boot=10
# independent draws (for Monte Carlo uncertainty quantification only).
# This matches the "10 bootstrap-expanded iterations for AGB/AGC" pattern
# used throughout this chapter, distinct from the 100-iteration Monte
# Carlo used for PA/morphology/SPC/carbon index.
#
# TUNING EFFICIENCY: hyperparameter search subsamples to 30% of the
# training data (slice_sample) and uses n_iter=3, fewer than other
# stages -- likely a deliberate runtime optimisation given AGB's
# national-scale training set size.
#
# PREDICTOR NOTE: env_vars is depth only, consistent with every other
# stage in this repository.
#
# Data availability: this script's inputs are not included in this
# repository.
# =============================================================================

library(dplyr); library(readr); library(randomForest)
library(caret);  library(ranger); library(purrr)
library(ggplot2); library(tibble); library(tidyr); library(maps)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "func_reg_rf.R"))


# Helper: normalise column types that vary across CSV files (prevents
# bind_rows type conflicts). Covers all columns present in master_raw,
# master_AGB, and master_AGC.
normalise_master_types <- function(df) {
  num_cols <- c(
    # GT response / reference columns
    "PA", "tSPC", "AGB_pred", "AGB_low", "AGB_up", "AGB_CIwidt",
    "AGC_pred", "AGC_low", "AGC_up", "AGC_CIwidt", "carbon_index",
    # Coordinates
    "xcoord", "ycoord",
    # Environmental predictors
    "depth", "distToLand", "rugosity", "slope",
    "elevation", "mean_wave_period", "sig_wave_height",
    # Species % cover columns
    "Ea_SPC","Th_SPC","Cr_SPC","Cs_SPC","Si_SPC","Hu_SPC",
    "Ho_SPC","Hp_SPC","Tc_SPC","Hm_SPC","Hs_SPC","Hd_SPC",
    # Predicted columns added by apply/model scripts
    "tSPC_pred", "tAGB_pred", "tAGC_pred",
    "PA_prob", "carbon_index_pred",
    "P_mixed_short_plus_mono_short", "P_mixed_long", "P_mono_Ea",
    # GSE spectral bands
    paste0("GSE_A", sprintf("%02d", 0:63))
  )
  char_cols <- c("OBJECTID", "gee_id", "compositio", "location",
                 "sg_morpho", "morph3")
  
  for (col in intersect(num_cols,  names(df)))
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  for (col in intersect(char_cols, names(df)))
    df[[col]] <- as.character(df[[col]])
  if ("year_gt" %in% names(df))
    df$year_gt <- suppressWarnings(as.integer(df$year_gt))
  if ("PA_pred" %in% names(df))
    df$PA_pred <- suppressWarnings(as.integer(df$PA_pred))
  df
}


source(here("R", "04_build_master_AGB_AGC.R"))   # provides prepare_agb_training()

result_dir <- here("result")
dir.create(result_dir, showWarnings = FALSE)

# =============================================================================
# STEP 1: Predictors
# =============================================================================
env_vars       <- "depth" # c("depth", "distToLand")
gse_vars       <- paste0("GSE_A", sprintf("%02d", 0:63))
agb_predictors <- c(env_vars, gse_vars)
response_var   <- "AGB_simulated"

# =============================================================================
# STEP 2: Load training data
# =============================================================================
training_df <- prepare_agb_training(filter = "noFilter")

training_df <- training_df %>%
  filter(!is.na(AGB_pred), AGB_pred > 0,
         !is.na(AGB_low), !is.na(AGB_up),
         AGB_low >= 0, AGB_up > AGB_low)

cat(sprintf("Rows with valid AGB CI: %d\n", nrow(training_df)))

# Bootstrap data (for Monte Carlo with error propagation)
n_boot <- 10; set.seed(42)
training_df_boot <- training_df %>%
  slice(rep(1:n(), each = n_boot)) %>%
  mutate(sim_id = rep(1:n_boot, times = nrow(training_df)),
         AGB_simulated = runif(n(), AGB_low, AGB_up)) %>%
  filter(!is.na(AGB_simulated))
cat(sprintf("Bootstrap data: %d rows (%d orig x %d draws)\n", nrow(training_df_boot), nrow(training_df), n_boot))

# Single draw for final model
set.seed(42)
training_df_final <- training_df %>% rowwise() %>%
  mutate(AGB_simulated = runif(1, AGB_low, AGB_up)) %>% ungroup()

# NA filter
n_before          <- nrow(training_df_final)
training_df_final <- training_df_final %>%
  filter(if_all(all_of(c(response_var, agb_predictors)), ~ !is.na(.)))
cat(sprintf("After NA filter: %d rows (removed %d)\n", nrow(training_df_final), n_before - nrow(training_df_final)))
cat("Samples per year:\n")
yr_col <- intersect(c("year","year_gt"), names(training_df_final))[1]
print(training_df_final %>% count(.data[[yr_col]], name = "n_samples"))

write_csv(training_df_final, file.path(result_dir, "training_AGB.csv"))
cat("Saved: training_AGB.csv\n\n")

# =============================================================================
# STEP 3: Spatial map
# =============================================================================
x_col <- ifelse("xcoord" %in% names(training_df_final), "xcoord", "lon")
y_col <- ifelse("ycoord" %in% names(training_df_final), "ycoord", "lat")
if (all(c(x_col, y_col) %in% names(training_df_final))) {
  map_df <- training_df_final %>%
    mutate(x = as.numeric(.data[[x_col]]), y = as.numeric(.data[[y_col]])) %>% filter(!is.na(x))
  print(ggplot() +
          geom_polygon(data = ggplot2::map_data("world", region = "Indonesia"),
                       aes(x = long, y = lat, group = group),
                       fill = "#d4e6c3", colour = "grey50", linewidth = 0.3) +
          geom_point(data = map_df, aes(x = x, y = y, colour = AGB_pred), size = 1, alpha = 0.5) +
          scale_colour_gradient(name = "AGB (gDW/m2)", low = "lightblue", high = "darkblue") +
          coord_quickmap(xlim = range(map_df$x) + c(-1,1), ylim = range(map_df$y) + c(-1,1)) +
          labs(title = sprintf("AGB Training Samples -- All Years (n = %d)", nrow(map_df)),
               x = "Longitude", y = "Latitude") +
          theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom"))
}

# =============================================================================
# STEP 4: Hyperparameter tuning (cached)
# =============================================================================
grid_rds <- file.path(result_dir, "grid_rf_tuning_AGB.rds")
tune_result <- if (file.exists(grid_rds)) {
  cat("Loading cached tuning...\n"); readRDS(grid_rds)
} else {
  cat("Running grid tuning...\n")
  tune_df <- training_df_final %>% slice_sample(prop = 0.3)
  res <- evaluate_rf_grid_ranger_regression(
    data = tune_df, covars = agb_predictors, response = response_var,
    ntree_values = c(300,500), mtry_values = c(5,8,10),
    nodesize_values = c(3,7), sample_fraction_values = c(0.7,0.8),
    n_iter = 3, train_frac = 0.7, seed = 42)
  saveRDS(res, grid_rds); res
}
write_csv(tune_result$grid_results, file.path(result_dir, "grid_rf_tuning_AGB.csv"))
write_csv(tune_result$best_params,  file.path(result_dir, "grid_rf_best_AGB.csv"))
best         <- tune_result$best_params
model_params <- list(ntree = best$ntree[1], mtry = best$mtry[1],
                     nodesize = best$nodesize[1], sample_fraction = best$sample_fraction[1])
cat("Best hyperparameters:\n"); print(model_params)

# =============================================================================
# STEP 5: Train final model
# =============================================================================
final_rf <- randomForest(
  formula = build_rf_formula(agb_predictors, response_var),
  data    = training_df_final[, c(response_var, agb_predictors)],
  ntree = model_params$ntree, mtry = model_params$mtry, nodesize = model_params$nodesize,
  sampsize = floor(model_params$sample_fraction * nrow(training_df_final)),
  replace = FALSE, importance = TRUE)
saveRDS(final_rf, file.path(result_dir, "final_rf_model_AGB.rds"))
cat("Saved: final_rf_model_AGB.rds\n")

# =============================================================================
# STEP 6: Predict + internal evaluation
# =============================================================================
training_df_final$tAGB_pred <- predict(final_rf, newdata = training_df_final)
cat("\n=== AGB Internal Evaluation (Training) ===\n")
cat(sprintf("RMSE: %.3f\n", caret::RMSE(training_df_final$tAGB_pred, training_df_final$AGB_pred)))
cat(sprintf("MAE : %.3f\n", caret::MAE(training_df_final$tAGB_pred,  training_df_final$AGB_pred)))
cat(sprintf("R2  : %.3f\n", caret::R2(training_df_final$tAGB_pred,   training_df_final$AGB_pred)))

# =============================================================================
# STEP 7: Variable importance + predictor structure
# =============================================================================
var_imp_df <- as.data.frame(final_rf$importance) %>%
  tibble::rownames_to_column("variable") %>%
  select(variable, importance = IncNodePurity) %>% arrange(desc(importance))
write_csv(var_imp_df, file.path(result_dir, "var_importance_rf_AGB.csv"))
predictor_structure <- list(
  classes = sapply(training_df_final[, agb_predictors], class),
  levels  = lapply(training_df_final[, agb_predictors], function(x) if(is.factor(x)) levels(x) else NULL))
saveRDS(predictor_structure, file.path(result_dir, "rf_predictor_structure_AGB.rds"))

# =============================================================================
# STEP 8: Monte Carlo (on bootstrap data)
# =============================================================================
training_df_boot <- training_df_boot %>%
  filter(if_all(all_of(c(response_var, agb_predictors)), ~ !is.na(.)))
cat(sprintf("\nBootstrap data for MC: %d rows\n", nrow(training_df_boot)))
cat("Running Monte Carlo (n=10, bootstrap)...\n")
mc_result  <- run_montecarlo_rf_regression(
  data = training_df_boot, covars = agb_predictors, response = response_var,
  n_iter = n_boot, ntree = model_params$ntree, mtry = model_params$mtry,
  nodesize = model_params$nodesize, sample_fraction = model_params$sample_fraction,
  train_frac = 0.7, seed = 42)
write_csv(mc_result, file.path(result_dir, "rf_agb_mc_results.csv"))
summary_ci <- summarize_ci_regression(mc_result)
write_csv(summary_ci, file.path(result_dir, "rf_agb_mc_summary_ci.csv"))
cat("=== AGB Monte Carlo CI ===\n"); print(summary_ci)

# =============================================================================
# STEP 9: Supplementary plots
# =============================================================================
print(ggplot(tune_result$grid_results,
             aes(x = factor(mtry), y = factor(nodesize), fill = rmse_median)) +
        geom_tile(colour = "white") + geom_text(aes(label = round(rmse_median, 2)), size = 2.5, colour = "white") +
        scale_fill_gradient(low = "steelblue", high = "red") +
        facet_wrap(~ ntree, labeller = label_both) +
        labs(title = "AGB Grid Tuning Heatmap", x = "mtry", y = "nodesize") +
        theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold")))

print(var_imp_df %>% slice_head(n = 20) %>%
        ggplot(aes(x = reorder(variable, importance), y = importance)) +
        geom_col(fill = "grey40") + coord_flip() +
        labs(title = "Top 20 Variable Importance -- AGB", x = NULL, y = "IncNodePurity") +
        theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

print(ggplot(training_df_final, aes(x = AGB_pred, y = tAGB_pred)) +
        geom_point(alpha = 0.3, size = 1) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
        labs(title = "AGB: Observed vs Predicted (Training)",
             x = "Observed AGB (gDW/m2)", y = "Predicted AGB (gDW/m2)") +
        theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

print(training_df_final %>% mutate(resid = tAGB_pred - AGB_pred) %>%
        ggplot(aes(x = resid)) + geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
        labs(title = "AGB Residual Distribution", x = "Residual", y = "Count") +
        theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

print(mc_result %>%
        pivot_longer(c(RMSE, MAE, R2), names_to = "metric", values_to = "value") %>%
        ggplot(aes(x = value)) + geom_histogram(bins = 10, fill = "steelblue", alpha = 0.7) +
        facet_wrap(~ metric, scales = "free") +
        labs(title = "AGB Monte Carlo Distributions", x = "Value", y = "Count") +
        theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold")))

cat("\nE1_rf_model_AGB.R COMPLETE\n")
