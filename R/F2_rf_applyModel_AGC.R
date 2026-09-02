# =============================================================================
# F2_rf_applyModel_AGC.R
# Apply AGC model to the full GEE grid -- final stage of the R pipeline.
#
# Input:  predicted_tAGB_YYYY.csv (from E2_rf_applyModel_AGB.R -- already
#         carries tSPC_pred forward from the accumulated-columns chain)
#         + predicted_CINDEX_YYYY.csv (from C2_rf_applyModel_cIndex.R --
#           carbon_index_pred joined in as "carbon_index" to match the
#           predictor name the model expects)
#         + predicted_MORPH3_probs_YYYY.csv (joined for completeness, but
#           NOT actually used as a predictor -- see note below)
# Output: predicted_tAGC_YYYY.csv (all rows retained; NA if predictor incomplete)
#
# NOTE: agc_predictors is read dynamically from the saved
# rf_predictor_structure_AGC.rds (F1_rf_model_AGC.R's confirmed leaner
# set: 64 GSE bands + tAGB_pred + carbon_index + tSPC_pred -- no depth,
# no PA_prob, no morphology probabilities). The carbon_index join below
# is functionally necessary (it's an active predictor). The P_morph3
# join is vestigial here -- it runs and adds the P_* columns to df, but
# since avail_preds filters down to only what's actually in
# agc_predictors, those joined columns never reach the prediction call.
# Consistent with the same leftover pattern noted in F1.
#
# Data availability: this script's inputs are not included in this
# repository.
# =============================================================================

library(dplyr); library(readr); library(randomForest)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "func_reg_rf.R"))

result_dir  <- here("result")
data_years  <- as.character(2017:2024)
morph_cols  <- c("P_mixed_short_plus_mono_short","P_mixed_long","P_mono_Ea")

model               <- readRDS(file.path(result_dir, "final_rf_model_AGC.rds"))
predictor_structure <- readRDS(file.path(result_dir, "rf_predictor_structure_AGC.rds"))
agc_predictors      <- names(predictor_structure$classes)
cat(sprintf("AGC model loaded | Predictors: %d\n\n", length(agc_predictors)))

for (yr in data_years) {
  cat(sprintf("\nPredicting tAGC for year: %s\n", yr))
  in_path     <- file.path(result_dir, paste0("predicted_tAGB_",        yr, ".csv"))
  cindex_path <- file.path(result_dir, paste0("predicted_CINDEX_",      yr, ".csv"))
  morph_path  <- file.path(result_dir, paste0("predicted_MORPH3_probs_", yr, ".csv"))
  out_path    <- file.path(result_dir, paste0("predicted_tAGC_",         yr, ".csv"))
  if (!file.exists(in_path)) { warning(sprintf("File not found: %s", in_path)); next }
  df <- read_csv(in_path, show_col_types = FALSE, guess_max = 100000)

  # Join carbon_index_pred if missing
  if (!"carbon_index" %in% names(df) && file.exists(cindex_path)) {
    df_ci <- read_csv(cindex_path, show_col_types = FALSE) %>%
      select(gee_id, carbon_index = carbon_index_pred)
    df    <- left_join(df, df_ci, by = "gee_id")
    cat(sprintf("carbon_index joined for year %s\n", yr))
  }

  # Join P_morph3 if missing
  missing_morph <- setdiff(morph_cols, names(df))
  if (length(missing_morph) > 0 && file.exists(morph_path)) {
    df_morph <- read_csv(morph_path, show_col_types = FALSE) %>%
      select(gee_id, any_of(morph_cols))
    df <- left_join(df, df_morph, by = "gee_id")
    cat(sprintf("P_morph3 joined for year %s\n", yr))
  }

  # Resolve duplicate P_* from successive joins
  if (any(grepl("^P_.*\\.y$", names(df)))) {
    df <- df %>%
      select(-matches("^P_.*\\.y$")) %>%
      rename_with(~ sub("\\.x$", "", .), matches("^P_.*\\.x$"))
  }

  for (mc in setdiff(agc_predictors, names(df))) df[[mc]] <- NA_real_
  avail_preds <- intersect(agc_predictors, names(df))
  valid_idx    <- df %>% select(all_of(avail_preds)) %>% complete.cases()
  df$tAGC_pred <- NA_real_
  cat(sprintf("Year %s: total=%d | valid=%d | NA-rows=%d\n",
              yr, nrow(df), sum(valid_idx), sum(!valid_idx)))

  if (sum(valid_idx) > 0) {
    preds <- apply_model_regression(df[valid_idx,], avail_preds, model, predictor_structure)
    df$tAGC_pred[valid_idx] <- pmax(preds, 0)
    cat(sprintf("tAGC_pred: min=%.3f max=%.3f mean=%.3f\n",
                min(df$tAGC_pred, na.rm=T), max(df$tAGC_pred, na.rm=T), mean(df$tAGC_pred, na.rm=T)))
  }

  write_csv(df, out_path)
  cat(sprintf("Saved: predicted_tAGC_%s.csv (%d rows)\n", yr, nrow(df)))
}

cat("\nF2_rf_applyModel_AGC.R COMPLETE\n")
