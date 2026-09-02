# =============================================================================
# B2_rf_applyModel_leafMorpho.R
# Apply MORPH3 model to the full GEE grid.
#
# Input:  predicted_PA_YYYY.csv (from A2_rf_applyModel_PA.R -- contains
#         GSE + env predictors, plus PA_prob/PA_pred)
# Output: predicted_MORPH3_probs_YYYY.csv (P_* probability columns)
#
# EXECUTION ORDER / MULTI-PASS NOTE: this script depends on
# A2_rf_applyModel_PA.R having already produced predicted_PA_<year>.csv.
# A2 in turn optionally merges in predicted_MORPH3_probs_<year>.csv IF
# that file already exists -- so on a first run, A2 runs without morph3
# data (that merge is simply skipped), this script (B2) runs next to
# produce it, and only a SECOND run of A2 would actually pick the morph3
# probabilities up. Same multi-pass pattern as the PA/morphology
# relationship in the Study 2 repository.
#
# Note: this script sources func_reg_rf.R, but doesn't appear to call any
# function from it (only base predict() is used below) -- likely a
# leftover from a shared template rather than a real dependency. Left
# as-is since it's harmless, just worth knowing.
#
# Data availability: this script's inputs are not included in this
# repository.
# =============================================================================

library(dplyr); library(readr); library(randomForest)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "func_reg_rf.R"))

result_dir  <- here("result")
data_years  <- as.character(2017:2024)

model               <- readRDS(file.path(result_dir, "final_rf_model_MORPH3.rds"))
predictor_structure <- readRDS(file.path(result_dir, "rf_predictor_structure_MORPH3.rds"))
morph_predictors    <- names(predictor_structure$classes)
cat(sprintf("MORPH3 model loaded | Predictors: %d\n\n", length(morph_predictors)))

for (yr in data_years) {
  cat(sprintf("\nPredicting MORPH3 for year: %s\n", yr))
  in_path  <- file.path(result_dir, paste0("predicted_PA_",          yr, ".csv"))
  out_path <- file.path(result_dir, paste0("predicted_MORPH3_probs_", yr, ".csv"))
  if (!file.exists(in_path)) { warning(sprintf("File not found: %s", in_path)); next }

  df <- read_csv(in_path, show_col_types = FALSE, guess_max = 100000)
  for (mc in setdiff(morph_predictors, names(df))) df[[mc]] <- NA_real_

  valid_idx <- df %>% select(all_of(morph_predictors)) %>% complete.cases()
  cat(sprintf("Year %s: total=%d | valid=%d | NA-rows=%d\n",
              yr, nrow(df), sum(valid_idx), sum(!valid_idx)))

  prob_cols <- c("P_mixed_short_plus_mono_short","P_mixed_long","P_mono_Ea")
  for (cc in prob_cols) df[[cc]] <- NA_real_

  if (sum(valid_idx) > 0) {
    prob_mat           <- predict(model, newdata = df[valid_idx, morph_predictors], type = "prob")
    prob_df            <- as.data.frame(prob_mat)
    colnames(prob_df)  <- paste0("P_", colnames(prob_df))
    df[valid_idx, intersect(names(prob_df), names(df))] <- prob_df[, intersect(names(prob_df), names(df))]
  }

  write_csv(df, out_path)
  cat(sprintf("Saved: predicted_MORPH3_probs_%s.csv (%d rows)\n", yr, nrow(df)))
}

cat("\nB2_rf_applyModel_leafMorpho.R COMPLETE\n")
