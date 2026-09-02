# =============================================================================
# D2_rf_applyModel_SPC.R
# Apply SPC model to the full GEE grid.
#
# Input:  predicted_PA_YYYY.csv (from A2_rf_applyModel_PA.R)
# Output: predicted_tSPC_YYYY.csv
#
# NOTE: unlike D1_rf_model_SPC.R's training step (which only checks
# whether predicted_PA already has the morph3 P_* columns), this apply
# script joins them itself directly from predicted_MORPH3_probs_<year>.csv
# (B2_rf_applyModel_leafMorpho.R's output) whenever they're missing --
# so it doesn't strictly depend on predicted_PA already having morph3
# merged in. Still needs B2 to have run at some point to produce that
# file in the first place.
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

model               <- readRDS(file.path(result_dir, "final_rf_model_PCT.rds"))
predictor_structure <- readRDS(file.path(result_dir, "rf_predictor_structure_PCT.rds"))
pct_predictors      <- names(predictor_structure$classes)
cat(sprintf("SPC model loaded | Predictors: %d\n\n", length(pct_predictors)))

for (yr in data_years) {
  cat(sprintf("\nPredicting tSPC for year: %s\n", yr))
  in_path  <- file.path(result_dir, paste0("predicted_PA_",   yr, ".csv"))
  out_path <- file.path(result_dir, paste0("predicted_tSPC_", yr, ".csv"))
  if (!file.exists(in_path)) { warning(sprintf("File not found: %s", in_path)); next }
  df <- read_csv(in_path, show_col_types = FALSE, guess_max = 100000)

  # Join P_morph3 if missing
  missing_morph <- setdiff(morph_cols, names(df))
  if (length(missing_morph) > 0) {
    morph_path <- file.path(result_dir, paste0("predicted_MORPH3_probs_", yr, ".csv"))
    if (file.exists(morph_path)) {
      df_morph <- read_csv(morph_path, show_col_types = FALSE) %>%
        select(gee_id, any_of(morph_cols))
      df <- left_join(df, df_morph, by = "gee_id")
    } else {
      for (mc in missing_morph) df[[mc]] <- NA_real_
    }
  }

  for (mc in setdiff(pct_predictors, names(df))) df[[mc]] <- NA_real_
  valid_idx    <- df %>% select(all_of(pct_predictors)) %>% complete.cases()
  df$tSPC_pred <- NA_real_
  cat(sprintf("Year %s: total=%d | valid=%d | NA-rows=%d\n",
              yr, nrow(df), sum(valid_idx), sum(!valid_idx)))

  if (sum(valid_idx) > 0) {
    preds <- apply_model_regression(df[valid_idx,], pct_predictors, model, predictor_structure)
    df$tSPC_pred[valid_idx] <- pmin(pmax(preds, 0), 100)
    cat(sprintf("tSPC_pred: min=%.2f max=%.2f mean=%.2f\n",
                min(df$tSPC_pred, na.rm=T), max(df$tSPC_pred, na.rm=T), mean(df$tSPC_pred, na.rm=T)))
  }

  write_csv(df, out_path)
  cat(sprintf("Saved: predicted_tSPC_%s.csv (%d rows)\n", yr, nrow(df)))
}

cat("\nD2_rf_applyModel_SPC.R COMPLETE\n")
