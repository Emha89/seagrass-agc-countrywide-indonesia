# =============================================================================
# E2_rf_applyModel_AGB.R
# Apply AGB model to the full GEE grid.
#
# Input:  predicted_tSPC_YYYY.csv (from D2_rf_applyModel_SPC.R)
# Output: predicted_tAGB_YYYY.csv (all rows retained; NA if predictor incomplete)
#
# NOTE: reads from predicted_tSPC (not predicted_PA directly), continuing
# the chain of accumulated columns from each upstream apply script, even
# though AGB's own predictors (depth + GSE bands only) don't actually
# need tSPC_pred -- this keeps every downstream stage able to read one
# file with everything accumulated so far, rather than re-joining from
# predicted_PA each time.
#
# Data availability: this script's inputs are not included in this
# repository.
# =============================================================================

library(dplyr); library(readr); library(randomForest)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "func_reg_rf.R"))

result_dir  <- here("result")
data_years  <- as.character(2017:2024)

model               <- readRDS(file.path(result_dir, "final_rf_model_AGB.rds"))
predictor_structure <- readRDS(file.path(result_dir, "rf_predictor_structure_AGB.rds"))
agb_predictors      <- names(predictor_structure$classes)
cat(sprintf("AGB model loaded | Predictors: %d\n\n", length(agb_predictors)))

for (yr in data_years) {
  cat(sprintf("\nPredicting tAGB for year: %s\n", yr))
  in_path  <- file.path(result_dir, paste0("predicted_tSPC_", yr, ".csv"))
  out_path <- file.path(result_dir, paste0("predicted_tAGB_", yr, ".csv"))
  if (!file.exists(in_path)) { warning(sprintf("File not found: %s", in_path)); next }
  
  df <- read_csv(in_path, show_col_types = FALSE, guess_max = 100000)
  for (mc in setdiff(agb_predictors, names(df))) df[[mc]] <- NA_real_
  
  valid_idx    <- df %>% select(all_of(agb_predictors)) %>% complete.cases()
  df$tAGB_pred <- NA_real_
  cat(sprintf("Year %s: total=%d | valid=%d | NA-rows=%d\n",
              yr, nrow(df), sum(valid_idx), sum(!valid_idx)))
  
  if (sum(valid_idx) > 0) {
    preds <- apply_model_regression(df[valid_idx,], agb_predictors, model, predictor_structure)
    df$tAGB_pred[valid_idx] <- pmax(preds, 0)
    cat(sprintf("tAGB_pred: min=%.3f max=%.3f mean=%.3f\n",
                min(df$tAGB_pred, na.rm=T), max(df$tAGB_pred, na.rm=T), mean(df$tAGB_pred, na.rm=T)))
  }
  write_csv(df, out_path)
  cat(sprintf("Saved: predicted_tAGB_%s.csv (%d rows)\n", yr, nrow(df)))
}

cat("\nE2_rf_applyModel_AGB.R COMPLETE\n")
