# =============================================================================
# C2_rf_applyModel_cIndex.R
# Apply carbon index model to the full GEE grid.
#
# Input:  predicted_PA_YYYY.csv (from A2_rf_applyModel_PA.R)
# Output: predicted_CINDEX_YYYY.csv (carbon_index_pred; all rows retained)
#
# Data availability: this script's inputs are not included in this
# repository.
# =============================================================================

library(dplyr); library(readr); library(randomForest)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "func_reg_rf.R"))

result_dir  <- here("result")
data_years  <- as.character(2017:2024)

model               <- readRDS(file.path(result_dir, "final_rf_model_CINDEX.rds"))
predictor_structure <- readRDS(file.path(result_dir, "rf_predictor_structure_CINDEX.rds"))
cindex_predictors   <- names(predictor_structure$classes)
cat(sprintf("CINDEX model loaded | Predictors: %d\n\n", length(cindex_predictors)))

for (yr in data_years) {
  cat(sprintf("\nPredicting carbon_index for year: %s\n", yr))
  in_path  <- file.path(result_dir, paste0("predicted_PA_",     yr, ".csv"))
  out_path <- file.path(result_dir, paste0("predicted_CINDEX_", yr, ".csv"))
  if (!file.exists(in_path)) { warning(sprintf("File not found: %s", in_path)); next }
  
  df <- read_csv(in_path, show_col_types = FALSE, guess_max = 100000)
  for (mc in setdiff(cindex_predictors, names(df))) df[[mc]] <- NA_real_
  
  valid_idx             <- df %>% select(all_of(cindex_predictors)) %>% complete.cases()
  df$carbon_index_pred  <- NA_real_
  cat(sprintf("Year %s: total=%d | valid=%d | NA-rows=%d\n",
              yr, nrow(df), sum(valid_idx), sum(!valid_idx)))
  
  if (sum(valid_idx) > 0) {
    preds <- apply_model_regression(df[valid_idx,], cindex_predictors, model, predictor_structure)
    df$carbon_index_pred[valid_idx] <- pmin(pmax(preds, 0), 1)
    n_below <- sum(preds < 0); n_above <- sum(preds > 1)
    if (n_below + n_above > 0)
      cat(sprintf("Clamped: %d below 0 | %d above 1\n", n_below, n_above))
    cat(sprintf("carbon_index_pred: min=%.4f max=%.4f mean=%.4f\n",
                min(df$carbon_index_pred, na.rm=T), max(df$carbon_index_pred, na.rm=T),
                mean(df$carbon_index_pred, na.rm=T)))
  }
  write_csv(df, out_path)
  cat(sprintf("Saved: predicted_CINDEX_%s.csv (%d rows)\n", yr, nrow(df)))
}

cat("\nC2_rf_applyModel_cIndex.R COMPLETE\n")
