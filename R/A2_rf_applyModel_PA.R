# =============================================================================
# A2_rf_applyModel_PA.R
# Apply PA model to the full GEE grid (all pixels per year).
#
# Input:  GSE_training_YYYY_CSV.csv (raw GEE image export)
# Output: predicted_PA_YYYY.csv (PA_prob, PA_pred for every pixel)
#
# NOTE ON THE INPUT FILE: this is the same GSE_training_<year>_CSV.csv
# file 01_build_master_raw.R reads, but used differently here. There, it
# gets INNER JOINed against the GT label file, which narrows it down to
# only the labelled training points. Here, it's read directly with no GT
# join, so every pixel in the file gets a prediction -- that's how this
# script covers the "full grid" rather than just training locations.
#
# THRESHOLD: best_threshold = 0.6 must match the threshold used in
# A1_rf_model_PA.R (same manually-fixed value, not the auto-tuned one).
#
# Data availability: this script's inputs are not included in this
# repository.
# =============================================================================

library(dplyr); library(readr); library(randomForest); library(purrr)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "func_reg_rf.R"))

result_dir <- here("result")

data_paths <- list(
  "2017" = here("data", "GSE_training_2017_CSV.csv"),
  "2018" = here("data", "GSE_training_2018_CSV.csv"),
  "2019" = here("data", "GSE_training_2019_CSV.csv"),
  "2020" = here("data", "GSE_training_2020_CSV.csv"),
  "2021" = here("data", "GSE_training_2021_CSV.csv"),
  "2022" = here("data", "GSE_training_2022_CSV.csv"),
  "2023" = here("data", "GSE_training_2023_CSV.csv"),
  "2024" = here("data", "GSE_training_2024_CSV.csv")
)

best_threshold <- 0.6   # must match threshold used in A1_rf_model_PA.R

model               <- readRDS(file.path(result_dir, "final_rf_model_PA.rds"))
predictor_structure <- readRDS(file.path(result_dir, "rf_predictor_structure_PA.rds"))
pa_predictors       <- names(predictor_structure$classes)
cat(sprintf("PA model loaded | Predictors: %d\n\n", length(pa_predictors)))

coerce_types <- function(df, cls_list, lvl_list) {
  for (nm in names(cls_list)) {
    if (!nm %in% names(df)) next
    cls <- cls_list[[nm]][1]
    if      (cls == "numeric")   df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
    else if (cls == "integer")   df[[nm]] <- suppressWarnings(as.integer(df[[nm]]))
    else if (cls == "factor")    df[[nm]] <- factor(df[[nm]], levels = lvl_list[[nm]])
    else if (cls == "character") df[[nm]] <- as.character(df[[nm]])
  }
  df
}

for (yr in names(data_paths)) {
  cat(sprintf("\nPredicting PA for year: %s\n", yr))
  in_path  <- data_paths[[yr]]
  out_path <- file.path(result_dir, paste0("predicted_PA_", yr, ".csv"))
  if (!file.exists(in_path)) { warning(sprintf("File not found: %s", in_path)); next }

  df_raw  <- read_csv(in_path, show_col_types = FALSE, guess_max = 100000) %>%
    filter(!is.na(xcoord), !is.na(ycoord))

  for (mc in setdiff(pa_predictors, names(df_raw))) df_raw[[mc]] <- NA_real_
  df_pred  <- df_raw %>% select(any_of(pa_predictors))
  df_pred  <- coerce_types(df_pred, predictor_structure$classes, predictor_structure$levels)

  valid_idx <- complete.cases(df_pred)
  cat(sprintf("Year %s: total=%d | valid=%d | NA-rows=%d\n",
              yr, nrow(df_raw), sum(valid_idx), sum(!valid_idx)))

  df_raw$PA_prob <- NA_real_; df_raw$PA_pred <- NA_integer_
  if (sum(valid_idx) > 0) {
    probs <- predict(model, newdata = df_pred[valid_idx, ], type = "prob")[, "1"]
    df_raw$PA_prob[valid_idx] <- probs
    df_raw$PA_pred[valid_idx] <- as.integer(probs > best_threshold)
  }

  # Join Morph3 probs if available
  morph_path <- file.path(result_dir, paste0("predicted_MORPH3_probs_", yr, ".csv"))
  if (file.exists(morph_path)) {
    df_morph  <- read_csv(morph_path, show_col_types = FALSE) %>%
      select(gee_id, starts_with("P_"))
    df_raw    <- left_join(df_raw, df_morph, by = "gee_id")
  }

  write_csv(df_raw, out_path)
  cat(sprintf("Saved: predicted_PA_%s.csv (%d rows)\n", yr, nrow(df_raw)))
}

cat("\nA2_rf_applyModel_PA.R COMPLETE\n")
