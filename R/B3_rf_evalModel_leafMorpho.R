# =============================================================================
# B3_rf_evalModel_leafMorpho.R
# Evaluate MORPH3 predictions per year (supplementary).
#
# Inputs:
#   predicted_MORPH3_probs_YYYY.csv  (from B2_rf_applyModel_leafMorpho.R)
#   master_raw_YYYY.csv              (ground truth morph3)
#
# Outputs:
#   morph3_performance_by_year.csv
#   morph3_confusion_by_year.csv
#   morph3_classwise_metrics_by_year.csv
#   Plots displayed in panel (supplementary paper)
#
# Data availability: this script's inputs are not included in this
# repository.
# =============================================================================

library(dplyr); library(readr); library(purrr)
library(caret);  library(tibble)
library(ggplot2); library(tidyr)
library(here)   # install.packages("here") if you don't have it yet


# Helper: normalise column types that vary across CSV files (prevents bind_rows type conflicts)
normalise_master_types <- function(df) {
  char_cols <- c("OBJECTID","compositio","location","sg_morpho","morph3",
                 grep("^[A-Za-z]{2}_SPC$", names(df), value = TRUE))
  for (col in intersect(char_cols, names(df))) df[[col]] <- as.character(df[[col]])
  if ("year_gt" %in% names(df)) df$year_gt <- suppressWarnings(as.integer(df$year_gt))
  num_cols <- c("PA","tSPC","AGB_pred","AGB_low","AGB_up","AGB_CIwidt",
                "AGC_pred","AGC_low","AGC_up","AGC_CIwidt",
                "carbon_index","xcoord","ycoord")
  for (col in intersect(num_cols, names(df))) df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  df
}


result_dir <- here("result")
data_years <- as.character(2017:2024)

# Helper: multi-class metrics
evaluate_metrics_multiclass <- function(actual, predicted) {
  actual    <- factor(actual)
  predicted <- factor(predicted, levels = levels(actual))
  cm        <- caret::confusionMatrix(predicted, actual)
  byclass   <- cm$byClass
  f1_mean   <- if (is.null(dim(byclass))) as.numeric(byclass["F1"]) else mean(byclass[,"F1"], na.rm = TRUE)
  data.frame(accuracy = round(as.numeric(cm$overall["Accuracy"]), 3),
             kappa    = round(as.numeric(cm$overall["Kappa"]),    3),
             macro_f1 = round(f1_mean, 3))
}

# =============================================================================
# BLOCK 1: Per-year evaluation
# =============================================================================
cat("Evaluating MORPH3 predictions per year...\n")

confusion_list <- list()
classwise_list <- list()

eval_results <- purrr::map_dfr(data_years, function(yr) {
  pred_path   <- file.path(result_dir, paste0("predicted_MORPH3_probs_", yr, ".csv"))
  master_path <- file.path(result_dir, paste0("master_raw_",             yr, ".csv"))
  if (!file.exists(pred_path) || !file.exists(master_path)) {
    message(sprintf("Year %s: missing file.", yr)); return(NULL)
  }
  df_pred  <- read_csv(pred_path,   show_col_types = FALSE)
  df_truth <- normalise_master_types(read_csv(master_path, show_col_types = FALSE)) %>% select(gee_id, morph3)
  
  df_eval <- inner_join(df_truth, df_pred, by = "gee_id") %>% filter(!is.na(morph3))
  if (nrow(df_eval) == 0) { message(sprintf("Year %s: 0 valid rows.", yr)); return(NULL) }
  
  # Predicted class = highest P_* probability
  prob_cols           <- grep("^P_", names(df_eval), value = TRUE)
  df_eval$morph3_pred <- sub("^P_", "", colnames(df_eval[, prob_cols])[max.col(df_eval[, prob_cols])])
  
  actual    <- df_eval$morph3
  predicted <- df_eval$morph3_pred
  
  if (length(unique(actual)) < 2) {
    warning(sprintf("Year %s: only one class -- skipped.", yr))
    return(data.frame(accuracy = NA, kappa = NA, macro_f1 = NA, year = yr, n = nrow(df_eval)))
  }
  
  all_levels <- sort(unique(c(as.character(actual), as.character(predicted))))
  cm <- caret::confusionMatrix(factor(predicted, levels = all_levels), factor(actual, levels = all_levels))
  
  # Store confusion matrix
  cm_df <- as.data.frame(cm$table)
  names(cm_df) <- c("Predicted","Reference","Count")
  confusion_list[[yr]] <<- cm_df %>% mutate(year = yr, n_test = nrow(df_eval))
  
  # Store class-wise metrics
  byclass_raw <- cm$byClass
  byclass_df  <- if (is.null(dim(byclass_raw)))
    as.data.frame(t(byclass_raw)) %>% mutate(Class = "overall")
  else
    as.data.frame(byclass_raw) %>% tibble::rownames_to_column("Class")
  classwise_list[[yr]] <<- byclass_df %>%
    select(any_of(c("Class","Precision","Recall","F1"))) %>%
    mutate(year = yr, n_test = nrow(df_eval))
  
  met      <- evaluate_metrics_multiclass(actual, predicted)
  met$year <- yr; met$n <- nrow(df_eval)
  cat(sprintf("Year %s: Accuracy=%.3f | Kappa=%.3f | Macro-F1=%.3f | n=%d\n",
              yr, met$accuracy, met$kappa, met$macro_f1, met$n))
  met
})

write_csv(eval_results,   file.path(result_dir, "morph3_performance_by_year.csv"))
cat("\n=== MORPH3 Model Performance per Year ===\n"); print(eval_results)

confusion_df <- bind_rows(confusion_list)
write_csv(confusion_df, file.path(result_dir, "morph3_confusion_by_year.csv"))

classwise_df <- bind_rows(classwise_list)
write_csv(classwise_df, file.path(result_dir, "morph3_classwise_metrics_by_year.csv"))
cat("Saved: morph3_performance, confusion, classwise CSVs\n")

# =============================================================================
# BLOCK 2: Supplementary plots
# =============================================================================

# P1: Performance metrics per year
print(eval_results %>%
        filter(!is.na(accuracy)) %>%
        select(year, accuracy, kappa, macro_f1) %>%
        pivot_longer(c(accuracy, kappa, macro_f1), names_to = "metric", values_to = "value") %>%
        mutate(year = factor(year, levels = data_years)) %>%
        ggplot(aes(x = year, y = value, fill = metric)) +
        geom_col(position = "dodge") +
        geom_text(aes(label = round(value, 3)), position = position_dodge(0.9), vjust = -0.4, size = 3) +
        labs(title = "MORPH3 Model Performance per Year", x = "Year", y = "Metric Value") +
        theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

# P2: Confusion matrix heatmap per year
print(ggplot(confusion_df, aes(x = Reference, y = Predicted, fill = Count)) +
        geom_tile(colour = "white") + geom_text(aes(label = Count), size = 3) +
        scale_fill_gradient(low = "white", high = "#2980b9") +
        facet_wrap(~ year) +
        labs(title = "Confusion Matrix per Year -- MORPH3", x = "Actual", y = "Predicted") +
        theme_minimal(base_size = 11) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1),
              plot.title  = element_text(face = "bold"), strip.text = element_text(face = "bold")))

# P3: Class-wise F1 per year
f1_df <- classwise_df %>% filter(!is.na(F1)) %>% mutate(year = factor(year, levels = data_years))
if (nrow(f1_df) > 0) {
  print(ggplot(f1_df, aes(x = year, y = F1, fill = Class)) +
          geom_col(position = "dodge") +
          geom_text(aes(label = round(F1, 2)), position = position_dodge(0.9), vjust = -0.4, size = 3) +
          labs(title = "Class-wise F1 Score per Year -- MORPH3", x = "Year", y = "F1 Score") +
          theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))
}

# P4: Probability boxplots per year (one plot per year)
for (yr in data_years) {
  pred_path   <- file.path(result_dir, paste0("predicted_MORPH3_probs_", yr, ".csv"))
  master_path <- file.path(result_dir, paste0("master_raw_",             yr, ".csv"))
  if (!file.exists(pred_path) || !file.exists(master_path)) next
  
  df_master <- normalise_master_types(read_csv(master_path, show_col_types = FALSE))
  if ("OBJECTID" %in% names(df_master)) df_master$OBJECTID <- as.character(df_master$OBJECTID)
  
  df_box <- inner_join(
    df_master %>% select(gee_id, morph3),
    read_csv(pred_path, show_col_types = FALSE),
    by = "gee_id"
  ) %>% filter(!is.na(morph3))
  
  prob_cols <- grep("^P_", names(df_box), value = TRUE)
  if (length(prob_cols) == 0) next
  
  # Build long format and filter out rows with NA probability or empty class label
  df_long <- df_box %>%
    select(morph3, all_of(prob_cols)) %>%
    pivot_longer(all_of(prob_cols), names_to = "class", values_to = "probability") %>%
    mutate(class = sub("^P_", "", class)) %>%
    filter(!is.na(probability), nchar(class) > 0, !is.na(morph3))
  
  # Skip year if no valid data after filtering
  if (nrow(df_long) == 0 || length(unique(df_long$class)) == 0) {
    message(sprintf("Year %s: no valid probability data for boxplot -- skipped.", yr))
    next
  }
  
  print(ggplot(df_long, aes(x = morph3, y = probability, fill = morph3)) +
          geom_boxplot(outlier.size = 0.3, alpha = 0.7, na.rm = TRUE) +
          facet_wrap(~ class) +
          labs(title = sprintf("Morphology Probability Boxplot -- Year %s", yr),
               x = "True Class", y = "Predicted Probability") +
          theme_minimal(base_size = 12) +
          theme(legend.position = "none", plot.title = element_text(face = "bold"),
                axis.text.x = element_text(angle = 30, hjust = 1)))
}

cat("\nB3_rf_evalModel_leafMorpho.R COMPLETE\n")
