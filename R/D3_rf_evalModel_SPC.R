# =============================================================================
# D3_rf_evalModel_SPC.R
# Evaluate SPC predictions per year (supplementary).
#
# Inputs:
#   predicted_tSPC_YYYY.csv  (from D2_rf_applyModel_SPC.R)
#   master_raw_YYYY.csv      (ground truth tSPC)
#
# Outputs:
#   spc_performance_by_year.csv
#   Plots displayed in panel (supplementary paper)
#
# Data availability: this script's inputs are not included in this
# repository.
# =============================================================================

library(dplyr); library(readr); library(purrr)
library(ggplot2); library(tidyr); library(maps)
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

rmse_fn <- function(a, p) sqrt(mean((a - p)^2, na.rm = TRUE))
mae_fn  <- function(a, p) mean(abs(a - p), na.rm = TRUE)
r2_fn   <- function(a, p) {
  ss_res <- sum((a - p)^2, na.rm = TRUE)
  ss_tot <- sum((a - mean(a, na.rm = TRUE))^2, na.rm = TRUE)
  if (ss_tot == 0) return(NA_real_); 1 - ss_res / ss_tot
}

# =============================================================================
# BLOCK 1: Per-year metrics
# =============================================================================
cat("Evaluating SPC predictions per year...\n")

results_yearly <- purrr::map_dfr(data_years, function(yr) {
  pred_path   <- file.path(result_dir, paste0("predicted_tSPC_", yr, ".csv"))
  master_path <- file.path(result_dir, paste0("master_raw_",     yr, ".csv"))
  if (!file.exists(pred_path) || !file.exists(master_path)) {
    message(sprintf("Year %s: missing file.", yr)); return(NULL)
  }
  df_pred <- read_csv(pred_path,   show_col_types = FALSE, guess_max = 100000) %>%
    select(gee_id, tSPC_pred)
  df_ref  <- normalise_master_types(read_csv(master_path, show_col_types = FALSE)) %>%
    select(gee_id, tSPC)
  
  df <- inner_join(df_ref, df_pred, by = "gee_id") %>%
    filter(!is.na(tSPC), !is.na(tSPC_pred), tSPC > 0, tSPC <= 100)
  if (nrow(df) == 0) { message(sprintf("Year %s: no valid rows.", yr)); return(NULL) }
  cat(sprintf("Year %s: %d matched samples.\n", yr, nrow(df)))
  
  tibble(year = yr, n = nrow(df),
         RMSE = round(rmse_fn(df$tSPC, df$tSPC_pred), 3),
         MAE  = round(mae_fn(df$tSPC,  df$tSPC_pred), 3),
         R2   = round(r2_fn(df$tSPC,   df$tSPC_pred), 3))
})

write_csv(results_yearly, file.path(result_dir, "spc_performance_by_year.csv"))
cat("\n=== SPC Model Performance per Year ===\n"); print(results_yearly)

# =============================================================================
# BLOCK 2: Load all years combined
# =============================================================================
data_all <- purrr::map_dfr(data_years, function(yr) {
  pred_path   <- file.path(result_dir, paste0("predicted_tSPC_", yr, ".csv"))
  master_path <- file.path(result_dir, paste0("master_raw_",     yr, ".csv"))
  if (!file.exists(pred_path) || !file.exists(master_path)) return(NULL)
  df_pred <- read_csv(pred_path,   show_col_types = FALSE, guess_max = 100000) %>% select(gee_id, tSPC_pred, any_of(c("xcoord","ycoord")))
  df_ref  <- normalise_master_types(read_csv(master_path, show_col_types = FALSE)) %>% select(gee_id, tSPC, any_of(c("xcoord","ycoord")))
  inner_join(df_ref, df_pred, by = "gee_id", suffix = c("",".p")) %>%
    filter(!is.na(tSPC), !is.na(tSPC_pred), tSPC > 0, tSPC <= 100) %>%
    mutate(year = yr)
})
cat(sprintf("\nTotal rows loaded: %d\n\n", nrow(data_all)))

# =============================================================================
# BLOCK 3: Supplementary plots
# =============================================================================

# P1: RMSE & MAE per year
print(results_yearly %>%
        select(year, RMSE, MAE) %>%
        pivot_longer(c(RMSE, MAE), names_to = "metric", values_to = "value") %>%
        mutate(year = factor(year, levels = data_years)) %>%
        ggplot(aes(x = year, y = value, fill = metric)) +
        geom_col(position = "dodge") +
        geom_text(aes(label = round(value, 2)), position = position_dodge(0.9), vjust = -0.4, size = 3) +
        scale_fill_manual(values = c("RMSE" = "grey30", "MAE" = "grey70")) +
        labs(title = "SPC Model: RMSE and MAE per Year", x = "Year", y = "Error (%)") +
        theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

# P2: R2 per year
print(ggplot(results_yearly, aes(x = factor(year, levels = data_years), y = R2)) +
        geom_col(fill = "grey50") + geom_text(aes(label = round(R2, 3)), vjust = -0.4, size = 3.5) +
        labs(title = "SPC Model R2 per Year", x = "Year", y = "R2") +
        theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

# P3: Scatter obs vs pred all years
print(ggplot(data_all, aes(x = tSPC, y = tSPC_pred, colour = factor(year, levels = data_years))) +
        geom_point(alpha = 0.3, size = 1) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        labs(title = "SPC: Observed vs Predicted (All Years)",
             x = "Observed tSPC (%)", y = "Predicted tSPC (%)", colour = "Year") +
        theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

# P4: Scatter faceted by year
print(ggplot(data_all, aes(x = tSPC, y = tSPC_pred)) +
        geom_point(alpha = 0.3, size = 1) + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        facet_wrap(~ year) +
        labs(title = "SPC: Observed vs Predicted per Year",
             x = "Observed tSPC (%)", y = "Predicted tSPC (%)") +
        theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold")))

# P5: Boxplot by predicted bin
print(data_all %>%
        mutate(Pred_Bin = cut(tSPC_pred, breaks = seq(0, 100, 10), include.lowest = TRUE)) %>%
        ggplot(aes(x = Pred_Bin, y = tSPC)) + geom_boxplot(fill = "grey85") +
        labs(title = "Actual tSPC by Predicted Interval", x = "Predicted tSPC Bin (%)", y = "Actual tSPC (%)") +
        theme_minimal(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold")))

# P6: Absolute error by actual bin
print(data_all %>%
        mutate(Error = abs(tSPC - tSPC_pred),
               Actual_Bin = cut(tSPC, breaks = seq(0, 100, 10), include.lowest = TRUE)) %>%
        ggplot(aes(x = Actual_Bin, y = Error)) + geom_boxplot(fill = "grey85") +
        labs(title = "Absolute Error by Actual tSPC Interval", x = "Actual tSPC Bin (%)", y = "Absolute Error") +
        theme_minimal(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold")))

# P7: Residual distribution
print(data_all %>% mutate(resid = tSPC_pred - tSPC) %>%
        ggplot(aes(x = resid)) +
        geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7, colour = "white") +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
        labs(title = "SPC Residual Distribution (Pred - Obs)", x = "Residual", y = "Count") +
        theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

# P8: Spatial map
x_col <- if ("xcoord" %in% names(data_all)) "xcoord" else "lon"
y_col <- if ("ycoord" %in% names(data_all)) "ycoord" else "lat"
if (all(c(x_col, y_col) %in% names(data_all))) {
  df_map <- data_all %>%
    mutate(x = suppressWarnings(as.numeric(.data[[x_col]])),
           y = suppressWarnings(as.numeric(.data[[y_col]]))) %>% filter(is.finite(x), is.finite(y))
  print(ggplot() +
          geom_polygon(data = ggplot2::map_data("world", region = "Indonesia"),
                       aes(x = long, y = lat, group = group), fill = "#d4e6c3", colour = "grey50", linewidth = 0.3) +
          geom_point(data = df_map, aes(x = x, y = y, colour = tSPC_pred), size = 1.2, alpha = 0.7) +
          scale_colour_gradient(name = "Predicted SPC (%)", low = "#ffffb2", high = "#006837") +
          coord_quickmap(xlim = range(df_map$x) + c(-1,1), ylim = range(df_map$y) + c(-1,1)) +
          labs(title = "Spatial Distribution - Predicted SPC (All Years)", x = "Longitude", y = "Latitude") +
          theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom"))
}

cat("\nD3_rf_evalModel_SPC.R COMPLETE\n")
