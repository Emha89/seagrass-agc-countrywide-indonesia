# =============================================================================
# E3_rf_evalModel_AGB.R
# Evaluate AGB predictions per year (supplementary).
#
# Inputs:
#   predicted_tAGB_YYYY.csv  (from E2_rf_applyModel_AGB.R)
#   master_AGB_YYYY.csv      (ground truth AGB_pred from 04_build_master_AGB_AGC.R)
#
# Outputs:
#   agb_performance_by_year.csv
#   Plots displayed in panel (supplementary paper)
#
# Data availability: this script's inputs are not included in this
# repository.
# =============================================================================

library(dplyr); library(readr); library(purrr)
library(ggplot2); library(tidyr); library(maps)
library(here)   # install.packages("here") if you don't have it yet


# Helper: normalise column types that vary across CSV files. Covers all
# columns present in master_raw, master_AGB, and master_AGC.
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
cat("Evaluating AGB predictions per year...\n")

results_yearly <- purrr::map_dfr(data_years, function(yr) {
  pred_path   <- file.path(result_dir, paste0("predicted_tAGB_", yr, ".csv"))
  master_path <- file.path(result_dir, paste0("master_AGB_",     yr, ".csv"))
  if (!file.exists(pred_path) || !file.exists(master_path)) {
    message(sprintf("Year %s: missing file.", yr)); return(NULL)
  }
  df_pred <- read_csv(pred_path,   show_col_types = FALSE, guess_max = 100000) %>% select(gee_id, tAGB_pred)
  df_ref  <- normalise_master_types(read_csv(master_path, show_col_types = FALSE)) %>% select(gee_id, AGB_pred, any_of(c("xcoord","ycoord")))

  df <- inner_join(df_ref, df_pred, by = "gee_id") %>%
    filter(!is.na(AGB_pred), !is.na(tAGB_pred), AGB_pred > 0)
  if (nrow(df) == 0) { message(sprintf("Year %s: no valid rows.", yr)); return(NULL) }
  cat(sprintf("Year %s: %d matched samples.\n", yr, nrow(df)))

  tibble(year = yr, n = nrow(df),
         RMSE = round(rmse_fn(df$AGB_pred, df$tAGB_pred), 3),
         MAE  = round(mae_fn(df$AGB_pred,  df$tAGB_pred), 3),
         R2   = round(r2_fn(df$AGB_pred,   df$tAGB_pred), 3))
})

write_csv(results_yearly, file.path(result_dir, "agb_performance_by_year.csv"))
cat("\n=== AGB Model Performance per Year ===\n"); print(results_yearly)

# =============================================================================
# BLOCK 2: Load all years combined
# =============================================================================
data_all <- purrr::map_dfr(data_years, function(yr) {
  pred_path   <- file.path(result_dir, paste0("predicted_tAGB_", yr, ".csv"))
  master_path <- file.path(result_dir, paste0("master_AGB_",     yr, ".csv"))
  if (!file.exists(pred_path) || !file.exists(master_path)) return(NULL)
  df_pred <- read_csv(pred_path,   show_col_types = FALSE, guess_max = 100000) %>% select(gee_id, tAGB_pred)
  df_ref  <- normalise_master_types(read_csv(master_path, show_col_types = FALSE)) %>% select(gee_id, AGB_pred, any_of(c("xcoord","ycoord")))
  inner_join(df_ref, df_pred, by = "gee_id") %>%
    filter(!is.na(AGB_pred), !is.na(tAGB_pred), AGB_pred > 0) %>% mutate(year = yr)
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
  labs(title = "AGB Model: RMSE and MAE per Year", x = "Year", y = "Error (gDW/m2)") +
  theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

# P2: R2 per year
print(ggplot(results_yearly, aes(x = factor(year, levels = data_years), y = R2)) +
  geom_col(fill = "grey50") + geom_text(aes(label = round(R2, 3)), vjust = -0.4, size = 3.5) +
  labs(title = "AGB Model R2 per Year", x = "Year", y = "R2") +
  theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

# P3: Scatter all years
print(ggplot(data_all, aes(x = AGB_pred, y = tAGB_pred, colour = factor(year, levels = data_years))) +
  geom_point(alpha = 0.3, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "AGB: Observed vs Predicted (All Years)",
       x = "Observed AGB (gDW/m2)", y = "Predicted AGB (gDW/m2)", colour = "Year") +
  theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

# P4: Scatter faceted
print(ggplot(data_all, aes(x = AGB_pred, y = tAGB_pred)) +
  geom_point(alpha = 0.3, size = 1) + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ year) +
  labs(title = "AGB: Observed vs Predicted per Year",
       x = "Observed AGB (gDW/m2)", y = "Predicted AGB (gDW/m2)") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold")))

# P5: Boxplot by predicted bin
print(data_all %>%
  mutate(Pred_Bin = cut(tAGB_pred, breaks = pretty(tAGB_pred, 10), include.lowest = TRUE)) %>%
  ggplot(aes(x = Pred_Bin, y = AGB_pred)) + geom_boxplot(fill = "grey85") +
  labs(title = "Actual AGB by Predicted Bin", x = "Predicted AGB Bin (gDW/m2)", y = "Actual AGB (gDW/m2)") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold")))

# P6: Absolute error by actual bin
print(data_all %>%
  mutate(Error = abs(AGB_pred - tAGB_pred),
         Actual_Bin = cut(AGB_pred, breaks = pretty(AGB_pred, 10), include.lowest = TRUE)) %>%
  ggplot(aes(x = Actual_Bin, y = Error)) + geom_boxplot(fill = "grey85") +
  labs(title = "Absolute Error by Actual AGB Bin", x = "Actual AGB Bin (gDW/m2)", y = "Absolute Error") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold")))

# P7: Residual distribution
print(data_all %>% mutate(resid = tAGB_pred - AGB_pred) %>%
  ggplot(aes(x = resid)) +
  geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7, colour = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  labs(title = "AGB Residual Distribution (Pred - Obs)", x = "Residual", y = "Count") +
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
    geom_point(data = df_map, aes(x = x, y = y, colour = tAGB_pred), size = 1.2, alpha = 0.7) +
    scale_colour_gradient(name = "Predicted AGB (gDW/m2)", low = "lightblue", high = "darkblue") +
    coord_quickmap(xlim = range(df_map$x) + c(-1,1), ylim = range(df_map$y) + c(-1,1)) +
    labs(title = "Spatial Distribution - Predicted AGB (All Years)", x = "Longitude", y = "Latitude") +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom"))
}

cat("\nE3_rf_evalModel_AGB.R COMPLETE\n")
