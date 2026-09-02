# =============================================================================
# F3_rf_evalModel_AGC.R
# Evaluate AGC predictions per year (supplementary) -- final script in the
# R-side proxy chain (A through F).
#
# Inputs:
#   predicted_tAGC_YYYY.csv  (from F2_rf_applyModel_AGC.R)
#   master_AGC_YYYY.csv      (ground truth AGC_pred from 04_build_master_AGB_AGC.R)
#
# Outputs:
#   agc_performance_by_year.csv                          (BLOCK 1, primary)
#   agc_performance_by_year_CI75.csv                      (BLOCK 4, supplementary)
#   agc_performance_comparison_ALL_vs_CI75.csv            (BLOCK 4)
#   fraction_pixels_above_60_by_year.csv                  (BLOCK 5)
#   fraction_pixels_above_60_overall.csv                  (BLOCK 5)
#   Plots displayed in panel (supplementary paper)
#
# BLOCK 5 is the source of the paper's stated finding that prediction
# accuracy is stable below 60 gC/m2 but declines above it due to sparse
# high-density training samples, with pixels above that threshold
# excluded from the total AGC computation -- this block computes the
# fraction of predicted pixels (across the full national grid, not just
# field-validation points) that fall above that 60 gC/m2 threshold, per
# year and overall.
#
# BLOCK 4 is explicitly supplementary (a CI-width<=75 filtered subset,
# for transparency) -- the primary reported evaluation is BLOCK 1 (all
# valid samples, no CI filter).
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
# BLOCK 1: Per-year metrics (primary reported evaluation)
# =============================================================================
cat("Evaluating AGC predictions per year...\n")

results_yearly <- purrr::map_dfr(data_years, function(yr) {
  pred_path   <- file.path(result_dir, paste0("predicted_tAGC_", yr, ".csv"))
  master_path <- file.path(result_dir, paste0("master_AGC_",     yr, ".csv"))
  if (!file.exists(pred_path) || !file.exists(master_path)) {
    message(sprintf("Year %s: missing file.", yr)); return(NULL)
  }
  df_pred <- read_csv(pred_path,   show_col_types = FALSE, guess_max = 100000) %>% select(gee_id, tAGC_pred)
  df_ref  <- normalise_master_types(read_csv(master_path, show_col_types = FALSE)) %>% select(gee_id, AGC_pred, any_of(c("xcoord","ycoord")))

  df <- inner_join(df_ref, df_pred, by = "gee_id") %>%
    filter(!is.na(AGC_pred), !is.na(tAGC_pred), AGC_pred > 0)
  if (nrow(df) == 0) { message(sprintf("Year %s: no valid rows.", yr)); return(NULL) }
  cat(sprintf("Year %s: %d matched samples.\n", yr, nrow(df)))

  tibble(year = yr, n = nrow(df),
         RMSE = round(rmse_fn(df$AGC_pred, df$tAGC_pred), 3),
         MAE  = round(mae_fn(df$AGC_pred,  df$tAGC_pred), 3),
         R2   = round(r2_fn(df$AGC_pred,   df$tAGC_pred), 3))
})

write_csv(results_yearly, file.path(result_dir, "agc_performance_by_year.csv"))
cat("\n=== AGC Model Performance per Year ===\n"); print(results_yearly)

# =============================================================================
# BLOCK 2: Load all years combined
# =============================================================================
data_all <- purrr::map_dfr(data_years, function(yr) {
  pred_path   <- file.path(result_dir, paste0("predicted_tAGC_", yr, ".csv"))
  master_path <- file.path(result_dir, paste0("master_AGC_",     yr, ".csv"))
  if (!file.exists(pred_path) || !file.exists(master_path)) return(NULL)
  df_pred <- read_csv(pred_path,   show_col_types = FALSE, guess_max = 100000) %>% select(gee_id, tAGC_pred)
  df_ref  <- normalise_master_types(read_csv(master_path, show_col_types = FALSE)) %>% select(gee_id, AGC_pred, any_of(c("xcoord","ycoord")))
  inner_join(df_ref, df_pred, by = "gee_id") %>%
    filter(!is.na(AGC_pred), !is.na(tAGC_pred), AGC_pred > 0) %>% mutate(year = yr)
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
  labs(title = "AGC Model: RMSE and MAE per Year", x = "Year", y = "Error (gC/m2)") +
  theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

# P2: R2 per year
print(ggplot(results_yearly, aes(x = factor(year, levels = data_years), y = R2)) +
  geom_col(fill = "grey50") + geom_text(aes(label = round(R2, 3)), vjust = -0.4, size = 3.5) +
  labs(title = "AGC Model R2 per Year", x = "Year", y = "R2") +
  theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

# P3: Scatter all years
print(ggplot(data_all, aes(x = AGC_pred, y = tAGC_pred, colour = factor(year, levels = data_years))) +
  geom_point(alpha = 0.3, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "AGC: Observed vs Predicted (All Years)",
       x = "Observed AGC (gC/m2)", y = "Predicted AGC (gC/m2)", colour = "Year") +
  theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

# P4: Scatter faceted
print(ggplot(data_all, aes(x = AGC_pred, y = tAGC_pred)) +
  geom_point(alpha = 0.3, size = 1) + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ year) +
  labs(title = "AGC: Observed vs Predicted per Year",
       x = "Observed AGC (gC/m2)", y = "Predicted AGC (gC/m2)") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold")))

# P5: Boxplot by predicted bin
print(data_all %>%
  mutate(Pred_Bin = cut(tAGC_pred, breaks = pretty(tAGC_pred, 10), include.lowest = TRUE)) %>%
  ggplot(aes(x = Pred_Bin, y = AGC_pred)) + geom_boxplot(fill = "grey85") +
  labs(title = "Actual AGC by Predicted Bin", x = "Predicted AGC Bin (gC/m2)", y = "Actual AGC (gC/m2)") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold")))

# P6: Absolute error by actual bin
print(data_all %>%
  mutate(Error = abs(AGC_pred - tAGC_pred),
         Actual_Bin = cut(AGC_pred, breaks = pretty(AGC_pred, 10), include.lowest = TRUE)) %>%
  ggplot(aes(x = Actual_Bin, y = Error)) + geom_boxplot(fill = "grey85") +
  labs(title = "Absolute Error by Actual AGC Bin", x = "Actual AGC Bin (gC/m2)", y = "Absolute Error") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold")))

# P7: Residual distribution
print(data_all %>% mutate(resid = tAGC_pred - AGC_pred) %>%
  ggplot(aes(x = resid)) +
  geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7, colour = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  labs(title = "AGC Residual Distribution (Pred - Obs)", x = "Residual", y = "Count") +
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
    geom_point(data = df_map, aes(x = x, y = y, colour = tAGC_pred), size = 1.2, alpha = 0.7) +
    scale_colour_gradient(name = "Predicted AGC (gC/m2)", low = "#ffffb2", high = "#006837") +
    coord_quickmap(xlim = range(df_map$x) + c(-1,1), ylim = range(df_map$y) + c(-1,1)) +
    labs(title = "Spatial Distribution - Predicted AGC (All Years)", x = "Longitude", y = "Latitude") +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom"))
}

cat("\nBLOCKS 1-3 COMPLETE\n")

# =============================================================================
# BLOCK 4: Supplementary comparison -- CI width <= 75 subset
# Evaluates model performance on samples matching the training CI filter
# (AGC_CIwidt <= 75) for transparency. Not the primary reported metric.
# Primary evaluation remains BLOCK 1 (all samples).
# =============================================================================
cat("\n=== BLOCK 4: Supplementary -- CI width <= 75 subset ===\n")

results_ci75 <- purrr::map_dfr(data_years, function(yr) {
  pred_path   <- file.path(result_dir, paste0("predicted_tAGC_", yr, ".csv"))
  master_path <- file.path(result_dir, paste0("master_AGC_",     yr, ".csv"))
  if (!file.exists(pred_path) || !file.exists(master_path)) return(NULL)
  
  df_pred <- read_csv(pred_path,   show_col_types = FALSE, guess_max = 100000) %>%
    select(gee_id, tAGC_pred)
  df_ref  <- normalise_master_types(read_csv(master_path, show_col_types = FALSE)) %>%
    select(gee_id, AGC_pred, AGC_CIwidt, any_of(c("xcoord","ycoord")))
  
  df <- inner_join(df_ref, df_pred, by = "gee_id") %>%
    filter(!is.na(AGC_pred), !is.na(tAGC_pred), AGC_pred > 0) %>%
    filter(!is.na(AGC_CIwidt), AGC_CIwidt <= 75)   # CI filter matching training
  
  if (nrow(df) < 2) {
    message(sprintf("Year %s: insufficient rows after CI filter (n=%d).", yr, nrow(df)))
    return(NULL)
  }
  cat(sprintf("Year %s: n=%d (CI<=75)\n", yr, nrow(df)))
  
  tibble(year = yr, n = nrow(df),
         RMSE = round(rmse_fn(df$AGC_pred, df$tAGC_pred), 3),
         MAE  = round(mae_fn(df$AGC_pred,  df$tAGC_pred), 3),
         R2   = round(r2_fn(df$AGC_pred,   df$tAGC_pred), 3))
})

write_csv(results_ci75, file.path(result_dir, "agc_performance_by_year_CI75.csv"))
cat("\n=== AGC Performance (CI width <= 75 subset) ===\n"); print(results_ci75)

# Comparison table: ALL vs CI75
cat("\n=== Comparison: ALL samples vs CI<=75 subset ===\n")
comparison <- left_join(
  results_yearly %>% select(year, n_all = n, RMSE_all = RMSE, MAE_all = MAE, R2_all = R2),
  results_ci75   %>% select(year, n_ci75 = n, RMSE_ci75 = RMSE, MAE_ci75 = MAE, R2_ci75 = R2),
  by = "year"
)
print(comparison)
write_csv(comparison, file.path(result_dir, "agc_performance_comparison_ALL_vs_CI75.csv"))

cat("\nBLOCK 4 COMPLETE\n")

# =============================================================================
# BLOCK 5: Fraction of national predictions above the 60 gC m-2 threshold
# Uses the FULL predicted_tAGC_YYYY.csv files (all predicted pixels), not the
# field-matched subset used in Blocks 1-4, since the >60 exclusion is applied
# to the entire national prediction surface. This is the source of the
# paper's stated finding about prediction reliability declining above
# 60 gC/m2 and those pixels being excluded from the total AGC computation.
# =============================================================================
cat("\n=== BLOCK 5: Fraction of predictions above 60 gC m-2 threshold ===\n")

threshold <- 60

results_above60 <- purrr::map_dfr(data_years, function(yr) {
  pred_path <- file.path(result_dir, paste0("predicted_tAGC_", yr, ".csv"))
  if (!file.exists(pred_path)) {
    message(sprintf("Year %s: missing prediction file.", yr))
    return(NULL)
  }
  
  df_pred <- read_csv(pred_path, show_col_types = FALSE, guess_max = 100000) %>%
    select(gee_id, tAGC_pred) %>%
    filter(!is.na(tAGC_pred))
  
  n_total   <- nrow(df_pred)
  n_above   <- sum(df_pred$tAGC_pred > threshold, na.rm = TRUE)
  pct_above <- round(100 * n_above / n_total, 3)
  
  cat(sprintf("Year %s: %d / %d pixels above %g gC/m2 (%.3f%%)\n",
              yr, n_above, n_total, threshold, pct_above))
  
  tibble(year = yr, n_total = n_total, n_above_60 = n_above, pct_above_60 = pct_above)
})

write_csv(results_above60, file.path(result_dir, "fraction_pixels_above_60_by_year.csv"))

# National aggregate across all years
overall_total <- sum(results_above60$n_total, na.rm = TRUE)
overall_above <- sum(results_above60$n_above_60, na.rm = TRUE)
overall_pct   <- round(100 * overall_above / overall_total, 3)

cat(sprintf("\n=== OVERALL (all years combined): %d / %d pixels above %g gC/m2 (%.3f%%) ===\n",
            overall_above, overall_total, threshold, overall_pct))

write_csv(
  tibble(n_total = overall_total, n_above_60 = overall_above, pct_above_60 = overall_pct),
  file.path(result_dir, "fraction_pixels_above_60_overall.csv")
)

cat("\nBLOCK 5 COMPLETE\n")
cat("\nF3_rf_evalModel_AGC.R COMPLETE\n")
