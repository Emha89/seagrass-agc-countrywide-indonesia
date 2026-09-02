# =============================================================================
# C3_rf_evalModel_cIndex.R
# Evaluate carbon index predictions per year (supplementary).
#
# Inputs:
#   predicted_CINDEX_YYYY.csv  (from C2_rf_applyModel_cIndex.R)
#   master_raw_YYYY.csv        (ground truth carbon_index)
#
# Outputs:
#   cindex_performance_by_year.csv
#   cindex_eval_points_all_years.csv
#   Plots displayed in panel (supplementary paper)
#
# Same "structural zeros excluded" filter as the training script
# (C1_rf_model_cIndex.R) is applied here too, for consistency between
# training and evaluation.
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

rmse_vec   <- function(p, o) sqrt(mean((p - o)^2, na.rm = TRUE))
mae_vec    <- function(p, o) mean(abs(p - o), na.rm = TRUE)
r2_vec     <- function(p, o) {
  ok <- is.finite(p) & is.finite(o); p <- p[ok]; o <- o[ok]
  if (length(o) < 2) return(NA_real_)
  ss_res <- sum((o - p)^2); ss_tot <- sum((o - mean(o))^2)
  if (ss_tot == 0) return(NA_real_); 1 - ss_res / ss_tot
}
bias_vec   <- function(p, o) mean(p - o, na.rm = TRUE)
safe_quant <- function(x, probs) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

# =============================================================================
# BLOCK 1: Per-year metrics
# =============================================================================
cat("Evaluating Carbon Index predictions per year...\n")

eval_points_all <- list()

perf_by_year <- purrr::map_dfr(data_years, function(yr) {
  pred_path   <- file.path(result_dir, paste0("predicted_CINDEX_", yr, ".csv"))
  master_path <- file.path(result_dir, paste0("master_raw_",       yr, ".csv"))
  if (!file.exists(pred_path) || !file.exists(master_path)) {
    message(sprintf("Year %s: missing file.", yr)); return(NULL)
  }
  df_pred  <- read_csv(pred_path,   show_col_types = FALSE, guess_max = 100000) %>%
    select(gee_id, carbon_index_pred, any_of(c("xcoord","ycoord")))
  df_truth <- normalise_master_types(read_csv(master_path, show_col_types = FALSE)) %>%
    select(gee_id, carbon_index, any_of(c("xcoord","ycoord")))
  
  if (!"carbon_index_pred" %in% names(df_pred)) {
    warning(sprintf("Year %s: carbon_index_pred missing.", yr)); return(NULL)
  }
  
  df_eval <- inner_join(df_truth, df_pred, by = "gee_id", suffix = c("","_p")) %>%
    mutate(carbon_index      = suppressWarnings(as.numeric(carbon_index)),
           carbon_index_pred = suppressWarnings(as.numeric(carbon_index_pred))) %>%
    filter(is.finite(carbon_index), is.finite(carbon_index_pred),
           carbon_index > 0, carbon_index <= 1) %>%
    mutate(year_eval = yr, resid = carbon_index_pred - carbon_index, abs_err = abs(resid))
  
  if (nrow(df_eval) == 0) { message(sprintf("Year %s: 0 valid rows.", yr)); return(NULL) }
  eval_points_all[[yr]] <<- df_eval
  
  q <- safe_quant(df_eval$abs_err, c(0.5, 0.75, 0.95))
  tibble::tibble(
    year = yr, n = nrow(df_eval),
    rmse        = round(rmse_vec(df_eval$carbon_index_pred, df_eval$carbon_index), 4),
    mae         = round(mae_vec(df_eval$carbon_index_pred,  df_eval$carbon_index), 4),
    r2          = round(r2_vec(df_eval$carbon_index_pred,   df_eval$carbon_index), 4),
    bias        = round(bias_vec(df_eval$carbon_index_pred, df_eval$carbon_index), 4),
    abs_err_p50 = round(q[1], 4), abs_err_p75 = round(q[2], 4), abs_err_p95 = round(q[3], 4)
  )
})

if (nrow(perf_by_year) == 0) stop("No evaluation results. Check input files.")
write_csv(perf_by_year, file.path(result_dir, "cindex_performance_by_year.csv"))
cat("\n=== Carbon Index Model Performance per Year ===\n"); print(perf_by_year)

points_df <- bind_rows(eval_points_all) %>% mutate(year_eval = factor(year_eval, levels = data_years))
write_csv(points_df, file.path(result_dir, "cindex_eval_points_all_years.csv"))
cat(sprintf("Saved: cindex_eval_points_all_years.csv (%d rows)\n\n", nrow(points_df)))

# =============================================================================
# BLOCK 2: Supplementary plots
# =============================================================================

# P1: Metrics by year
print(perf_by_year %>%
        select(year, rmse, mae, r2, bias) %>%
        pivot_longer(c(rmse, mae, r2, bias), names_to = "metric", values_to = "value") %>%
        mutate(year = factor(year, levels = data_years)) %>%
        ggplot(aes(x = year, y = value, fill = metric)) +
        geom_col(position = "dodge") +
        geom_text(aes(label = value), position = position_dodge(0.9), vjust = -0.35, size = 3) +
        labs(title = "Carbon Index Model Performance per Year", x = "Year", y = "Metric Value") +
        theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

# P2: Scatter all years
print(ggplot(points_df, aes(x = carbon_index, y = carbon_index_pred, colour = year_eval)) +
        geom_point(alpha = 0.4, size = 1.5) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        labs(title = "Carbon Index: Observed vs Predicted (All Years)",
             x = "Observed", y = "Predicted", colour = "Year") +
        theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

# P3: Scatter faceted by year
print(ggplot(points_df, aes(x = carbon_index, y = carbon_index_pred)) +
        geom_point(alpha = 0.35, size = 1.2) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        coord_equal(xlim = c(0,1), ylim = c(0,1)) + facet_wrap(~ year_eval) +
        labs(title = "Carbon Index: Observed vs Predicted per Year", x = "Observed", y = "Predicted") +
        theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold")))

# P4: Residual distribution
print(ggplot(points_df, aes(x = resid)) +
        geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7, colour = "white") +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
        labs(title = "Residual Distribution (Pred - Obs) - All Years", x = "Residual", y = "Count") +
        theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

# P5: Residuals vs observed (heteroscedasticity check)
print(ggplot(points_df, aes(x = carbon_index, y = resid)) +
        geom_point(alpha = 0.35, size = 1.4) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
        labs(title = "Residuals vs Observed Carbon Index", x = "Observed", y = "Residual (Pred - Obs)") +
        theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

# P6: Absolute error boxplot per year
print(ggplot(points_df, aes(x = year_eval, y = abs_err)) +
        geom_boxplot(outlier.size = 0.6, alpha = 0.7) +
        labs(title = "Absolute Error per Year - Carbon Index", x = "Year", y = "|Residual|") +
        theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

# P7: Density overlay obs vs pred
print(ggplot(points_df) +
        geom_density(aes(x = carbon_index,      colour = "Observed"),  linewidth = 1.1) +
        geom_density(aes(x = carbon_index_pred, colour = "Predicted"), linewidth = 1.1) +
        labs(title = "Distribution: Observed vs Predicted Carbon Index",
             x = "Carbon Index", y = "Density", colour = "") +
        theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

# P8: Calibration plot (binned means)
print(points_df %>%
        mutate(bin = cut(carbon_index, breaks = seq(0, 1, 0.1), include.lowest = TRUE)) %>%
        group_by(bin) %>%
        summarise(obs_mean = mean(carbon_index, na.rm = TRUE),
                  pred_mean = mean(carbon_index_pred, na.rm = TRUE), n = n(), .groups = "drop") %>%
        ggplot(aes(x = obs_mean, y = pred_mean)) +
        geom_point(size = 2.5) + geom_line() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
        coord_equal(xlim = c(0,1), ylim = c(0,1)) +
        labs(title = "Calibration Plot (Binned Means) - Carbon Index",
             x = "Mean Observed", y = "Mean Predicted") +
        theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold")))

# P9: Spatial map
x_col <- if ("xcoord" %in% names(points_df)) "xcoord" else "lon"
y_col <- if ("ycoord" %in% names(points_df)) "ycoord" else "lat"
if (all(c(x_col, y_col) %in% names(points_df))) {
  df_map <- points_df %>%
    mutate(x = suppressWarnings(as.numeric(.data[[x_col]])),
           y = suppressWarnings(as.numeric(.data[[y_col]]))) %>% filter(is.finite(x), is.finite(y))
  print(ggplot() +
          geom_polygon(data = ggplot2::map_data("world", region = "Indonesia"),
                       aes(x = long, y = lat, group = group), fill = "#d4e6c3", colour = "grey50", linewidth = 0.3) +
          geom_point(data = df_map, aes(x = x, y = y, colour = carbon_index_pred), size = 1.2, alpha = 0.8) +
          scale_colour_gradient(name = "Predicted Carbon Index", low = "#ffffb2", high = "#006837") +
          coord_quickmap(xlim = range(df_map$x) + c(-1,1), ylim = range(df_map$y) + c(-1,1)) +
          labs(title = "Spatial Distribution - Predicted Carbon Index (All Years)", x = "Longitude", y = "Latitude") +
          theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom"))
}

cat("\nC3_rf_evalModel_cIndex.R COMPLETE\n")
