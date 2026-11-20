# =============================================================================
# P4_05_tier_uncertainty.R
# Paper 4, design v3: uncertainty terms, Tier 1 national values and the Tier 2
# ecoregion strata with their leave-one-year-out check.
#
# Reads : P4_samples.csv, P4_oof_agc_cluster.csv, P4_oof_agc_pixel.csv
# Writes: P4_uncertainty_params.csv     ICC, DEFF, N_eff, SD_model, SD_data, bias SE
#         P4_sd_model_by_year.csv
#         P4_tier1_national.csv         per year and pooled, DEFF-corrected CI
#                                       (area CFG$extent_km2; final totals in P4_07)
#         P4_tier1_agc_bounds.csv       pooled Tier 1 mean from AGC_pred, AGC_low, AGC_up
#         P4_tier2_eco_strata.csv       ecoregion means, n, n_eff; stratum = n >= CFG$min_n_eco
#         P4_tier2_eco_loyo.csv         the same, leaving one survey year out
#         P4_year_composition.csv
# Tier 2 totals need the ecoregion areas of the GEE extent: P4_07 STEP 3.
# =============================================================================

source("R/P4_00_config_functions.R")

samples <- read_samples()
agc <- samples |> dplyr::filter(is_agc)
oof_cl  <- read_oof(designs = "cluster")
oof_pix <- read_oof(designs = "pixel")

# =============================================================================
# STEP 1: Design effect from out-of-fold residuals grouped by DBSCAN cluster
# (clustered points only, as in the thesis definition)
# =============================================================================
icc_from <- function(oof, label) {
  oof |>
    dplyr::filter(grepl("^c", g_cluster)) |>
    dplyr::mutate(resid = obs - pred_pixel) |>
    dplyr::group_by(rep) |>
    dplyr::group_modify(~ icc_oneway(.x$resid, .x$g_cluster)) |>
    dplyr::ungroup() |>
    dplyr::summarise(dplyr::across(c(N, k, n0, icc, deff), mean)) |>
    dplyr::mutate(source = label, .before = 1)
}
icc_tbl <- dplyr::bind_rows(icc_from(oof_cl, "cluster design (used)"),
                            icc_from(oof_pix, "pixel design (sensitivity)"))
print_tbl(icc_tbl)
DEFF <- icc_tbl$deff[1]
cat(sprintf("Design effect used: %.2f | N_eff of all AGC samples: %.0f\n", DEFF, nrow(agc) / DEFF))

# =============================================================================
# STEP 2: SD_model (cluster design RMSE), SD_data (RMS of reference SD)
# =============================================================================
sd_model <- oof_cl |>
  dplyr::group_by(rep) |>
  dplyr::summarise(RMSE = sqrt(mean((obs - pred_pixel)^2)), .groups = "drop") |>
  dplyr::summarise(SD_model = mean(RMSE)) |>
  dplyr::pull(SD_model)

sd_model_year <- oof_cl |>
  dplyr::group_by(rep, year) |>
  dplyr::summarise(n = dplyr::n(), RMSE = sqrt(mean((obs - pred_pixel)^2)), .groups = "drop") |>
  dplyr::group_by(year) |>
  dplyr::summarise(n = mean(n), SD_model = mean(RMSE), .groups = "drop") |>
  dplyr::mutate(use_year_value = n >= 30)
readr::write_csv(sd_model_year, out_file("P4_sd_model_by_year.csv"))
print_tbl(sd_model_year)

sd_data <- sqrt(mean(agc$sd_ref^2))

# Diagnostic: bias of the mean prediction with a cluster-robust SE.
# This is the error that does not shrink with the number of pixels.
bias_tbl <- oof_cl |>
  dplyr::mutate(resid = obs - pred_pixel) |>
  dplyr::group_by(rep) |>
  dplyr::group_modify(~ cluster_se_mean(.x$resid, .x$g_cluster)) |>
  dplyr::ungroup() |>
  dplyr::summarise(dplyr::across(c(mean_resid, se_cluster, n_groups), mean))

unc <- tibble::tribble(
  ~parameter, ~value, ~note,
  "ICC", icc_tbl$icc[1], "out-of-fold residuals, cluster design, clustered points only",
  "adjusted_mean_cluster_size", icc_tbl$n0[1], "n0 of the one-way ANOVA",
  "DEFF", DEFF, "1 + (n0 - 1) x ICC; applied to Tier 1, Tier 2 and pixel-based CIs",
  "N_agc", nrow(agc), "AGC samples",
  "N_eff", nrow(agc) / DEFF, "N_agc / DEFF",
  "SD_model_gC_m2", sd_model, "RMSE of the cluster design (between meadows)",
  "SD_data_gC_m2", sd_data, "root mean square of (AGC_up - AGC_low) / 3.92",
  "SD_total_gC_m2", sqrt(sd_model^2 + sd_data^2), "sqrt(SD_model^2 + SD_data^2)",
  "mean_resid_gC_m2", bias_tbl$mean_resid, "diagnostic: mean (obs - pred), cluster design",
  "se_mean_resid_gC_m2", bias_tbl$se_cluster, "diagnostic: cluster-robust SE of the mean residual",
  "bias_halfwidth_national_tC", 1.96 * bias_tbl$se_cluster * CFG$extent_km2,
  "diagnostic: 1.96 x SE x extent (km2) = +/- on the national pixel-based total"
)
print_tbl(unc)
readr::write_csv(unc, out_file("P4_uncertainty_params.csv"))

# =============================================================================
# STEP 3: Tier 1 (per survey year and pooled), DEFF-corrected CI
# =============================================================================
tier1 <- dplyr::bind_rows(
  agc |> dplyr::group_by(year) |>
    dplyr::group_modify(~ tier1_summary(.x$AGC_pred, DEFF, CFG$extent_km2)) |>
    dplyr::ungroup() |> dplyr::mutate(year = as.character(year)),
  tier1_summary(agc$AGC_pred, DEFF, CFG$extent_km2) |> dplyr::mutate(year = "pooled")
) |>
  dplyr::mutate(ci_note = ifelse(is.na(ci_half_gC_m2), "n_eff < 2, CI not defined", ""))
print_tbl(tier1, 2)
readr::write_csv(tier1, out_file("P4_tier1_national.csv"))

# Pooled Tier 1 from the lower and upper 95% bounds of the cover-based AGC of each
# sample: how much the error of the reference AGC moves the national mean
bounds <- purrr::map_dfr(c("AGC_pred", "AGC_low", "AGC_up"), function(v)
  tier1_summary(agc[[v]], DEFF, CFG$extent_km2) |> dplyr::mutate(response = v, .before = 1))
cat("\nPooled Tier 1 from AGC_pred and its 95% bounds (gC m-2):\n")
print_tbl(bounds |> dplyr::select(response, n, mean_gC_m2, sd_gC_m2, ci_half_gC_m2), 2)
readr::write_csv(bounds, out_file("P4_tier1_agc_bounds.csv"))

# =============================================================================
# STEP 4: Tier 2 ecoregion strata (IPCC-style geographic strata)
# stratum = ecoregion with >= CFG$min_n_eco AGC samples; the other ecoregions take
# the national mean. Survey year and region are confounded in this dataset, so
# STEP 5 checks how much one survey year moves each stratum mean.
# =============================================================================
strata <- function(d) {
  d |>
    dplyr::group_by(ecoregion) |>
    dplyr::summarise(n = dplyr::n(), n_years = dplyr::n_distinct(year),
                     mean_gC_m2 = mean(AGC_pred), sd_gC_m2 = sd(AGC_pred), .groups = "drop") |>
    dplyr::mutate(n_eff = n / DEFF, stratum = n >= CFG$min_n_eco,
                  ci_half_gC_m2 = ifelse(stratum & n_eff >= 2, 1.96 * sd_gC_m2 / sqrt(n_eff), NA_real_))
}
eco_strata <- strata(agc) |> dplyr::arrange(dplyr::desc(n))
cat(sprintf("\nTier 2 strata (n >= %d): %d of %d ecoregions with AGC samples; national mean %.2f gC m-2\n",
            CFG$min_n_eco, sum(eco_strata$stratum), nrow(eco_strata), mean(agc$AGC_pred)))
print_tbl(eco_strata, 2)
readr::write_csv(eco_strata, out_file("P4_tier2_eco_strata.csv"))

# =============================================================================
# STEP 5: Leave-one-survey-year-out
# =============================================================================
loyo <- purrr::map_dfr(sort(unique(agc$year)), function(y) {
  d <- agc |> dplyr::filter(year != y)
  dplyr::bind_rows(strata(d),
                   tibble::tibble(ecoregion = "National", n = nrow(d), n_years = dplyr::n_distinct(d$year),
                                  mean_gC_m2 = mean(d$AGC_pred), sd_gC_m2 = sd(d$AGC_pred),
                                  n_eff = nrow(d) / DEFF, stratum = TRUE, ci_half_gC_m2 = NA_real_)) |>
    dplyr::mutate(left_out_year = y, .before = 1)
})
readr::write_csv(loyo, out_file("P4_tier2_eco_loyo.csv"))
loyo_range <- loyo |>
  dplyr::filter(stratum) |>
  dplyr::group_by(ecoregion) |>
  dplyr::summarise(n_years_left_out = dplyr::n(), min_mean = min(mean_gC_m2), max_mean = max(mean_gC_m2),
                   .groups = "drop")
cat("\nStratum means when one survey year is left out (gC m-2):\n")
print_tbl(dplyr::left_join(eco_strata |> dplyr::filter(stratum) |> dplyr::select(ecoregion, mean_gC_m2) |>
                             dplyr::bind_rows(tibble::tibble(ecoregion = "National",
                                                             mean_gC_m2 = mean(agc$AGC_pred))),
                           loyo_range, by = "ecoregion"), 2)

# =============================================================================
# STEP 6: Survey composition by year and ecoregion (for the time-series figure)
# =============================================================================
comp <- agc |>
  dplyr::count(year, ecoregion, name = "n") |>
  dplyr::group_by(year) |>
  dplyr::mutate(share = n / sum(n)) |>
  dplyr::ungroup()
readr::write_csv(comp, out_file("P4_year_composition.csv"))

p_comp <- ggplot(comp, aes(x = factor(year), y = share, fill = ecoregion)) +
  geom_col(colour = "white", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Survey composition behind each Tier 1 annual mean",
       x = "Survey year", y = "Share of AGC samples", fill = "Ecoregion") +
  theme_minimal(base_size = 11)
print(p_comp)
ggsave(out_file("P4_fig_year_composition.png"), p_comp, width = 8, height = 5, dpi = 300)

cat("\nP4_05_tier_uncertainty.R COMPLETE\n")
