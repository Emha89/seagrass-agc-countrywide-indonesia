# =============================================================================
# P4_00_config_functions.R
# Paper 4, design v3: configuration and shared functions.
# Script set version 2026-11-14 (v3c) (printed when this file is sourced; all P4_
# scripts must come from the same set).
#
# Open the RStudio project (.Rproj) and run every script from the project root:
#   source("R/P4_01_build_training.R")  ...  source("R/P4_07_gee_results.R")
# Every path below is relative to the project root. Edit the CONFIG block only.
#
# Run order (README.md):
#   P4_01_build_training.R   samples, group IDs, ecoregion, DBSCAN, Sentinel-2 join
#   P4_02_tuning.R           RF hyperparameters (pixel-grouped CV): GSE, and GSE and
#                            Sentinel-2 on the samples with a clear Sentinel-2 pixel
#   P4_03_validation.R       leak-free CV for five designs, four methods (long run)
#   P4_04_metrics.R          sample-level and inventory-level metrics, shrinkage
#   P4_05_tier_uncertainty.R DEFF, SD terms, Tier 1, Tier 2 ecoregion strata
#   P4_06_gee_prepare.R      GEE upload table, MEOW shapefile, eco_id lookup
#   P4_07_gee_results.R      national and ecoregion totals from the GEE exports
#
# v3 (14 Nov 2026): Tier 2 = ecoregion strata (the cover model is removed);
# Sentinel-2 percentile predictors as a benchmark model (validation only);
# inventory-level validation for sites, meadows and ecoregions. The ROI step is
# removed. v3b/v3c: the Sentinel-2 annual mosaic is built like the GSE annual
# layer from L2A surface reflectance (P4_G11),
# a median-only Sentinel-2 check, and the local mean (meadow mean) benchmark. The sample build is unchanged, so P4_samples.csv keeps the row order of
# the P4v2 GEE training asset (P4_07 checks this) and the GSE results are the same.
# =============================================================================

# Packages (install once):
# install.packages(c("dplyr", "readr", "tidyr", "purrr", "ggplot2", "scales",
#                    "ranger", "sf", "dbscan", "FNN"))
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ranger)
  library(sf)
  library(dbscan)
  library(FNN)
})

# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------
CFG <- list(
  # ---- Paths (relative to the project root) ----------------------------------
  result_dir   = "data/master_raw",                # master_raw_YYYY.csv (field data, not shared)
  out_dir      = "output/P4v3",
  meow_path    = "data/MEOW/WCMC-036-MEOW-PPOW-2007-2012-NoCoast.shp",
  gee_dir      = "data/gee",                       # CSV exports of the GEE scripts (P4v2 maps)
  # GSE re-extracted from the current collection (P4_G00_extract_gse.js; columns gee_id,
  # year, GSE_A00..A63). These values replace the master_raw GSE values in P4_01, so the
  # training vectors and the GEE maps use the same dataset version (1.1 since 2025-11-17,
  # regenerated 2017 layer). NA = keep master_raw values.
  gse_path     = "data/gee/P4_gee_gse_extract.csv",

  # ---- Sentinel-2 benchmark (P4_G11_extract_s2.js) ----------------------------
  # One CSV per survey year: rid, gee_id, year, S2_<band>_p<pct>, S2_nobs (-9999 = no
  # clear observation). Product: L2A surface reflectance (the best-corrected level).
  # Annual mosaic built like the GSE annual layer: calendar year, every scene with a
  # Cloud Score+ pixel mask, one observation per pass, bands B2, B3, B4, B8, B11,
  # percentiles 10-90 (25 predictors). Samples without a clear L2A pixel in their
  # survey year (L2A gap in 2017-2018) have has_s2 = FALSE. Every method of the
  # Sentinel-2 comparison (Tier 1, Tier 2, local mean, GSE, Sentinel-2) is trained
  # and tested on the has_s2 samples only (P4_03).
  s2_dir       = "data/gee",
  s2_pattern   = "^P4_gee_s2_extract_\\d{4}\\.csv$",
  s2_bands     = c("B2", "B3", "B4", "B8", "B11"),
  s2_pcts      = c(10, 25, 50, 75, 90),
  s2_median_check = TRUE, # TRUE = also fit a median-only Sentinel-2 model (p50 of each
                          # band; the composite baseline of the AlphaEarth paper)
  s2_min_obs   = 1,       # clear observations needed in the survey year (S2_nobs)
  s2_years     = 2019:2023,  # survey years of the Sentinel-2 comparison (NULL = every year).
                          # P4_G11 console, 14 Nov 2026: L2A scenes over the samples were
                          # 9 of 143 L1C scenes in 2017 and 89 of 839 in 2018, but equal to
                          # L1C from 2019 on. The few 2017-2018 L2A samples would be
                          # scattered and thinly observed, so whole years are used.
  s2_equal_training = TRUE,   # TRUE = GSE refitted on the same Sentinel-2 subset in every
                              # fold, so GSE and Sentinel-2 see the same training samples
  run_s2       = TRUE,    # FALSE = skip the Sentinel-2 models in P4_02/P4_03

  # AGC measured in the field (input data of the cover-AGB model, "actual" rows), linked to gee_id by
  # P4_01a. measured_path = NA reads P4_measured_agc.csv from out_dir when it exists.
  p1_input_path    = "data/cover_agb/predict_tAGB2.csv",
  # Files for tracing the actual AGB (P4_01b): training table of the cover-AGB
  # model (01_model_dev.R) and the quadrat output of the apply script
  p1_train_path    = "data/cover_agb/data_long_format.rds",
  p1_summary_path  = "data/cover_agb/RF_total_AGB_AGC_summary_with_ratio_v4.csv",
  measured_methods = NULL,   # `method` values of that input that count as measured (NULL = all)
  link_max_dist_m  = 10,     # largest distance between that input and master_raw coordinates
  link_year_tol    = 1,      # composition_id links accepted when the years differ by up to this
                             # (one quadrat id = one quadrat; e.g. Rote 2021 vs 2022 in master_raw)
  measured_path    = NA_character_,
  use_measured     = FALSE,  # 15 Nov 2026: the "actual" values are not used in P4. They were
                             # compiled from available surveys, partly estimated (shoot density)
                             # and trained the cover-AGB model, so they are not an independent
                             # reference. P4_01a/P4_01b stay as optional diagnostics.

  # ---- Data ----------------------------------------------------------------
  years_train     = 2017:2023,
  gse_vars        = sprintf("GSE_A%02d", 0:63),
  # Environmental predictors: none (30 Oct 2026). distToLand was NA at Derawan (445 AGC
  # samples) and ACA depth was NA at 378 AGC samples (27% of 2018). Depth stays only in
  # the GEE extent mask (0-5 m), as before. Add names here to test them again.
  env_vars        = character(0),
  diag_cols       = c("depth", "distToLand"),   # kept for diagnostics only, never as predictors
  exclude_ids_path = NA_character_,             # optional CSV with column gee_id: samples to drop
                                                # (e.g. points on land from surveyor error)
  agc_response    = "AGC_pred",     # field point estimate (not a random draw)
  # Cover samples are kept in the sample build only, so P4_samples.csv and the GEE
  # training asset keep the same rows as in P4v2. No cover model is fitted in v3.
  cover_response  = "tSPC",
  cover_train_set = "present",      # "present" = tSPC > 0 only; "all" = include tSPC = 0
  cover_cap       = 100,            # tSPC above this (19 samples, up to 200 %) is set to 100
  gse_round       = 6,              # decimals used to find quadrats sharing one 10-m pixel

  # ---- Validation ----------------------------------------------------------
  designs       = c("pixel", "b100", "b1000", "cluster", "ecoregion"),
  block_sizes   = c(b100 = 100, b1000 = 1000),   # metres
  k_folds       = 10,
  n_repeats     = 10,
  n_seeds_eco   = 3,
  n_agc_strata  = 4,       # quantile strata of group mean AGC for fold assignment
  dbscan_eps    = 5000,    # metres
  dbscan_minPts = 5,
  min_n_eco     = 30,      # ecoregion: Tier 2 stratum, metrics and inventory unit above this n
  min_n_site    = 10,      # site (100 m block): AGC samples in its main survey year
  min_n_meadow  = 10,      # meadow (DBSCAN cluster): AGC samples
  min_n_local   = 10,      # local mean: training AGC samples in the test sample's meadow
                           # (DBSCAN cluster) needed; otherwise the Tier 2 value is used
  min_sites_eco = 5,       # site metrics per ecoregion only with at least this many sites
  overwrite     = FALSE,   # FALSE = skip designs whose OOF files already exist

  # ---- Random forest -------------------------------------------------------
  num_trees   = 500,
  num_threads = max(1, parallel::detectCores() - 1),
  tune_folds  = 5,
  grid_mtry   = c(5, 8, 12, 20),
  grid_node   = c(5, 10, 20),
  grid_frac   = c(0.5, 0.7, 0.9),

  # ---- National totals -----------------------------------------------------
  extent_km2 = 6895.80,    # thesis value (counted 100 m cells, about 2x too large; 2 Nov 2026).
                           # Used only by P4_05 (Tier 1 and bias half-width diagnostics);
                           # the final national tables come from P4_07 with the GEE pixel area.
  seed       = 42
)
# Optional override (testing or alternative paths): define CFG_USER before sourcing
if (exists("CFG_USER")) CFG <- modifyList(CFG, CFG_USER)

CFG$predictors <- c(CFG$gse_vars, CFG$env_vars)
CFG$s2_vars    <- sprintf("S2_%s_p%d", rep(CFG$s2_bands, each = length(CFG$s2_pcts)),
                          rep(CFG$s2_pcts, times = length(CFG$s2_bands)))
CFG$s2_med_vars <- sprintf("S2_%s_p50", CFG$s2_bands)   # median-only check
dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

EARTH_R <- 6371008.8   # mean Earth radius (m)

# Column that holds the group ID of each validation design
DESIGN_COL <- c(
  pixel = "g_pixel",
  setNames(paste0("g_", names(CFG$block_sizes)), names(CFG$block_sizes)),
  cluster   = "g_cluster",
  ecoregion = "ecoregion"
)

# Estimation methods: OOF column -> label (P4_04)
METHOD_LEVELS <- c("Tier 1", "Tier 2", "Local mean", "Pixel-based GSE", "Pixel-based S2",
                   "Pixel-based S2 (median)")

# -----------------------------------------------------------------------------
# File helpers
# -----------------------------------------------------------------------------
out_file <- function(name) file.path(CFG$out_dir, name)

# Read master_raw_YYYY.csv (from 01_build_master_raw.R) and keep only the
# columns this workflow needs, so years with differing column types still bind.
read_master <- function(years = CFG$years_train, dir = CFG$result_dir) {
  keep <- c("gee_id", "xcoord", "ycoord", "PA", "tSPC",
            "AGC_pred", "AGC_low", "AGC_up", CFG$predictors)
  purrr::map_dfr(years, function(yr) {
    p <- file.path(dir, paste0("master_raw_", yr, ".csv"))
    if (!file.exists(p)) { message("Missing: ", basename(p)); return(NULL) }
    df <- readr::read_csv(p, show_col_types = FALSE, guess_max = 100000)
    miss <- setdiff(keep, names(df))
    if (length(miss) > 0) stop(sprintf("%s lacks columns: %s", basename(p), paste(miss, collapse = ", ")))
    diag <- setdiff(CFG$diag_cols, keep)
    for (v in setdiff(diag, names(df))) df[[v]] <- NA_real_
    df <- df[, c(keep, diag)]
    for (v in setdiff(c(keep, diag), "gee_id")) df[[v]] <- suppressWarnings(as.numeric(df[[v]]))
    # -9999 = GEE fill code for "no pixel value" (unmask(-9999) in the extraction script):
    # treat it as missing so such samples are dropped like any other missing predictor
    for (v in intersect(c(CFG$predictors, diag), names(df))) df[[v]][!is.na(df[[v]]) & df[[v]] <= -9999] <- NA_real_
    df$gee_id <- as.character(df$gee_id)
    df$year   <- as.integer(yr)
    df
  })
}

read_samples <- function() {
  readr::read_csv(
    out_file("P4_samples.csv"), show_col_types = FALSE, guess_max = 100000,
    col_types = readr::cols(gee_id = "c", g_pixel = "c", g_b100 = "c", g_b1000 = "c",
                            g_cluster = "c", ecoregion = "c", has_s2 = "l",
                            .default = readr::col_guess())
  )
}

read_oof <- function(designs = CFG$designs) {
  purrr::map_dfr(designs, function(d) {
    p <- out_file(sprintf("P4_oof_agc_%s.csv", d))
    if (!file.exists(p)) { message("Missing: ", basename(p)); return(NULL) }
    readr::read_csv(p, show_col_types = FALSE, guess_max = 100000,
                    col_types = readr::cols(gee_id = "c", ecoregion = "c", g_cluster = "c",
                                            g_b100 = "c", design = "c", has_s2 = "l",
                                            t2_fallback = "l", local_fallback = "l",
                                            .default = readr::col_guess()))
  })
}

# model = "AGC" (GSE, all AGC samples), "AGC_GSE_SUB" (GSE, Sentinel-2 subset),
# "AGC_S2" (Sentinel-2 percentiles) or "AGC_S2MED" (Sentinel-2 median), both on the
# Sentinel-2 subset
read_best_params <- function(model = "AGC") {
  b <- readr::read_csv(out_file("P4_best_params.csv"), show_col_types = FALSE)
  b <- b[b$model == model, ]
  if (nrow(b) != 1) stop(sprintf("P4_best_params.csv has no row for model %s: run P4_02.", model))
  list(mtry = b$mtry, min_node_size = b$min_node_size, sample_fraction = b$sample_fraction)
}

# Print a table with numeric columns rounded (console only)
print_tbl <- function(x, digits = 3) {
  x <- as.data.frame(x)
  num <- vapply(x, is.numeric, logical(1))
  x[num] <- lapply(x[num], function(v) round(v, digits))
  print(x, row.names = FALSE)
  invisible(x)
}

# -----------------------------------------------------------------------------
# Geometry
# -----------------------------------------------------------------------------
lonlat_to_xyz <- function(lon, lat) {
  lo <- lon * pi / 180; la <- lat * pi / 180
  cbind(x = EARTH_R * cos(la) * cos(lo), y = EARTH_R * cos(la) * sin(lo), z = EARTH_R * sin(la))
}
chord_to_arc <- function(d) 2 * EARTH_R * asin(pmin(1, d / (2 * EARTH_R)))

# Great-circle distance (m) from each test point to its nearest training point
nn_distance_m <- function(test_lon, test_lat, train_lon, train_lat) {
  nn <- FNN::get.knnx(lonlat_to_xyz(train_lon, train_lat), lonlat_to_xyz(test_lon, test_lat), k = 1)
  chord_to_arc(nn$nn.dist[, 1])
}

# Square blocks on an equirectangular metric grid (x scale error < 2% within 11 deg of the equator)
block_id <- function(lon, lat, size_m) {
  x <- EARTH_R * lon * pi / 180
  y <- EARTH_R * lat * pi / 180
  paste(floor(x / size_m), floor(y / size_m), sep = "_")
}

# Quadrats of the same survey year with identical GSE values share one 10-m pixel
pixel_key <- function(df) {
  m <- round(as.matrix(df[, CFG$gse_vars]), CFG$gse_round)
  paste(df$year, apply(m, 1, paste, collapse = "_"), sep = "|")
}

# Assign each point to a MEOW ecoregion (NoCoast version); nearest polygon if none intersects
join_ecoregion <- function(df, meow_path = CFG$meow_path) {
  sf::sf_use_s2(FALSE)
  meow <- sf::st_read(meow_path, quiet = TRUE)
  meow <- meow[meow$TYPE == "MEOW", "ECOREGION"]
  meow <- sf::st_transform(meow, 4326)
  meow <- suppressWarnings(suppressMessages(
    sf::st_crop(meow, xmin = 94, ymin = -12, xmax = 142, ymax = 7)))
  pts <- sf::st_as_sf(as.data.frame(df[, c("gee_id", "xcoord", "ycoord")]),
                      coords = c("xcoord", "ycoord"), crs = 4326, remove = FALSE)
  j <- suppressMessages(sf::st_join(pts, meow, join = sf::st_intersects))
  j <- j[!duplicated(j$gee_id), ]
  miss <- is.na(j$ECOREGION)
  if (any(miss)) {
    nf <- suppressMessages(sf::st_nearest_feature(j[miss, ], meow))
    j$ECOREGION[miss] <- meow$ECOREGION[nf]
  }
  message(sprintf("Ecoregion join: %d points, %d assigned to the nearest polygon", nrow(j), sum(miss)))
  setNames(j$ECOREGION, j$gee_id)
}

# -----------------------------------------------------------------------------
# Folds: size-balanced group k-fold, stratified by group mean AGC
# groups = data.frame(group, n_agc, mean_agc, n_cov)
# Groups with AGC samples are balanced on AGC counts; cover-only groups on cover counts.
# (Cover-only groups are kept so the folds of the AGC groups stay those of P4v2.)
# -----------------------------------------------------------------------------
assign_folds <- function(groups, k, seed, n_strata = CFG$n_agc_strata) {
  set.seed(seed)
  fold <- setNames(integer(nrow(groups)), groups$group)
  load_agc <- numeric(k); load_cov <- numeric(k)
  pick <- function(load) { m <- which(load == min(load)); if (length(m) > 1) sample(m, 1) else m }

  agc <- groups[groups$n_agc > 0, , drop = FALSE]
  if (nrow(agc) > 0) {
    ns <- max(1, min(n_strata, floor(nrow(agc) / k)))
    agc$stratum <- dplyr::ntile(agc$mean_agc, ns)
    for (s in sort(unique(agc$stratum))) {
      idx <- which(agc$stratum == s)
      idx <- idx[sample.int(length(idx))]
      for (i in idx) {
        f <- pick(load_agc)
        fold[agc$group[i]] <- f
        load_agc[f] <- load_agc[f] + agc$n_agc[i]
      }
    }
  }
  cov <- groups[groups$n_agc == 0 & groups$n_cov > 0, , drop = FALSE]
  if (nrow(cov) > 0) {
    for (i in sample.int(nrow(cov))) {
      f <- pick(load_cov)
      fold[cov$group[i]] <- f
      load_cov[f] <- load_cov[f] + cov$n_cov[i]
    }
  }
  fold
}

# -----------------------------------------------------------------------------
# Models
# -----------------------------------------------------------------------------
fit_rf <- function(df, response, p, seed, vars = CFG$predictors) {
  ranger::ranger(
    x = as.data.frame(df[, vars]), y = df[[response]],
    num.trees = CFG$num_trees, mtry = p$mtry, min.node.size = p$min_node_size,
    sample.fraction = p$sample_fraction, replace = FALSE,
    num.threads = CFG$num_threads, seed = seed
  )
}
# The predictor set is read from the fitted model
predict_rf <- function(model, df) {
  vars <- model$forest$independent.variable.names
  predict(model, data = as.data.frame(df[, vars]), num.threads = CFG$num_threads)$predictions
}

# Tier 2: mean AGC per ecoregion where the ecoregion has at least min_n samples,
# otherwise the national mean (IPCC-style geographic strata)
tier2_fit <- function(y, eco, min_n = CFG$min_n_eco) {
  n <- tapply(y, eco, length); m <- tapply(y, eco, mean)
  m <- setNames(as.numeric(m), names(m))   # plain named vector, not a 1-d array
  list(national = mean(y), strata = m[as.numeric(n) >= min_n])
}
tier2_predict <- function(fit, eco) {
  v <- as.numeric(unname(fit$strata[as.character(eco)]))
  fb <- is.na(v)
  v[fb] <- fit$national
  list(pred = v, fallback = fb)
}

# Local mean: mean AGC of the training samples of the same meadow (DBSCAN cluster,
# ids starting with "c") when the meadow has >= min_n training samples; otherwise
# the Tier 2 value (ecoregion stratum, else national mean)
local_predict <- function(y_tr, g_tr, g_te, t2_te, min_n = CFG$min_n_local) {
  ok <- grepl("^c", g_tr)
  n <- tapply(y_tr[ok], g_tr[ok], length); m <- tapply(y_tr[ok], g_tr[ok], mean)
  m <- setNames(as.numeric(m), names(m))[as.numeric(n) >= min_n]
  v <- as.numeric(unname(m[as.character(g_te)]))
  fb <- is.na(v)
  v[fb] <- t2_te[fb]
  list(pred = v, fallback = fb)
}

# -----------------------------------------------------------------------------
# Metrics
# -----------------------------------------------------------------------------
reg_metrics <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]
  sse <- sum((obs - pred)^2); sst <- sum((obs - mean(obs))^2)
  tibble::tibble(
    n        = length(obs),
    RMSE     = sqrt(mean((obs - pred)^2)),
    MAE      = mean(abs(obs - pred)),
    R2       = if (length(obs) > 2 && sd(pred) > 0) cor(obs, pred)^2 else NA_real_,
    MEC      = if (sst > 0) 1 - sse / sst else NA_real_,
    bias_pct = (mean(pred) - mean(obs)) / mean(obs) * 100,
    slope    = if (length(obs) > 2 && sst > 0) stats::cov(obs, pred) / stats::var(obs) else NA_real_
  )
}

# One-way ANOVA intraclass correlation and design effect (cluster sampling)
icc_oneway <- function(y, g) {
  ok <- is.finite(y) & !is.na(g)
  y <- y[ok]; g <- as.character(g[ok])
  N <- length(y); ni <- table(g); k <- length(ni)
  mu <- tapply(y, g, mean); grand <- mean(y)
  msb <- sum(ni * (mu[names(ni)] - grand)^2) / (k - 1)
  msw <- sum((y - mu[g])^2) / (N - k)
  n0  <- (N - sum(ni^2) / N) / (k - 1)
  icc <- (msb - msw) / (msb + (n0 - 1) * msw)
  tibble::tibble(N = N, k = k, n0 = n0, icc = icc, deff = 1 + (n0 - 1) * max(icc, 0))
}

# Mean of residuals with a cluster-robust standard error
cluster_se_mean <- function(r, g) {
  ok <- is.finite(r) & !is.na(g)
  r <- r[ok]; g <- as.character(g[ok])
  N <- length(r); rbar <- mean(r)
  s <- tapply(r - rbar, g, sum); G <- length(s)
  tibble::tibble(mean_resid = rbar, se_cluster = sqrt(G / (G - 1) * sum(s^2)) / N, n_groups = G)
}

# Tier 1: mean x area with a clustering-corrected CI (tC = gC m-2 x km2)
tier1_summary <- function(y, deff, area_km2) {
  n <- length(y); m <- mean(y); s <- sd(y); neff <- n / deff
  half <- if (neff >= 2) 1.96 * s / sqrt(neff) else NA_real_
  tibble::tibble(n = n, n_eff = neff, mean_gC_m2 = m, sd_gC_m2 = s, ci_half_gC_m2 = half,
                 total_tC = m * area_km2,
                 ci_low_tC = (m - half) * area_km2, ci_up_tC = (m + half) * area_km2)
}

# Tier 2 total: ecoregion strata x ecoregion area with a stratified CI.
# areas: data.frame(ecoregion, area_km2). Strata = ecoregions with >= min_n samples;
# the other ecoregions use the national mean, and their areas form one pooled term
# of the variance. n_eff = n / DEFF in every stratum (national DEFF).
tier2_total <- function(areas, y, eco, deff, min_n = CFG$min_n_eco) {
  st <- tibble::tibble(ecoregion = eco, y = y) |>
    dplyr::group_by(ecoregion) |>
    dplyr::summarise(n = dplyr::n(), mean_gC_m2 = mean(y), sd_gC_m2 = sd(y), .groups = "drop")
  nat <- tibble::tibble(n = length(y), mean_gC_m2 = mean(y), sd_gC_m2 = sd(y))
  d <- areas |>
    dplyr::left_join(st, by = "ecoregion") |>
    dplyr::mutate(n = dplyr::coalesce(n, 0L),
                  stratum = n >= min_n,
                  dens_gC_m2 = ifelse(stratum, mean_gC_m2, nat$mean_gC_m2),
                  total_tC = dens_gC_m2 * area_km2,
                  n_eff = n / deff,
                  ci_half_tC = ifelse(stratum, 1.96 * sd_gC_m2 / sqrt(n_eff) * area_km2,
                                      1.96 * nat$sd_gC_m2 / sqrt(nat$n / deff) * area_km2))
  var_strata <- sum((d$area_km2^2 * d$sd_gC_m2^2 / d$n_eff)[d$stratum])
  a_fb <- sum(d$area_km2[!d$stratum])
  var_fb <- a_fb^2 * nat$sd_gC_m2^2 / (nat$n / deff)
  list(total = tibble::tibble(total_tC = sum(d$total_tC),
                              ci_half_tC = 1.96 * sqrt(var_strata + var_fb),
                              n_strata = sum(d$stratum),
                              area_strata_pct = 100 * sum(d$area_km2[d$stratum]) / sum(d$area_km2),
                              area_fallback_km2 = a_fb),
       by_ecoregion = d)
}

P4_SCRIPT_VERSION <- "2026-11-14 (v3c)"
cat(sprintf("P4 config loaded | script set %s | out_dir: %s | threads: %d | gse_path: %s | S2 files: %d\n",
            P4_SCRIPT_VERSION, CFG$out_dir, CFG$num_threads,
            if (is.na(CFG$gse_path)) "none (master_raw GSE)" else
              paste0(basename(CFG$gse_path), if (!file.exists(CFG$gse_path)) " (NOT FOUND)" else ""),
            length(list.files(CFG$s2_dir, pattern = CFG$s2_pattern))))
