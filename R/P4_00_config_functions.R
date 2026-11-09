# =============================================================================
# P4_00_config_functions.R
# Paper 4, revised design: configuration and shared functions.
# Script set version 2026-10-31 (printed when this file is sourced; all P4_
# scripts must come from the same set).
#
# Sourced at the top of P4_01 ... P4_07. Edit the CONFIG block only.
# Run order: see README.md of the repository
#   P4_01_build_training.R   training tables, group IDs, ecoregion, DBSCAN
#   P4_02_tuning.R           RF hyperparameters (pixel-grouped CV)
#   P4_03_validation.R       leak-free CV for five designs (long run)
#   P4_04_metrics.R          metrics, figures, cover accuracy, shrinkage
#   P4_05_tier_uncertainty.R Tier 1 / Tier 2 tables, N_eff, SD terms
#   P4_06_gee_prepare.R      upload table for GEE, MEOW shapefile
#   P4_07_gee_results.R      tables from the GEE exports
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
  # ---- Paths ---------------------------------------------------------------
  # Relative to the working directory (the R/ folder of the repository).
  # Edit these, or define CFG_USER with your own paths before sourcing this file.
  result_dir   = "../result",                # inputs: master_raw_<year>.csv (field data, not shared)
  out_dir      = "../result/P4v2",           # all P4_ outputs
  meow_path    = "../data/MEOW/WCMC-036-MEOW-PPOW-2007-2012-NoCoast.shp",   # Marine Ecoregions of the World (UNEP-WCMC)
  roi_path     = NA_character_,   # optional ROI polygons (shp/gpkg/geojson); NA = skip ROI check
  roi_id_field = "loc_id",

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
  # GSE re-extracted from the current collection (P4_G00_extract_gse.js; columns gee_id, year,
  # GSE_A00..A63). When set, these values replace the master_raw GSE values in P4_01, so the
  # training vectors and the GEE maps use the same dataset version (1.1 since 2025-11-17,
  # regenerated 2017 layer). NA = keep master_raw values.
  gse_path        = NA_character_,
  agc_response    = "AGC_pred",     # field point estimate (not a random draw)
  cover_response  = "tSPC",
  cover_train_set = "present",      # "present" = tSPC > 0 only; "all" = include tSPC = 0
  cover_cap       = 100,            # tSPC above this (19 samples, up to 200 %) is set to 100
                                    # in P4_01 (31 Oct 2026); the GEE cover model trains on
                                    # the same capped values via P4_06
  gse_round       = 6,              # decimals used to find quadrats sharing one 10-m pixel

  # ---- Cover classes (KepMenLH 200/2004) -----------------------------------
  cover_breaks = c(-Inf, 30, 60, Inf),
  cover_labels = c("<30", "30-59.9", ">=60"),

  # ---- Validation ----------------------------------------------------------
  designs       = c("pixel", "b100", "b1000", "cluster", "ecoregion"),
  block_sizes   = c(b100 = 100, b1000 = 1000),   # metres
  k_folds       = 10,
  n_repeats     = 10,
  n_seeds_eco   = 3,
  n_agc_strata  = 4,       # quantile strata of group mean AGC for fold assignment
  dbscan_eps    = 5000,    # metres
  dbscan_minPts = 5,
  min_n_eco     = 30,      # per-ecoregion metrics reported only above this n
  min_n_cell    = 20,      # Tier 2 class x year table: cells flagged below this n
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
dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

EARTH_R <- 6371008.8   # mean Earth radius (m)

# Column that holds the group ID of each validation design
DESIGN_COL <- c(
  pixel = "g_pixel",
  setNames(paste0("g_", names(CFG$block_sizes)), names(CFG$block_sizes)),
  cluster   = "g_cluster",
  ecoregion = "ecoregion"
)

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
                            g_cluster = "c", ecoregion = "c", .default = readr::col_guess())
  )
}

read_oof <- function(kind = c("agc", "cov"), designs = CFG$designs) {
  kind <- match.arg(kind)
  purrr::map_dfr(designs, function(d) {
    p <- out_file(sprintf("P4_oof_%s_%s.csv", kind, d))
    if (!file.exists(p)) { message("Missing: ", basename(p)); return(NULL) }
    ct <- if (kind == "agc") {
      readr::cols(gee_id = "c", ecoregion = "c", g_cluster = "c", cov_class = "c",
                  design = "c", .default = readr::col_guess())
    } else {
      readr::cols(gee_id = "c", ecoregion = "c", design = "c", .default = readr::col_guess())
    }
    readr::read_csv(p, show_col_types = FALSE, guess_max = 100000, col_types = ct)
  })
}

read_best_params <- function(model = c("AGC", "COVER")) {
  model <- match.arg(model)
  b <- readr::read_csv(out_file("P4_best_params.csv"), show_col_types = FALSE)
  b <- b[b$model == model, ]
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
fit_rf <- function(df, response, p, seed) {
  ranger::ranger(
    x = as.data.frame(df[, CFG$predictors]), y = df[[response]],
    num.trees = CFG$num_trees, mtry = p$mtry, min.node.size = p$min_node_size,
    sample.fraction = p$sample_fraction, replace = FALSE,
    num.threads = CFG$num_threads, seed = seed
  )
}
predict_rf <- function(model, df) {
  predict(model, data = as.data.frame(df[, CFG$predictors]), num.threads = CFG$num_threads)$predictions
}
cover_class <- function(x) {
  cut(x, breaks = CFG$cover_breaks, labels = CFG$cover_labels, right = FALSE)
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
    bias_pct = (mean(pred) - mean(obs)) / mean(obs) * 100
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

P4_SCRIPT_VERSION <- "2026-10-31"
cat(sprintf("P4 config loaded | script set %s | out_dir: %s | threads: %d | gse_path: %s\n",
            P4_SCRIPT_VERSION, CFG$out_dir, CFG$num_threads,
            if (is.na(CFG$gse_path)) "none (master_raw GSE)" else basename(CFG$gse_path)))
