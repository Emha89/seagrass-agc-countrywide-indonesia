# =============================================================================
# P4_01_build_training.R
# Paper 4, design v3: GSE training tables, validation groups and the
# Sentinel-2 predictors of the benchmark model.
#
# Reads : data/master_raw/master_raw_YYYY.csv (2017-2023) from 01_build_master_raw.R
#         MEOW shapefile (NoCoast version)
#         CFG$gse_path  GSE of the current dataset version (P4_G00)
#         CFG$s2_dir    P4_gee_s2_extract_YYYY.csv (P4_G11), optional
#         P4_measured_agc.csv  measured AGC linked to gee_id (P4_01a; optional)
# Writes (CFG$out_dir):
#   P4_predictor_na.csv         samples that missing predictors would remove
#   P4_agc_without_env.csv      AGC samples without ACA depth or distToLand (QGIS check)
#   P4_samples.csv              all samples of the sample build, with group IDs for
#                               every validation design, Sentinel-2 predictors and has_s2
#   (GEE upload table: P4_06_gee_prepare.R)
#   P4_counts_ecoregion_year.csv
#   P4_pixel_sharing.csv
#   P4_cluster_summary.csv
#   P4_gse_norm_check.csv
#   P4_gse_version_check.csv    (when CFG$gse_path is set) master_raw vs current GSE per year
#   P4_s2_availability.csv      AGC samples with a clear Sentinel-2 pixel, by year and ecoregion
#
# v3: STEPS 1-7 are those of P4v2, so P4_samples.csv keeps its rows and row order
# (the GEE training asset and every GEE export are joined on that order). STEP 8
# adds the Sentinel-2 columns at the end; STEP 9 the optional measured AGC.
# =============================================================================

source("R/P4_00_config_functions.R")

# =============================================================================
# STEP 1: Read master_raw and keep complete predictors
# =============================================================================
raw <- read_master()
cat(sprintf("Rows read from master_raw (%s): %d\n",
            paste(range(CFG$years_train), collapse = "-"), nrow(raw)))
cat(sprintf("Predictors: %d GSE bands + %s\n", length(CFG$gse_vars),
            if (length(CFG$env_vars) == 0) "no environmental layers" else paste(CFG$env_vars, collapse = ", ")))

# ---- Cover above 100 % (31 Oct 2026): capped at CFG$cover_cap -----------------
over_cap <- !is.na(raw$tSPC) & raw$tSPC > CFG$cover_cap
if (any(over_cap)) {
  cat(sprintf("tSPC > %g %%: %d samples (max %.1f) -> set to %g\n", CFG$cover_cap,
              sum(over_cap), max(raw$tSPC[over_cap]), CFG$cover_cap))
  print(as.data.frame(dplyr::count(raw[over_cap, ], year, name = "n_capped")), row.names = FALSE)
  raw$tSPC[over_cap] <- CFG$cover_cap
}

# =============================================================================
# STEP 2: Flag AGC and cover samples
# AGC: valid point estimate and 95% interval (one record per sample, no draws)
# Cover: tSPC > 0 ("present") or any tSPC ("all"), see CFG$cover_train_set
# =============================================================================
raw <- raw |>
  dplyr::mutate(
    is_agc = !is.na(AGC_pred) & AGC_pred > 0 & !is.na(AGC_low) & !is.na(AGC_up) &
             AGC_low >= 0 & AGC_up > AGC_low,
    is_cov = if (CFG$cover_train_set == "present") !is.na(tSPC) & tSPC > 0 else !is.na(tSPC),
    sd_ref = ifelse(is_agc, (AGC_up - AGC_low) / 3.92, NA_real_)
  ) |>
  dplyr::filter(is_agc | is_cov)

# ---- GSE from the current dataset version (CFG$gse_path, from P4_G00) -----------
# Samples not in the file get no GSE (they are dropped below), so every training
# vector comes from the same dataset version as the maps.
if (!is.na(CFG$gse_path) && !file.exists(CFG$gse_path))
  warning(sprintf("CFG$gse_path is set but %s was not found: master_raw GSE values are used.",
                  CFG$gse_path))
if (!is.na(CFG$gse_path) && file.exists(CFG$gse_path)) {
  gnew <- readr::read_csv(CFG$gse_path, show_col_types = FALSE,
                          col_types = readr::cols(gee_id = "c", .default = readr::col_double()))
  gnew <- gnew[, c("gee_id", "year", CFG$gse_vars)]
  for (v in CFG$gse_vars) gnew[[v]][!is.na(gnew[[v]]) & gnew[[v]] <= -9999] <- NA_real_
  key_new <- paste(gnew$gee_id, gnew$year)
  if (anyDuplicated(key_new) > 0) stop("Duplicated gee_id x year in ", basename(CFG$gse_path))
  m  <- match(paste(raw$gee_id, raw$year), key_new)
  ok <- !is.na(m)
  old <- as.matrix(raw[ok, CFG$gse_vars])
  new <- as.matrix(gnew[m[ok], CFG$gse_vars])
  gse_cmp <- tibble::tibble(year = raw$year[ok],
                            max_abs_diff = apply(abs(old - new), 1, max),
                            dot = rowSums(old * new))
  cat(sprintf("GSE from %s: %d of %d samples matched; %d not in the file (no GSE -> dropped)\n",
              basename(CFG$gse_path), sum(ok), nrow(raw), sum(!ok)))
  gse_year <- gse_cmp |>
    dplyr::group_by(year) |>
    dplyr::summarise(n = dplyr::n(), changed = sum(max_abs_diff > 1e-4, na.rm = TRUE),
                     min_dot = suppressWarnings(min(dot, na.rm = TRUE)),
                     median_dot_changed = median(dot[max_abs_diff > 1e-4], na.rm = TRUE),
                     .groups = "drop")
  cat("master_raw vs current GSE (changed = max |diff| > 1e-4):\n")
  print(as.data.frame(gse_year), row.names = FALSE)
  readr::write_csv(gse_year, out_file("P4_gse_version_check.csv"))
  raw[ok, CFG$gse_vars]  <- gnew[m[ok], CFG$gse_vars]
  raw[!ok, CFG$gse_vars] <- NA_real_
}

# ---- Optional exclusion list (documented data cleaning, e.g. points on land) ----
if (!is.na(CFG$exclude_ids_path) && file.exists(CFG$exclude_ids_path)) {
  excl <- readr::read_csv(CFG$exclude_ids_path, show_col_types = FALSE,
                          col_types = readr::cols(gee_id = "c", .default = readr::col_guess()))
  drop <- raw$gee_id %in% excl$gee_id
  cat(sprintf("Excluded by %s: %d samples (%d AGC, %d cover)\n", basename(CFG$exclude_ids_path),
              sum(drop), sum(drop & raw$is_agc), sum(drop & raw$is_cov)))
  print(as.data.frame(dplyr::count(raw[drop, ], year, name = "n_excluded")), row.names = FALSE)
  raw <- raw[!drop, ]
}

# ---- Diagnostic: which predictor would remove samples --------------------
raw$complete <- !is.na(raw$xcoord) & !is.na(raw$ycoord) &
  stats::complete.cases(as.data.frame(raw[, CFG$predictors]))
na_tbl <- tibble::tibble(
  check = c("coordinates", "any GSE band", CFG$env_vars, "any predictor or coordinate"),
  agc_samples_lost = c(
    sum(raw$is_agc & (is.na(raw$xcoord) | is.na(raw$ycoord))),
    sum(raw$is_agc & !stats::complete.cases(as.data.frame(raw[, CFG$gse_vars]))),
    vapply(CFG$env_vars, function(v) sum(raw$is_agc & is.na(raw[[v]])), numeric(1)),
    sum(raw$is_agc & !raw$complete)),
  cover_samples_lost = c(
    sum(raw$is_cov & (is.na(raw$xcoord) | is.na(raw$ycoord))),
    sum(raw$is_cov & !stats::complete.cases(as.data.frame(raw[, CFG$gse_vars]))),
    vapply(CFG$env_vars, function(v) sum(raw$is_cov & is.na(raw[[v]])), numeric(1)),
    sum(raw$is_cov & !raw$complete))
)
cat("\nSamples that missing predictors would remove:\n")
print(as.data.frame(na_tbl), row.names = FALSE)
readr::write_csv(na_tbl, out_file("P4_predictor_na.csv"))

lost_year <- raw |>
  dplyr::group_by(year) |>
  dplyr::summarise(agc_valid = sum(is_agc), agc_lost = sum(is_agc & !complete),
                   cov_valid = sum(is_cov), cov_lost = sum(is_cov & !complete), .groups = "drop")
print(as.data.frame(lost_year), row.names = FALSE)
if (sum(lost_year$agc_lost) > 0)
  cat("WARNING: AGC samples are lost to missing predictors. Check P4_predictor_na.csv before going on.\n")

samples <- raw |>
  dplyr::filter(complete) |>
  dplyr::select(gee_id, year, xcoord, ycoord, PA, tSPC, AGC_pred, AGC_low, AGC_up, sd_ref,
                is_agc, is_cov, dplyr::all_of(CFG$predictors), dplyr::any_of(CFG$diag_cols))

cat(sprintf("AGC samples: %d | cover samples: %d (%s) | union: %d\n",
            sum(samples$is_agc), sum(samples$is_cov), CFG$cover_train_set, nrow(samples)))
print(samples |> dplyr::group_by(year) |>
        dplyr::summarise(n_agc = sum(is_agc), n_cov = sum(is_cov),
                         mean_agc = mean(AGC_pred[is_agc]), .groups = "drop"))

# =============================================================================
# STEP 3: GSE vector length check (unit length => dot product = cosine similarity)
# =============================================================================
gse_norm <- sqrt(rowSums(as.matrix(samples[, CFG$gse_vars])^2))
norm_tbl <- tibble::tibble(min = min(gse_norm), q01 = quantile(gse_norm, 0.01),
                           median = median(gse_norm), q99 = quantile(gse_norm, 0.99),
                           max = max(gse_norm))
print(norm_tbl)
if (norm_tbl$min < 0.9 || norm_tbl$max > 1.1)
  stop("GSE vector lengths outside 0.9-1.1: check for fill codes (-9999) in the GSE columns.")
cat(if (all(abs(norm_tbl[, c("q01", "q99")] - 1) < 0.01))
  "GSE vectors are unit length: use the dot product for the inter-annual similarity test.\n" else
  "GSE vectors are NOT unit length: use explicit cosine similarity in GEE.\n")
readr::write_csv(norm_tbl, out_file("P4_gse_norm_check.csv"))

# =============================================================================
# STEP 4: Pixel groups (identical GSE values in the same year)
# =============================================================================
pk <- pixel_key(samples)
samples$g_pixel <- paste0("p", match(pk, unique(pk)))

pix_share <- samples |>
  dplyr::filter(is_agc) |>
  dplyr::add_count(g_pixel, name = "n_in_pixel") |>
  dplyr::summarise(n_agc = dplyr::n(), n_pixels = dplyr::n_distinct(g_pixel),
                   max_per_pixel = max(n_in_pixel),
                   pct_sharing = 100 * mean(n_in_pixel > 1))
print(pix_share)
readr::write_csv(pix_share, out_file("P4_pixel_sharing.csv"))

# =============================================================================
# STEP 5: Spatial blocks (100 m and 1 km)
# =============================================================================
for (nm in names(CFG$block_sizes)) {
  samples[[paste0("g_", nm)]] <- paste0(nm, "_", block_id(samples$xcoord, samples$ycoord,
                                                           CFG$block_sizes[[nm]]))
}

# =============================================================================
# STEP 6: DBSCAN clusters on AGC samples (eps 5 km), the N_eff unit
# Noise points = their own group. Cover-only samples join the group of the
# nearest AGC sample within eps, otherwise their own group.
# =============================================================================
agc_idx <- which(samples$is_agc)
xyz_agc <- lonlat_to_xyz(samples$xcoord[agc_idx], samples$ycoord[agc_idx])
db <- dbscan::dbscan(xyz_agc, eps = CFG$dbscan_eps, minPts = CFG$dbscan_minPts)

samples$dbscan_cluster <- NA_integer_
samples$dbscan_cluster[agc_idx] <- db$cluster
g_agc <- ifelse(db$cluster > 0, paste0("c", db$cluster), paste0("n", samples$gee_id[agc_idx]))
samples$g_cluster <- NA_character_
samples$g_cluster[agc_idx] <- g_agc

oth <- which(!samples$is_agc)
if (length(oth) > 0) {
  nn <- FNN::get.knnx(xyz_agc, lonlat_to_xyz(samples$xcoord[oth], samples$ycoord[oth]), k = 1)
  d  <- chord_to_arc(nn$nn.dist[, 1])
  samples$g_cluster[oth] <- ifelse(d <= CFG$dbscan_eps, g_agc[nn$nn.index[, 1]],
                                   paste0("o", samples$gee_id[oth]))
}

cl_sizes <- table(db$cluster[db$cluster > 0])
n_clustered <- sum(cl_sizes)
n0 <- if (length(cl_sizes) > 1) (n_clustered - sum(cl_sizes^2) / n_clustered) / (length(cl_sizes) - 1) else NA
cl_summary <- tibble::tibble(
  eps_m = CFG$dbscan_eps, minPts = CFG$dbscan_minPts,
  n_agc = length(agc_idx), n_clusters = length(cl_sizes), n_clustered = n_clustered,
  n_noise = sum(db$cluster == 0), adjusted_mean_cluster_size = n0,
  largest_cluster = max(cl_sizes)
)
print(cl_summary)
cat("Thesis version (4,728 AGC samples): 33 clusters, 4,250 clustered, 478 noise, adjusted mean size 108.0\n")
readr::write_csv(cl_summary, out_file("P4_cluster_summary.csv"))

# =============================================================================
# STEP 7: MEOW ecoregion (LOLO unit)
# =============================================================================
eco <- join_ecoregion(samples)
samples$ecoregion <- unname(eco[samples$gee_id])

counts_eco_year <- samples |>
  dplyr::filter(is_agc) |>
  dplyr::count(ecoregion, year, name = "n_agc") |>
  tidyr::pivot_wider(names_from = year, values_from = n_agc, values_fill = 0,
                     names_sort = TRUE) |>
  dplyr::mutate(total = rowSums(dplyr::across(-ecoregion))) |>
  dplyr::arrange(dplyr::desc(total))
print(as.data.frame(counts_eco_year), row.names = FALSE)
readr::write_csv(counts_eco_year, out_file("P4_counts_ecoregion_year.csv"))

# ---- Diagnostic: AGC samples without ACA depth or distance to land ----------
# These stay in the GSE-only model. A sample without ACA depth lies where the ACA
# bathymetry has no value (no coverage, or a point on land); the GEE extent mask
# (depth 0-5 m) also drops such pixels.
if (all(c("depth", "distToLand") %in% names(samples))) {
  env_diag <- samples |>
    dplyr::filter(is_agc) |>
    dplyr::group_by(ecoregion) |>
    dplyr::summarise(n_agc = dplyr::n(), no_depth = sum(is.na(depth)),
                     no_distToLand = sum(is.na(distToLand)), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(no_depth))
  cat("\nAGC samples without depth or distToLand (kept, not used as predictors):\n")
  print(as.data.frame(env_diag), row.names = FALSE)
  nod <- samples |> dplyr::filter(is_agc, is.na(depth) | is.na(distToLand))
  readr::write_csv(nod |> dplyr::transmute(gee_id, year, xcoord, ycoord, ecoregion, AGC_pred, tSPC,
                                           no_depth = is.na(depth), no_distToLand = is.na(distToLand)),
                   out_file("P4_agc_without_env.csv"))
  cat(sprintf("Saved: P4_agc_without_env.csv (%d AGC samples; open in QGIS to find points on land)\n", nrow(nod)))
}

eco_small <- counts_eco_year$ecoregion[counts_eco_year$total < CFG$min_n_eco]
if (length(eco_small) > 0)
  cat(sprintf("Ecoregions below n = %d (pooled metrics only): %s\n",
              CFG$min_n_eco, paste(eco_small, collapse = ", ")))

# =============================================================================
# STEP 8: Sentinel-2 predictors at the AGC samples (P4_G11), survey-year composite
# has_s2 = AGC sample with all percentile features and >= CFG$s2_min_obs clear
# observations in its survey year. Only these samples enter the Sentinel-2 comparison.
# =============================================================================
s2_files <- list.files(CFG$s2_dir, pattern = CFG$s2_pattern, full.names = TRUE)
if (length(s2_files) == 0) {
  cat(sprintf("\nSentinel-2: no %s files in %s (run P4_G11): has_s2 = FALSE for all samples\n",
              CFG$s2_pattern, CFG$s2_dir))
  samples[c(CFG$s2_vars, "S2_nobs")] <- NA_real_
  samples$has_s2 <- FALSE
} else {
  s2 <- purrr::map_dfr(s2_files, readr::read_csv, show_col_types = FALSE,
                       col_types = readr::cols(gee_id = "c", .default = readr::col_double()))
  miss <- setdiff(c("gee_id", "year", CFG$s2_vars, "S2_nobs"), names(s2))
  if (length(miss) > 0) stop("Sentinel-2 files lack columns: ", paste(miss, collapse = ", "))
  s2 <- s2[, c("gee_id", "year", CFG$s2_vars, "S2_nobs")]
  for (v in CFG$s2_vars) s2[[v]][!is.na(s2[[v]]) & s2[[v]] <= -9999] <- NA_real_
  s2$S2_nobs[is.na(s2$S2_nobs) | s2$S2_nobs < 0] <- 0
  key_s2 <- paste(s2$gee_id, s2$year)
  if (anyDuplicated(key_s2) > 0) stop("Duplicated gee_id x year in the Sentinel-2 files")
  m <- match(paste(samples$gee_id, samples$year), key_s2)
  samples[c(CFG$s2_vars, "S2_nobs")] <- s2[m, c(CFG$s2_vars, "S2_nobs")]
  samples$has_s2 <- samples$is_agc & !is.na(m) & dplyr::coalesce(samples$S2_nobs, 0) >= CFG$s2_min_obs &
    stats::complete.cases(as.data.frame(samples[, CFG$s2_vars]))
  if (!is.null(CFG$s2_years)) {
    samples$has_s2 <- samples$has_s2 & samples$year %in% CFG$s2_years
    cat(sprintf("Sentinel-2 comparison limited to the survey years %s (CFG$s2_years)\n",
                paste(CFG$s2_years, collapse = ", ")))
  }
  cat(sprintf("\nSentinel-2 (%d files): %d of %d AGC samples in the files, %d with a clear pixel (has_s2)\n",
              length(s2_files), sum(samples$is_agc & !is.na(m)), sum(samples$is_agc), sum(samples$has_s2)))
}

s2_avail <- dplyr::bind_rows(
  samples |> dplyr::filter(is_agc) |>
    dplyr::group_by(year = as.character(year), ecoregion) |>
    dplyr::summarise(n_agc = dplyr::n(), n_has_s2 = sum(has_s2),
                     median_nobs = suppressWarnings(stats::median(S2_nobs, na.rm = TRUE)), .groups = "drop"),
  samples |> dplyr::filter(is_agc) |>
    dplyr::group_by(year = as.character(year)) |>
    dplyr::summarise(ecoregion = "All", n_agc = dplyr::n(), n_has_s2 = sum(has_s2),
                     median_nobs = suppressWarnings(stats::median(S2_nobs, na.rm = TRUE)), .groups = "drop"),
  samples |> dplyr::filter(is_agc) |>
    dplyr::summarise(year = "All", ecoregion = "All", n_agc = dplyr::n(), n_has_s2 = sum(has_s2),
                     median_nobs = suppressWarnings(stats::median(S2_nobs, na.rm = TRUE)))
) |>
  dplyr::mutate(pct_has_s2 = 100 * n_has_s2 / n_agc)
cat("AGC samples with a clear Sentinel-2 pixel, by survey year:\n")
print_tbl(s2_avail |> dplyr::filter(ecoregion == "All"), 1)
readr::write_csv(s2_avail, out_file("P4_s2_availability.csv"))

# =============================================================================
# STEP 9 (optional): AGC measured directly in the field (CFG$measured_path)
# =============================================================================
meas_path <- if (is.na(CFG$measured_path)) out_file("P4_measured_agc.csv") else CFG$measured_path
if (isTRUE(CFG$use_measured) && file.exists(meas_path)) {
  meas <- readr::read_csv(meas_path, show_col_types = FALSE,
                          col_types = readr::cols_only(gee_id = "c", year = "i", AGC_measured = "d"))
  m <- match(paste(samples$gee_id, samples$year), paste(meas$gee_id, meas$year))
  samples$AGC_measured <- meas$AGC_measured[m]
  cat(sprintf("Measured AGC: %d AGC samples matched\n", sum(samples$is_agc & !is.na(samples$AGC_measured))))
}

# =============================================================================
# STEP 10: Save
# =============================================================================
readr::write_csv(samples, out_file("P4_samples.csv"))
cat(sprintf("Saved: P4_samples.csv (%d rows x %d cols)\n", nrow(samples), ncol(samples)))

# The GEE upload table (with GSE values and a row id) is written by P4_06_gee_prepare.R

# =============================================================================
# STEP 11: Plots (printed only)
# =============================================================================
comp <- samples |>
  dplyr::filter(is_agc) |>
  dplyr::count(year, ecoregion, name = "n")
print(
  ggplot(comp, aes(x = factor(year), y = n, fill = ecoregion)) +
    geom_col(position = "fill", colour = "white", linewidth = 0.2) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "AGC samples by survey year and MEOW ecoregion",
         x = "Survey year", y = "Share of AGC samples", fill = "Ecoregion") +
    theme_minimal(base_size = 11)
)

map_df <- samples |> dplyr::filter(is_agc)
p_map <- ggplot(map_df, aes(x = xcoord, y = ycoord, colour = ecoregion)) +
  geom_point(size = 0.8, alpha = 0.6) +
  coord_quickmap() +
  labs(title = sprintf("AGC samples (n = %d) by MEOW ecoregion", nrow(map_df)),
       x = "Longitude", y = "Latitude", colour = "Ecoregion") +
  theme_bw(base_size = 10)
if (requireNamespace("maps", quietly = TRUE)) {
  p_map <- p_map + geom_polygon(data = ggplot2::map_data("world", region = "Indonesia"),
                                aes(x = long, y = lat, group = group), inherit.aes = FALSE,
                                fill = NA, colour = "grey50", linewidth = 0.2)
}
print(p_map)
ggsave(out_file("P4_fig_samples_map.png"), p_map, width = 10, height = 5, dpi = 300)

cat("\nP4_01_build_training.R COMPLETE\n")
