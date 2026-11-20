# =============================================================================
# P4_01a_link_measured_agc.R
# Paper 4, design v3: give the field-measured AGB and AGC of the input data of the
# cover-AGB model the gee_id of the P4 samples.
#
# Chain of the sample IDs:
#   predict_tAGB2.csv (cover-AGB model input; composition_id, X, Y, year, information,
#   total_AGB, total_AGC)
#     -> cover-AGB RF -> RF_total_AGB_AGC_summary_with_ratio_v4.csv
#     -> + total_SPC, duplicates of (X, Y) removed (largest total_SPC kept)
#        -> RF_AGC_summary_v4_clean_for_QGIS(_merged).csv, where gee_id was added
#     -> master_raw_YYYY.csv (01_build_master_raw.R) keeps gee_id and compositio
#        (composition_id, field name cut to 10 characters by the shapefile)
# So composition_id = compositio links each measured quadrat to its gee_id.
# Year and coordinates are checked for every link.
#
# Rules
#   - measured rows: information == "actual" with total_AGC, survey years
#     CFG$years_train (GSE starts in 2017; earlier surveys cannot enter P4)
#   - link "id"   : composition_id found in master_raw, same year, <= CFG$link_max_dist_m
#   - link "id, year +-1": the same, survey years differ by <= CFG$link_year_tol (the
#                   P4 sample keeps its master_raw year)
#   - link "coord": composition_id not in master_raw, but a master_raw sample of the
#                   same year lies within CFG$link_max_dist_m and has no compositio
#   - not linked  : same position but another compositio (the quadrat was removed by
#                   the (X, Y) de-duplication), or no master_raw sample nearby
#   - CFG$measured_methods limits which `method` values of the cover-AGB model input count as measured
#     (NULL = every "actual" row; decide from P4_measured_link_report.csv)
#
# Reads : CFG$p1_input_path (predict_tAGB2.csv), data/master_raw/master_raw_YYYY.csv,
#         P4_samples.csv (optional: marks links that are P4 AGC samples)
# Writes: P4_measured_agc.csv          gee_id, year, AGC_measured, AGB_measured, ...
#                                      (read by P4_01 STEP 9)
#         P4_measured_link_report.csv  per location x year x method: rows and link result
#         P4_measured_link_detail.csv  every "actual" row with its link result
# Run   : before P4_01 (or after it, then rerun P4_01 and P4_04).
# OPTIONAL DIAGNOSTIC (15 Nov 2026): not used for the paper results. The "actual"
# values of the cover-AGB model input are not an independent reference (see CFG$use_measured).
# =============================================================================

source("R/P4_00_config_functions.R")

haversine_m <- function(lon1, lat1, lon2, lat2) {
  to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad; dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat / 2)^2 + cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  2 * EARTH_R * asin(pmin(1, sqrt(a)))
}

# =============================================================================
# STEP 1: input data of the cover-AGB model, rows with field-measured biomass
# =============================================================================
if (!file.exists(CFG$p1_input_path)) stop("Not found: ", CFG$p1_input_path, " (CFG$p1_input_path)")
p1 <- readr::read_csv(
  CFG$p1_input_path, na = c("", "NA", "-", "NaN"), show_col_types = FALSE,
  col_types = readr::cols(location = "c", region = "c", ownner = "c", method = "c",
                          year = "i", X = "d", Y = "d", original_id = "c", information = "c",
                          sample_id = "c", composition_id = "c", .default = readr::col_double()))
need <- c("composition_id", "year", "X", "Y", "information", "total_AGB", "total_AGC")
if (length(setdiff(need, names(p1))) > 0)
  stop("predict_tAGB2.csv lacks: ", paste(setdiff(need, names(p1)), collapse = ", "))
if (anyDuplicated(p1$composition_id[!is.na(p1$composition_id)]) > 0)
  warning("composition_id is not unique in predict_tAGB2.csv: links use the first row of each id")

act <- p1 |>
  dplyr::filter(tolower(trimws(information)) == "actual", is.finite(total_AGC)) |>
  dplyr::mutate(in_period = year %in% CFG$years_train,
                method_ok = if (is.null(CFG$measured_methods)) TRUE else method %in% CFG$measured_methods)
cat(sprintf("cover-AGB model input: %d rows, %d 'actual' with total_AGC, %d of them in %d-%d\n",
            nrow(p1), nrow(act), sum(act$in_period), min(CFG$years_train), max(CFG$years_train)))

# =============================================================================
# STEP 2: master_raw IDs (gee_id, compositio, coordinates, AGC_pred)
# =============================================================================
mr <- purrr::map_dfr(CFG$years_train, function(yr) {
  p <- file.path(CFG$result_dir, paste0("master_raw_", yr, ".csv"))
  if (!file.exists(p)) return(NULL)
  d <- readr::read_csv(p, show_col_types = FALSE, guess_max = 100000,
                       col_types = readr::cols(gee_id = "c", compositio = "c", .default = readr::col_guess()))
  if (!"compositio" %in% names(d)) d$compositio <- NA_character_
  tibble::tibble(gee_id = d$gee_id, mr_year = as.integer(yr), compositio = d$compositio,
                 xcoord = as.numeric(d$xcoord), ycoord = as.numeric(d$ycoord),
                 AGB_pred = if ("AGB_pred" %in% names(d)) suppressWarnings(as.numeric(d$AGB_pred)) else NA_real_,
                 AGC_pred = suppressWarnings(as.numeric(d$AGC_pred)))
})
cat(sprintf("master_raw: %d rows, %d with compositio\n", nrow(mr), sum(!is.na(mr$compositio))))
if (all(is.na(mr$compositio)))
  stop("master_raw has no compositio values: only coordinate links would be possible. ",
       "Check the label file of 01_build_master_raw.R.")

# =============================================================================
# STEP 3: Link by composition_id, then by position for the rest
# =============================================================================
by_id <- act |>
  dplyr::filter(in_period) |>
  dplyr::left_join(mr |> dplyr::filter(!is.na(compositio)) |> dplyr::distinct(compositio, .keep_all = TRUE),
                   by = c("composition_id" = "compositio")) |>
  dplyr::mutate(dist_m = haversine_m(X, Y, xcoord, ycoord),
                link = dplyr::case_when(
                  is.na(gee_id) ~ NA_character_,
                  mr_year == year & dist_m <= CFG$link_max_dist_m ~ "id",
                  abs(mr_year - year) <= CFG$link_year_tol & dist_m <= CFG$link_max_dist_m ~ "id, year +-1",
                  mr_year != year ~ "id, other year (not linked)",
                  TRUE ~ sprintf("id, %.0f m apart (not linked)", dist_m)))

# position search for rows whose composition_id is not in master_raw
rest <- by_id |> dplyr::filter(is.na(gee_id)) |> dplyr::select(dplyr::all_of(names(act)))
pos <- purrr::map_dfr(seq_len(nrow(rest)), function(i) {
  r <- rest[i, ]
  cand <- mr |> dplyr::filter(mr_year == r$year)
  if (nrow(cand) == 0) return(dplyr::mutate(r, link = "no master_raw sample in that year"))
  d <- haversine_m(r$X, r$Y, cand$xcoord, cand$ycoord)
  j <- which.min(d)
  if (d[j] > CFG$link_max_dist_m)
    return(dplyr::mutate(r, link = "not in master_raw (no sample nearby)", dist_m = d[j]))
  c1 <- cand[j, ]
  if (!is.na(c1$compositio))
    return(dplyr::mutate(r, link = "same position, other quadrat kept (X, Y de-duplication)",
                         dist_m = d[j], other_compositio = c1$compositio))
  dplyr::bind_cols(r, c1 |> dplyr::select(gee_id, mr_year, xcoord, ycoord, AGB_pred, AGC_pred)) |>
    dplyr::mutate(dist_m = d[j], link = "coord")
})

detail <- dplyr::bind_rows(
  by_id |> dplyr::filter(!is.na(gee_id)),
  pos,
  act |> dplyr::filter(!in_period) |> dplyr::mutate(link = "survey year outside P4 years")
) |>
  dplyr::mutate(linked = link %in% c("id", "id, year +-1", "coord") & method_ok)

if (file.exists(out_file("P4_samples.csv"))) {
  s <- read_samples()
  detail$in_P4_agc <- detail$gee_id %in% s$gee_id[s$is_agc]
} else {
  detail$in_P4_agc <- NA
}

# one link per gee_id (two measured quadrats on one gee_id would be ambiguous)
dup <- detail |> dplyr::filter(linked) |> dplyr::count(gee_id) |> dplyr::filter(n > 1)
if (nrow(dup) > 0) {
  detail$link[detail$linked & detail$gee_id %in% dup$gee_id] <- "two measured quadrats on one gee_id"
  detail$linked[detail$gee_id %in% dup$gee_id] <- FALSE
}
readr::write_csv(detail, out_file("P4_measured_link_detail.csv"))

# =============================================================================
# STEP 4: Report and output
# =============================================================================
report <- detail |>
  dplyr::group_by(location, year, method) |>
  dplyr::summarise(n_actual = dplyr::n(), n_linked_id = sum(linked & link %in% c("id", "id, year +-1")),
                   n_linked_coord = sum(linked & link == "coord"),
                   n_in_P4_agc = sum(linked & in_P4_agc %in% TRUE),
                   not_linked = paste(unique(link[!linked]), collapse = "; "),
                   mean_AGC_measured = mean(total_AGC),
                   # same quadrats: cover-AGB model output in master_raw vs the input value
                   mean_AGC_pred_master = mean(AGC_pred[linked], na.rm = TRUE),
                   ratio_AGC_pred_measured = stats::median((AGC_pred / total_AGC)[linked], na.rm = TRUE),
                   mean_AGB_measured = mean(total_AGB, na.rm = TRUE),
                   mean_AGB_pred_master = mean(AGB_pred[linked], na.rm = TRUE),
                   ratio_AGB_pred_measured = stats::median((AGB_pred / total_AGB)[linked], na.rm = TRUE),
                   ratio_AGC_AGB_measured = stats::median(total_AGC / total_AGB, na.rm = TRUE),
                   ratio_AGC_AGB_master = stats::median((AGC_pred / AGB_pred)[linked], na.rm = TRUE),
                   .groups = "drop") |>
  dplyr::arrange(year, location)
cat("\nLink result per location x year x method:\n")
print_tbl(report, 1)
readr::write_csv(report, out_file("P4_measured_link_report.csv"))

meas <- detail |>
  dplyr::filter(linked) |>
  dplyr::transmute(gee_id, year_p1 = year, year = mr_year, AGC_measured = total_AGC,
                   AGB_measured = total_AGB, composition_id, location, method, link, dist_m,
                   AGB_pred_master = AGB_pred, AGC_pred_master = AGC_pred, in_P4_agc)
readr::write_csv(meas, out_file("P4_measured_agc.csv"))

# Unit and plausibility check: AGB_pred and AGC_pred of master_raw vs the measured AGB
# and AGC on the same quadrat (report columns ratio_*). For these quadrats the cover-AGB
# model was fitted with them, so the ratios should be near 1; a ratio far from 1 means
# the two values are not the same quantity (unit, quadrat area or carbon fraction).
# ratio_AGC_AGB_measured vs ratio_AGC_AGB_master shows whether the carbon fraction differs.
if (nrow(meas) >= 3) {
  ok <- is.finite(meas$AGC_pred_master)
  cat(sprintf(paste0("\nLinked: %d quadrats (%d by id, %d by position), %d of them P4 AGC samples\n",
                     "AGC_pred / AGC_measured: median %.2f | r = %.2f (fit check, not validation)\n"),
              nrow(meas), sum(meas$link %in% c("id", "id, year +-1")), sum(meas$link == "coord"),
              sum(meas$in_P4_agc %in% TRUE),
              median(meas$AGC_pred_master[ok] / meas$AGC_measured[ok]),
              cor(meas$AGC_pred_master[ok], meas$AGC_measured[ok])))
  r <- median(meas$AGC_pred_master[ok] / meas$AGC_measured[ok])
  if (r < 0.5 || r > 2) warning("AGC_pred and total_AGC differ by more than a factor of 2: check the units.")
}
cat(sprintf("Saved: P4_measured_agc.csv, P4_measured_link_report.csv, P4_measured_link_detail.csv\n"))
cat("Next: P4_01 (STEP 9 reads P4_measured_agc.csv), then P4_04 (STEP 8).\n")
cat("\nP4_01a_link_measured_agc.R COMPLETE\n")
