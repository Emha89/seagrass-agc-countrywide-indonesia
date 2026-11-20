# =============================================================================
# P4_01b_trace_actual_agb.R
# Paper 4, design v3: trace the actual AGB of each quadrat through the three
# files of the cover-AGB model, before the measured AGB/AGC are used in P4.
#
#   1. predict_tAGB2.csv                 "actual" rows: total_AGB (Excel: per species
#                                        shoot biomass x shoot density, summed, per m2)
#   2. data_long_format.rds              training table of the cover-AGB model
#                                        (01_model_dev.R): AGB per species, cover 0-1
#   3. RF_total_AGB_AGC_summary_with_ratio_v4.csv
#                                        apply output: total_AGB_pred = sum of the
#                                        species AGB predicted from species cover
#   4. master_raw_YYYY.csv               AGB_pred and AGC_pred (predictions, from file 3)
#
# Columns: ACTUAL = total_AGB, total_AGC (file 1) and AGB (file 2);
#          PREDICTED = total_AGB_pred, total_AGC_pred (file 3) = AGB_pred, AGC_pred (file 4).
#
# Questions per location x year x method:
#   - Is the quadrat in the training table (sample_id)?
#   - Is total_AGB of the input file the sum of the species AGB of the training
#     table (ratio_input_train ~ 1)? Is the total cover the same?
#   - How far is the model output from the actual AGB (ratio_pred_actual)? For
#     quadrats used in training a random forest should be close to 1.
# Writes: P4_actual_agb_trace.csv        per location x year x method
#         P4_actual_agb_trace_detail.csv per quadrat
# OPTIONAL DIAGNOSTIC (15 Nov 2026): not used for the paper results. The "actual"
# values of the cover-AGB model input are not an independent reference (see CFG$use_measured).
# =============================================================================

source("R/P4_00_config_functions.R")

# Required: the input file and the apply output. The training table is optional:
# data_long_format.rds (01_model_dev.R) or, without it, total_agb_actual_v4.csv
# (sample_id, total_AGB; also written by 01_model_dev.R). Both are in the output folder of the cover-AGB model.
for (p in c(CFG$p1_input_path, CFG$p1_summary_path))
  if (!file.exists(p)) stop("Not found: ", p, " (copy it from the output folder of the cover-AGB model into data/cover_agb/)")

p1 <- readr::read_csv(
  CFG$p1_input_path, na = c("", "NA", "-", "NaN"), show_col_types = FALSE,
  col_types = readr::cols(location = "c", region = "c", ownner = "c", method = "c",
                          year = "i", X = "d", Y = "d", original_id = "c", information = "c",
                          sample_id = "c", composition_id = "c", .default = readr::col_double()))
act <- p1 |> dplyr::filter(tolower(trimws(information)) == "actual")

# training table: one row per sample x species (or one row per sample as fallback)
train_total_path <- file.path(dirname(CFG$p1_train_path), "total_agb_actual_v4.csv")
if (file.exists(CFG$p1_train_path)) {
  dl <- readRDS(CFG$p1_train_path) |> dplyr::mutate(sample_id = as.character(sample_id))
  train <- dl |>
    dplyr::group_by(sample_id) |>
    dplyr::summarise(train_location = dplyr::first(as.character(location)),
                     n_species_train = dplyr::n(),
                     AGB_train = sum(AGB, na.rm = TRUE),
                     cover_train_pct = 100 * dplyr::first(coverTotal), .groups = "drop")
  cat(sprintf("Training table: %s (%d samples)\n", basename(CFG$p1_train_path), nrow(train)))
} else if (file.exists(train_total_path)) {
  train <- readr::read_csv(train_total_path, show_col_types = FALSE,
                           col_types = readr::cols(sample_id = "c", .default = readr::col_double())) |>
    dplyr::group_by(sample_id) |>
    dplyr::summarise(AGB_train = dplyr::first(total_AGB), .groups = "drop") |>
    dplyr::mutate(train_location = NA_character_, n_species_train = NA_integer_,
                  cover_train_pct = NA_real_)
  cat(sprintf("Training table: %s (%d samples; no cover, so cover_train_pct is NA)\n",
              basename(train_total_path), nrow(train)))
} else {
  cat(sprintf("No training table (%s or %s): the training columns are NA\n",
              CFG$p1_train_path, train_total_path))
  train <- tibble::tibble(sample_id = character(0), train_location = character(0),
                          n_species_train = integer(0), AGB_train = numeric(0),
                          cover_train_pct = numeric(0))
}

summ <- readr::read_csv(CFG$p1_summary_path, show_col_types = FALSE,
                        col_types = readr::cols(composition_id = "c", .default = readr::col_guess())) |>
  dplyr::distinct(composition_id, .keep_all = TRUE) |>
  dplyr::select(composition_id, total_AGB_pred, total_AGC_pred)

d <- act |>
  dplyr::select(location, year, method, sample_id, composition_id, total_SPC, total_AGB, total_AGC) |>
  dplyr::left_join(train, by = "sample_id") |>
  dplyr::left_join(summ, by = "composition_id") |>
  dplyr::mutate(in_training = !is.na(AGB_train),
                ratio_input_train = total_AGB / AGB_train,
                ratio_pred_actual = total_AGB_pred / total_AGB)
# predictions in master_raw = predictions in the apply output? (compositio = composition_id)
mr <- purrr::map_dfr(CFG$years_train, function(yr) {
  f <- file.path(CFG$result_dir, paste0("master_raw_", yr, ".csv"))
  if (!file.exists(f)) return(NULL)
  x <- readr::read_csv(f, show_col_types = FALSE, guess_max = 100000,
                       col_types = readr::cols(compositio = "c", .default = readr::col_guess()))
  if (!all(c("compositio", "AGB_pred", "AGC_pred") %in% names(x))) return(NULL)
  tibble::tibble(composition_id = x$compositio, AGB_pred_master = as.numeric(x$AGB_pred),
                 AGC_pred_master = as.numeric(x$AGC_pred))
})
if (nrow(mr) > 0) {
  mr <- mr |> dplyr::filter(!is.na(composition_id)) |> dplyr::distinct(composition_id, .keep_all = TRUE)
  d <- dplyr::left_join(d, mr, by = "composition_id")
} else {
  cat("master_raw has no compositio/AGB_pred columns: prediction check skipped\n")
  d$AGB_pred_master <- NA_real_; d$AGC_pred_master <- NA_real_
}
k <- is.finite(d$AGB_pred_master) & is.finite(d$total_AGB_pred)
if (any(k)) {
  cat(sprintf(paste0("Check of the prediction columns (%d actual quadrats in master_raw):
",
                     "  master_raw AGB_pred = apply total_AGB_pred within 0.01: %.1f%%
",
                     "  master_raw AGC_pred = apply total_AGC_pred within 0.01: %.1f%%
"),
              sum(k), 100 * mean(abs(d$AGB_pred_master[k] - d$total_AGB_pred[k]) < 0.01),
              100 * mean(abs(d$AGC_pred_master[k] - d$total_AGC_pred[k]) < 0.01, na.rm = TRUE)))
}
readr::write_csv(d, out_file("P4_actual_agb_trace_detail.csv"))

med <- function(x) suppressWarnings(stats::median(x, na.rm = TRUE))
trace <- d |>
  dplyr::group_by(location, year, method) |>
  dplyr::summarise(n_actual = dplyr::n(), n_in_training = sum(in_training),
                   AGB_actual = med(total_AGB), AGB_train = med(AGB_train),
                   ratio_input_train = med(ratio_input_train),
                   cover_input_pct = med(total_SPC), cover_train_pct = med(cover_train_pct),
                   AGB_pred = med(total_AGB_pred),
                   # before ratio_pred_actual: summarise() evaluates in order
                   ratio_pred_actual_in_training = med(ratio_pred_actual[in_training]),
                   ratio_pred_actual = med(ratio_pred_actual),
                   .groups = "drop") |>
  dplyr::arrange(year, location)
cat("\nActual AGB through the files of the cover-AGB model (medians; AGB in g DW m-2, cover in %):\n")
print_tbl(trace, 2)
readr::write_csv(trace, out_file("P4_actual_agb_trace.csv"))

cat("\nReading the table:\n",
    " n_in_training < n_actual       some actual quadrats did not train the cover-AGB model\n",
    " ratio_input_train far from 1   total_AGB of the input file is not the training AGB\n",
    " cover_input vs cover_train     the model was applied to another cover than it learned from\n",
    " ratio_pred_actual far from 1   the model output differs from the actual AGB at the same quadrat\n")
cat("\nP4_01b_trace_actual_agb.R COMPLETE\n")
