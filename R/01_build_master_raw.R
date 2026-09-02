# =============================================================================
# 01_build_master_raw.R
# Build master raw dataset -- one file per year
#
# Replaces: dataPrep_func_PA.R + dataPrep_composSP.R + dataPrep_cIndex.R
# (an earlier, multi-script version of this same data-prep step -- this
# single script now does all three jobs at once: GT+GEE join, species
# composition join, and carbon_index/morph3 derivation).
#
# Join logic:
#   GT label data  INNER JOIN  GEE image data    (on gee_id)
#   + sg_compos                                  LEFT JOIN  (on compositio)
#   + derived: carbon_index, morph3              computed from species cols
#
# Output per year:
#   master_raw_YYYY.csv
#   - All columns from all sources; NA where data absent (species, derived)
#   - dataPrep_func_AGB.R and dataPrep_func_AGC.R now read from this file
#
# Column groups in master_raw_YYYY:
#   ID/spatial  : gee_id, year, xcoord, ycoord, OBJECTID, location, compositio
#   GT labels   : PA, tSPC, AGB_pred, AGB_low, AGB_up, AGB_CIwidt,
#                 AGC_pred, AGC_low, AGC_up, AGC_CIwidt, sg_morpho
#   GEE spectral: depth, distToLand, rugosity, slope, elevation,
#                 mean_wave_period, sig_wave_height, GSE_A00..GSE_A63
#   Species     : Ea_SPC, Th_SPC, Cr_SPC, Cs_SPC, Si_SPC, Hu_SPC,
#                 Ho_SPC, Hp_SPC, Tc_SPC, Hm_SPC, Hs_SPC, Hd_SPC
#   Derived     : carbon_index  (0-1), morph3  (3-class factor)
#
# FIX APPLIED: removed a stray line, `df <- normalise_master_types(df)`,
# that referenced an undefined variable `df` (the script never creates a
# plain `df` -- only `gt_df`). This would have stopped the script with
# "object 'df' not found" on every run. Looks like a leftover from an
# earlier draft before the variable was renamed to gt_df.
#
# xcoord/ycoord are explicitly taken from the GT label file, not the GEE
# training file (the GEE file's own coordinate columns, if present, are
# dropped before the join) -- GT's coordinates are treated as
# authoritative.
#
# Data availability: this script's inputs are not included in this
# repository.
# =============================================================================

library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(ggplot2)
library(maps)
library(here)   # install.packages("here") if you don't have it yet


# Helper: normalise column types -- called after column selection
# so only known columns remain; no ambiguous metadata columns
normalise_master_types <- function(df) {
  num_cols <- intersect(
    c("PA", "tSPC",
      "AGB_pred", "AGB_low", "AGB_up", "AGB_CIwidt",
      "AGC_pred", "AGC_low", "AGC_up", "AGC_CIwidt",
      "carbon_index", "xcoord", "ycoord",
      "depth", "distToLand",
      "Ea_SPC","Th_SPC","Cr_SPC","Cs_SPC","Si_SPC","Hu_SPC",
      "Ho_SPC","Hp_SPC","Tc_SPC","Hm_SPC","Hs_SPC","Hd_SPC",
      paste0("GSE_A", sprintf("%02d", 0:63))),
    names(df))
  char_cols <- intersect(
    c("gee_id", "compositio", "sg_morpho", "morph3"),
    names(df))
  int_cols  <- intersect(c("year_gt"), names(df))
  
  for (col in num_cols)  df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  for (col in char_cols) df[[col]] <- as.character(df[[col]])
  for (col in int_cols)  df[[col]] <- suppressWarnings(as.integer(df[[col]]))
  df
}


# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
result_dir <- here("result")
label_path <- here("data", "RF_AGC_summary_v4_clean_for_QGIS_merged.csv")
species_path <- here("data", "sg_compos.csv")

gee_paths <- list(
  "2017" = here("data", "GSE_training_2017_CSV.csv"),
  "2018" = here("data", "GSE_training_2018_CSV.csv"),
  "2019" = here("data", "GSE_training_2019_CSV.csv"),
  "2020" = here("data", "GSE_training_2020_CSV.csv"),
  "2021" = here("data", "GSE_training_2021_CSV.csv"),
  "2022" = here("data", "GSE_training_2022_CSV.csv"),
  "2023" = here("data", "GSE_training_2023_CSV.csv"),
  "2024" = here("data", "GSE_training_2024_CSV.csv")
)
years <- names(gee_paths)
dir.create(result_dir, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Carbon index parameters (Hafiz data)
# -----------------------------------------------------------------------------
carbon_weights <- c(
  Ea_SPC = 0.829, Th_SPC = 0.321, Cr_SPC = 0.152, Cs_SPC = 0.190,
  Si_SPC = 0.053, Hu_SPC = 0.014, Ho_SPC = 0.017, Hp_SPC = 0.031,
  Tc_SPC = 0.108, Hm_SPC = 0.000, Hs_SPC = 0.000, Hd_SPC = 0.000
)
max_bonus    <- 1
morph_levels <- c("mixed_short_plus_mono_short", "mixed_long", "mono_Ea")
species_cols <- names(carbon_weights)

# =============================================================================
# STEP 1: Load GT label data (all years combined, deduped on gee_id)
# =============================================================================
cat("Loading GT label data...\n")
gt_df <- read_csv(label_path, show_col_types = FALSE) %>%
  mutate(year = as.integer(year)) %>%
  distinct(gee_id, .keep_all = TRUE)
gt_df <- normalise_master_types(gt_df)

if (!all(c("xcoord", "ycoord") %in% names(gt_df)))
  stop("xcoord / ycoord missing from GT label file.")

gt_valid <- gt_df %>% filter(!is.na(year))
cat(sprintf("GT rows loaded: %d | years: %s\n\n",
            nrow(gt_valid), paste(sort(unique(gt_valid$year)), collapse = ", ")))

# Check duplicates
dup_ids <- gt_valid %>% group_by(gee_id) %>% filter(n() > 1) %>% ungroup()
if (nrow(dup_ids) > 0) {
  cat(sprintf("WARNING: %d duplicate gee_id in GT data -- using first occurrence.\n", nrow(dup_ids)))
  write_csv(dup_ids, file.path(result_dir, "diagnostic_duplicate_gee_ids.csv"))
}

# =============================================================================
# STEP 2: Load species composition (sg_compos.csv)
# =============================================================================
cat("Loading species composition data...\n")
species_df <- read_csv(species_path, show_col_types = FALSE)
cat(sprintf("Species rows loaded: %d\n\n", nrow(species_df)))


# =============================================================================
# Define columns to retain in master_raw_YYYY
# Only columns used by at least one model as predictor, response, or join key
# =============================================================================
id_cols      <- c("gee_id", "year_gt", "xcoord", "ycoord", "compositio")
response_cols <- c("PA", "tSPC",
                   "AGB_pred", "AGB_low", "AGB_up", "AGB_CIwidt",
                   "AGC_pred", "AGC_low", "AGC_up", "AGC_CIwidt",
                   "sg_morpho", "morph3", "carbon_index")
keep_species_cols <- c("Ea_SPC","Th_SPC","Cr_SPC","Cs_SPC","Si_SPC","Hu_SPC",
                       "Ho_SPC","Hp_SPC","Tc_SPC","Hm_SPC","Hs_SPC","Hd_SPC")
keep_env_cols     <- c("depth", "distToLand")
keep_gse_cols     <- paste0("GSE_A", sprintf("%02d", 0:63))

master_keep_cols  <- c(id_cols, response_cols, keep_species_cols, keep_env_cols, keep_gse_cols)

# =============================================================================
# STEP 3: Build master_raw_YYYY per year
# =============================================================================
cat("Building master_raw_YYYY files...\n\n")

master_list <- list()

for (yr in years) {
  cat(sprintf("--- Year %s ---\n", yr))
  yr_int     <- as.integer(yr)
  gee_path   <- gee_paths[[yr]]
  
  if (!file.exists(gee_path)) {
    cat(sprintf("GEE image file not found for %s -- skipped.\n\n", yr)); next
  }
  
  # ---- GT rows for this year ----
  gt_yr <- gt_valid %>% filter(year == yr_int)
  cat(sprintf("GT rows for %s: %d\n", yr, nrow(gt_yr)))
  
  # ---- GEE image data ----
  gee_df   <- read_csv(gee_path, show_col_types = FALSE)
  gse_cols <- grep("^GSE_", names(gee_df), value = TRUE)
  env_cols <- intersect(c("distToLand","rugosity","slope","depth",
                          "elevation","mean_wave_period","sig_wave_height"),
                        names(gee_df))
  
  gee_slim <- gee_df %>%
    select(any_of(c("gee_id", env_cols, gse_cols))) %>%
    # drop xcoord/ycoord from GEE file -- use GT's authoritative coords
    select(-any_of(c("xcoord","ycoord","lat","lon")))
  
  # ---- INNER JOIN: GT x GEE image ----
  merged <- inner_join(gt_yr, gee_slim, by = "gee_id") %>%
    distinct(gee_id, .keep_all = TRUE) %>%
    rename(year_gt = year)
  
  cat(sprintf("After inner join GT x GEE: %d rows\n", nrow(merged)))
  
  # ---- LEFT JOIN: species composition (on compositio) ----
  if ("compositio" %in% names(merged) && "compositio" %in% names(species_df)) {
    sp_slim  <- species_df %>%
      select(compositio, any_of(c("sg_morpho", species_cols)))
    
    # Resolve column conflicts -- keep GT version for shared cols
    overlap  <- setdiff(intersect(names(merged), names(sp_slim)), "compositio")
    sp_join  <- sp_slim %>% select(-any_of(overlap))
    
    merged   <- left_join(merged, sp_join, by = "compositio") %>%
      select(-matches("\\.y$")) %>%
      rename_with(~ sub("\\.x$", "", .), matches("\\.x$"))
    
    n_sp <- sum(!is.na(merged[[species_cols[1]]]))
    cat(sprintf("Species data joined: %d/%d rows have species cols\n", n_sp, nrow(merged)))
  } else {
    cat("compositio column missing -- species join skipped, cols set to NA.\n")
    for (sc in species_cols) merged[[sc]] <- NA_real_
  }
  
  # ---- Fill NA in species cols (0 where PA=1 but no record -> keep as NA) ----
  # Only replace NA with 0 for rows that HAVE species data (at least one col non-NA)
  has_species <- rowSums(!is.na(merged[, intersect(species_cols, names(merged))])) > 0
  for (sc in intersect(species_cols, names(merged))) {
    merged[[sc]] <- ifelse(has_species & is.na(merged[[sc]]), 0, merged[[sc]])
  }
  
  # ---- Derive morph3 ----
  if ("sg_morpho" %in% names(merged)) {
    merged <- merged %>%
      mutate(
        sg_morpho = as.character(sg_morpho),
        morph3 = case_when(
          sg_morpho == "mixed_long"               ~ "mixed_long",
          sg_morpho == "mono"  & Ea_SPC > 0       ~ "mono_Ea",
          sg_morpho %in% c("mixed_short", "mono") ~ "mixed_short_plus_mono_short",
          TRUE ~ NA_character_
        ),
        morph3 = factor(morph3, levels = morph_levels)
      )
  } else {
    merged$morph3 <- factor(NA_character_, levels = morph_levels)
  }
  
  # ---- Derive carbon_index ----
  sp_available <- intersect(species_cols, names(merged))
  if (length(sp_available) > 0 && any(has_species)) {
    n_total_sp <- length(sp_available)
    merged <- merged %>%
      mutate(
        carbon_base = rowSums(
          sweep(across(all_of(sp_available)),
                2, carbon_weights[sp_available], `*`) / 100,
          na.rm = TRUE
        ),
        n_sp_present  = rowSums(across(all_of(sp_available)) > 0, na.rm = TRUE),
        species_factor = ifelse(
          n_total_sp > 1,
          1 + max_bonus * ((n_sp_present - 1) / (n_total_sp - 1)),
          1
        ),
        carbon_index = pmin(1, carbon_base * species_factor)
      ) %>%
      # Set carbon_index to NA for rows with no species data
      mutate(
        carbon_index = ifelse(has_species, carbon_index, NA_real_)
      ) %>%
      select(-carbon_base, -n_sp_present, -species_factor)
  } else {
    merged$carbon_index <- NA_real_
  }
  
  cat(sprintf("carbon_index: %d non-NA | morph3: %d non-NA\n",
              sum(!is.na(merged$carbon_index)),
              sum(!is.na(merged$morph3))))
  
  # ---- Save ----
  out_path <- file.path(result_dir, paste0("master_raw_", yr, ".csv"))
  # Select only columns needed by models (drop GT metadata not used anywhere)
  merged <- merged %>% dplyr::select(dplyr::any_of(master_keep_cols))
  
  write_csv(merged, out_path)
  cat(sprintf("Saved: master_raw_%s.csv (%d rows x %d cols)\n\n",
              yr, nrow(merged), ncol(merged)))
  
  master_list[[yr]] <- merged
}

# =============================================================================
# STEP 4: Per-year and overall summary
# =============================================================================
cat("=== Master Raw Data Summary ===\n")

summary_tbl <- purrr::map_dfr(names(master_list), function(yr) {
  df <- master_list[[yr]]
  tibble::tibble(
    year         = yr,
    n_total      = nrow(df),
    n_PA1        = sum(df$PA == 1, na.rm = TRUE),
    n_PA0        = sum(df$PA == 0, na.rm = TRUE),
    n_tSPC       = sum(!is.na(df$tSPC)),
    n_AGB        = sum(!is.na(df$AGB_pred)),
    n_AGC        = sum(!is.na(df$AGC_pred)),
    n_species    = sum(!is.na(df[[species_cols[1]]])),
    n_carbon_idx = sum(!is.na(df$carbon_index)),
    n_morph3     = sum(!is.na(df$morph3))
  )
})

print(summary_tbl)
write_csv(summary_tbl, file.path(result_dir, "master_raw_summary.csv"))
cat("Saved: master_raw_summary.csv\n\n")

# =============================================================================
# STEP 5: NA diagnostics across all predictors (for PA model)
# =============================================================================
gse_bands     <- paste0("GSE_A", sprintf("%02d", 0:63))
env_diag_vars <- c("depth","distToLand","rugosity","slope","elevation",
                   "mean_wave_period","sig_wave_height")
all_pred      <- c(gse_bands, env_diag_vars)

cat("NA diagnostic per predictor (all years combined):\n")
all_years_df <- bind_rows(master_list)
na_diag <- all_years_df %>%
  summarise(across(all_of(intersect(all_pred, names(all_years_df))),
                   ~ sum(is.na(.)), .names = "na_{.col}")) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_NA") %>%
  mutate(variable = sub("^na_", "", variable)) %>%
  filter(n_NA > 0) %>%
  arrange(desc(n_NA))

if (nrow(na_diag) == 0) {
  cat("No NA values in GEE predictors.\n")
} else {
  cat(sprintf("%d predictors have NA values:\n", nrow(na_diag)))
  print(na_diag, n = Inf)
}

# =============================================================================
# STEP 6: Spatial coverage map (all years, all PA rows)
# =============================================================================
cat("\nGenerating spatial coverage map...\n")

map_df <- all_years_df %>%
  mutate(
    x        = as.numeric(xcoord),
    y        = as.numeric(ycoord),
    PA_label = factor(PA, levels = c(0, 1), labels = c("Absent", "Present"))
  ) %>%
  filter(!is.na(x), !is.na(y))

indonesia_map <- ggplot2::map_data("world", region = "Indonesia")

# Resolve facet column: master_raw saves as 'year_gt', but all_years_df may have 'year'
facet_col <- if ("year" %in% names(map_df)) "year" else "year_gt"

print(
  ggplot() +
    geom_polygon(data = indonesia_map,
                 aes(x = long, y = lat, group = group),
                 fill = "#d4e6c3", colour = "grey50", linewidth = 0.3) +
    geom_point(data = map_df,
               aes(x = x, y = y, colour = PA_label),
               size = 0.8, alpha = 0.4) +
    scale_colour_manual(values = c("Absent" = "#e74c3c", "Present" = "#2980b9"),
                        name = "Seagrass PA") +
    facet_wrap(facets = vars(.data[[facet_col]]), ncol = 4) +
    coord_quickmap(xlim = range(map_df$x) + c(-1, 1),
                   ylim = range(map_df$y) + c(-1, 1)) +
    labs(title = sprintf("Master Raw Data -- All Years (n = %d)", nrow(map_df)),
         x = "Longitude", y = "Latitude") +
    theme_bw(base_size = 10) +
    theme(plot.title      = element_text(face = "bold"),
          legend.position = "bottom",
          strip.text      = element_text(face = "bold"))
)

cat("\n01_build_master_raw.R COMPLETE\n")