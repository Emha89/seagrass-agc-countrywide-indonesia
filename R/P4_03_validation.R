# =============================================================================
# P4_03_validation.R
# Paper 4 revised design: leak-free validation in five designs, ordered by
# distance between test and training data.
#
#   pixel     : 10-m pixel groups (identical GSE values in the same year)
#   b100      : 100-m blocks (one Seagrass-Watch site)
#   b1000     : 1-km blocks (between sites)
#   cluster   : DBSCAN clusters, eps 5 km (between meadows; the N_eff unit)
#   ecoregion : MEOW ecoregions, leave-one-ecoregion-out (region without field data)
#
# In every fold, three predictions are made for the same test samples:
#   Tier 1      mean AGC of the training samples
#   Tier 2      mean AGC of the training samples per cover class; the cover model
#               is refitted in the fold and samples enter the class predicted at
#               their own pixel (never their observed cover)
#   Pixel-based GSE-only AGC model
# No information from a test group reaches either model.
#
# Reads : P4_samples.csv, P4_best_params.csv
# Writes: P4_oof_agc_<design>.csv, P4_oof_cov_<design>.csv, P4_fold_balance_<design>.csv
# Runtime: about 1-2 h on a laptop with the default settings (see the log).
# =============================================================================

# setwd("path/to/seagrass-agc-countrywide-indonesia/R")   # folder that holds these scripts
source("P4_00_config_functions.R")

samples <- read_samples()
p_agc <- read_best_params("AGC")
p_cov <- read_best_params("COVER")
cat(sprintf("AGC params: mtry=%d node=%d frac=%.2f | COVER params: mtry=%d node=%d frac=%.2f\n",
            p_agc$mtry, p_agc$min_node_size, p_agc$sample_fraction,
            p_cov$mtry, p_cov$min_node_size, p_cov$sample_fraction))

validate_design <- function(design) {
  gcol  <- DESIGN_COL[[design]]
  f_agc <- out_file(sprintf("P4_oof_agc_%s.csv", design))
  f_cov <- out_file(sprintf("P4_oof_cov_%s.csv", design))
  f_bal <- out_file(sprintf("P4_fold_balance_%s.csv", design))
  if (!CFG$overwrite && file.exists(f_agc) && file.exists(f_cov)) {
    cat(sprintf("\n[%s] output exists, skipped (set CFG$overwrite = TRUE to rerun)\n", design))
    return(invisible(NULL))
  }

  groups <- samples |>
    dplyr::group_by(group = .data[[gcol]]) |>
    dplyr::summarise(n_agc = sum(is_agc),
                     mean_agc = if (any(is_agc)) mean(AGC_pred[is_agc]) else NA_real_,
                     n_cov = sum(is_cov), .groups = "drop") |>
    as.data.frame()

  if (design == "ecoregion") {
    reps <- seq_len(CFG$n_seeds_eco)
    eco_fold <- setNames(as.integer(factor(groups$group)), groups$group)
  } else {
    reps <- seq_len(CFG$n_repeats)
  }
  cat(sprintf("\n[%s] %d groups (%d with AGC), %d repeats\n",
              design, nrow(groups), sum(groups$n_agc > 0), length(reps)))

  out_agc <- list(); out_cov <- list(); bal <- list()
  t0 <- Sys.time()
  for (r in reps) {
    seed_r <- CFG$seed + 1000 * r
    folds  <- if (design == "ecoregion") eco_fold else assign_folds(groups, CFG$k_folds, seed_r)
    s_fold <- unname(folds[samples[[gcol]]])

    for (fi in sort(unique(s_fold))) {
      te     <- s_fold == fi
      tr_agc <- !te & samples$is_agc
      te_agc <-  te & samples$is_agc
      tr_cov <- !te & samples$is_cov
      te_cov <-  te & samples$is_cov
      seed_f <- seed_r + fi

      m_cov <- fit_rf(samples[tr_cov, ], CFG$cover_response, p_cov, seed_f)

      if (any(te_agc)) {
        m_agc  <- fit_rf(samples[tr_agc, ], CFG$agc_response, p_agc, seed_f)
        y_tr   <- samples$AGC_pred[tr_agc]
        cls_tr <- cover_class(predict_rf(m_cov, samples[tr_agc, ]))
        cov_te <- predict_rf(m_cov, samples[te_agc, ])
        cls_te <- cover_class(cov_te)
        class_mean <- tapply(y_tr, cls_tr, mean)
        t1 <- mean(y_tr)
        t2 <- as.numeric(class_mean[as.character(cls_te)])
        fb <- is.na(t2); t2[fb] <- t1

        out_agc[[length(out_agc) + 1]] <- tibble::tibble(
          design = design, rep = r, fold = fi,
          gee_id = samples$gee_id[te_agc], year = samples$year[te_agc],
          ecoregion = samples$ecoregion[te_agc], g_cluster = samples$g_cluster[te_agc],
          obs = samples$AGC_pred[te_agc], sd_ref = samples$sd_ref[te_agc],
          pred_pixel = predict_rf(m_agc, samples[te_agc, ]),
          pred_t1 = t1, pred_t2 = t2, t2_fallback = fb,
          cov_pred = cov_te, cov_class = as.character(cls_te),
          dist_m = nn_distance_m(samples$xcoord[te_agc], samples$ycoord[te_agc],
                                 samples$xcoord[tr_agc], samples$ycoord[tr_agc])
        )
      }
      if (any(te_cov)) {
        out_cov[[length(out_cov) + 1]] <- tibble::tibble(
          design = design, rep = r, fold = fi,
          gee_id = samples$gee_id[te_cov], year = samples$year[te_cov],
          ecoregion = samples$ecoregion[te_cov],
          obs = samples[[CFG$cover_response]][te_cov],
          pred = predict_rf(m_cov, samples[te_cov, ])
        )
      }
      bal[[length(bal) + 1]] <- tibble::tibble(
        design = design, rep = r, fold = fi,
        n_test_agc = sum(te_agc), n_train_agc = sum(tr_agc),
        mean_test_agc = if (any(te_agc)) mean(samples$AGC_pred[te_agc]) else NA_real_,
        n_test_cov = sum(te_cov), n_train_cov = sum(tr_cov)
      )
    }
    cat(sprintf("  [%s] repeat %d/%d done (%.1f min)\n", design, r, length(reps),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }

  readr::write_csv(dplyr::bind_rows(out_agc), f_agc)
  readr::write_csv(dplyr::bind_rows(out_cov), f_cov)
  readr::write_csv(dplyr::bind_rows(bal), f_bal)
  cat(sprintf("  [%s] saved %s, %s, %s\n", design, basename(f_agc), basename(f_cov), basename(f_bal)))
}

for (d in CFG$designs) validate_design(d)

cat("\nP4_03_validation.R COMPLETE\n")
