# =============================================================================
# P4_03_validation.R
# Paper 4, design v3: leak-free validation in five designs, ordered by
# distance between test and training data.
#
#   pixel     : 10-m pixel groups (identical GSE values in the same year)
#   b100      : 100-m blocks (one monitoring site)
#   b1000     : 1-km blocks (between sites)
#   cluster   : DBSCAN clusters, eps 5 km (between meadows; the N_eff unit)
#   ecoregion : MEOW ecoregions, leave-one-ecoregion-out (region without field data)
#
# In every fold, the same test samples get one prediction per method:
#   Tier 1           mean AGC of the training samples
#   Tier 2           mean AGC of the training samples of the test sample's ecoregion
#                    when that ecoregion has >= CFG$min_n_eco training samples,
#                    otherwise the Tier 1 value (t2_fallback = TRUE). In the ecoregion
#                    design the test ecoregion has no training data, so Tier 2 = Tier 1.
#   Local mean       mean AGC of the training samples of the test sample's meadow
#                    (DBSCAN cluster) when it has >= CFG$min_n_local of them, otherwise
#                    the Tier 2 value (local_fallback = TRUE). It tests whether the map
#                    adds more than the field samples nearby. In the cluster and
#                    ecoregion designs the meadow is held out, so Local mean = Tier 2.
#   Pixel-based GSE  RF on the 64 GSE bands, all training AGC samples (pred_pixel)
#   Sentinel-2 benchmark, only for test samples with a clear Sentinel-2 pixel:
#     pred_s2        RF on the Sentinel-2 percentiles, training samples with has_s2
#     pred_s2med     RF on the Sentinel-2 medians (CFG$s2_median_check), same samples
#     pred_gse_sub   RF on GSE, the same training samples as pred_s2
#                    (CFG$s2_equal_training = TRUE)
#     pred_t1_sub, pred_t2_sub, pred_local_sub  Tier 1, Tier 2 and local mean from
#                    the same training samples, so every method of the Sentinel-2
#                    comparison uses the same, smaller, sample set
#     dist_s2_m      distance to the nearest training sample of that set
# No information from a test group reaches any model.
#
# The folds and seeds are those of P4v2 (cover-only groups are kept when the folds
# are drawn), so pred_pixel, pred_t1 and dist_m reproduce the P4v2 values.
#
# Reads : P4_samples.csv, P4_best_params.csv
# Writes: P4_oof_agc_<design>.csv, P4_fold_balance_<design>.csv
# Runtime: about 2-3.5 h on a laptop with the default settings (four RF fits per fold).
# =============================================================================

source("R/P4_00_config_functions.R")

samples <- read_samples()
p_agc <- read_best_params("AGC")
S2_ON <- isTRUE(CFG$run_s2) && "has_s2" %in% names(samples) && sum(samples$has_s2) > 0
p_s2med <- NULL
if (S2_ON) {
  p_s2   <- read_best_params("AGC_S2")
  if (isTRUE(CFG$s2_median_check)) p_s2med <- read_best_params("AGC_S2MED")
  p_gsub <- if (isTRUE(CFG$s2_equal_training)) read_best_params("AGC_GSE_SUB") else NULL
}
cat(sprintf("AGC params: mtry=%d node=%d frac=%.2f | Sentinel-2 models: %s\n",
            p_agc$mtry, p_agc$min_node_size, p_agc$sample_fraction,
            if (S2_ON) sprintf("on (%d AGC samples with has_s2)", sum(samples$has_s2)) else "off"))

validate_design <- function(design) {
  gcol  <- DESIGN_COL[[design]]
  f_agc <- out_file(sprintf("P4_oof_agc_%s.csv", design))
  f_bal <- out_file(sprintf("P4_fold_balance_%s.csv", design))
  if (!CFG$overwrite && file.exists(f_agc)) {
    cat(sprintf("\n[%s] output exists, skipped (set CFG$overwrite = TRUE to rerun)\n", design))
    return(invisible(NULL))
  }

  # all samples of the sample build, so the folds of the AGC groups are those of P4v2
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

  out_agc <- list(); bal <- list()
  t0 <- Sys.time()
  for (r in reps) {
    seed_r <- CFG$seed + 1000 * r
    folds  <- if (design == "ecoregion") eco_fold else assign_folds(groups, CFG$k_folds, seed_r)
    s_fold <- unname(folds[samples[[gcol]]])

    for (fi in sort(unique(s_fold))) {
      te     <- s_fold == fi
      tr_agc <- !te & samples$is_agc
      te_agc <-  te & samples$is_agc
      seed_f <- seed_r + fi
      if (!any(te_agc)) next

      m_agc <- fit_rf(samples[tr_agc, ], CFG$agc_response, p_agc, seed_f)
      y_tr  <- samples$AGC_pred[tr_agc]
      t1    <- mean(y_tr)
      t2    <- tier2_predict(tier2_fit(y_tr, samples$ecoregion[tr_agc]), samples$ecoregion[te_agc])
      loc   <- local_predict(y_tr, samples$g_cluster[tr_agc], samples$g_cluster[te_agc], t2$pred)

      res <- tibble::tibble(
        design = design, rep = r, fold = fi,
        gee_id = samples$gee_id[te_agc], year = samples$year[te_agc],
        ecoregion = samples$ecoregion[te_agc], g_cluster = samples$g_cluster[te_agc],
        g_b100 = samples$g_b100[te_agc],
        has_s2 = if (S2_ON) samples$has_s2[te_agc] else FALSE,
        obs = samples$AGC_pred[te_agc], sd_ref = samples$sd_ref[te_agc],
        pred_pixel = predict_rf(m_agc, samples[te_agc, ]),
        pred_t1 = t1, pred_t2 = t2$pred, t2_fallback = t2$fallback,
        pred_local = loc$pred, local_fallback = loc$fallback,
        pred_gse_sub = NA_real_, pred_s2 = NA_real_, pred_s2med = NA_real_,
        pred_t1_sub = NA_real_, pred_t2_sub = NA_real_, pred_local_sub = NA_real_,
        dist_s2_m = NA_real_,
        dist_m = nn_distance_m(samples$xcoord[te_agc], samples$ycoord[te_agc],
                               samples$xcoord[tr_agc], samples$ycoord[tr_agc])
      )

      n_tr_s2 <- 0L; n_te_s2 <- 0L
      if (S2_ON) {
        tr_s2 <- tr_agc & samples$has_s2
        te_s2 <- te_agc & samples$has_s2
        n_tr_s2 <- sum(tr_s2); n_te_s2 <- sum(te_s2)
        if (n_te_s2 > 0 && n_tr_s2 >= 20) {
          k <- res$has_s2                         # rows of te_s2, in the same order
          y_s2 <- samples$AGC_pred[tr_s2]
          t2s  <- tier2_predict(tier2_fit(y_s2, samples$ecoregion[tr_s2]), samples$ecoregion[te_s2])
          res$pred_t1_sub[k]    <- mean(y_s2)
          res$pred_t2_sub[k]    <- t2s$pred
          res$pred_local_sub[k] <- local_predict(y_s2, samples$g_cluster[tr_s2],
                                                 samples$g_cluster[te_s2], t2s$pred)$pred
          res$dist_s2_m[k] <- nn_distance_m(samples$xcoord[te_s2], samples$ycoord[te_s2],
                                            samples$xcoord[tr_s2], samples$ycoord[tr_s2])
          m_s2 <- fit_rf(samples[tr_s2, ], CFG$agc_response, p_s2, seed_f, CFG$s2_vars)
          res$pred_s2[k] <- predict_rf(m_s2, samples[te_s2, ])
          if (!is.null(p_s2med)) {
            m_s2med <- fit_rf(samples[tr_s2, ], CFG$agc_response, p_s2med, seed_f, CFG$s2_med_vars)
            res$pred_s2med[k] <- predict_rf(m_s2med, samples[te_s2, ])
          }
          if (!is.null(p_gsub)) {
            m_gsub <- fit_rf(samples[tr_s2, ], CFG$agc_response, p_gsub, seed_f)
            res$pred_gse_sub[k] <- predict_rf(m_gsub, samples[te_s2, ])
          } else {
            res$pred_gse_sub[k] <- res$pred_pixel[k]
          }
        }
      }
      out_agc[[length(out_agc) + 1]] <- res

      bal[[length(bal) + 1]] <- tibble::tibble(
        design = design, rep = r, fold = fi,
        n_test_agc = sum(te_agc), n_train_agc = sum(tr_agc),
        mean_test_agc = mean(samples$AGC_pred[te_agc]),
        n_test_s2 = n_te_s2, n_train_s2 = n_tr_s2,
        pct_tier2_fallback = 100 * mean(t2$fallback),
        pct_local_fallback = 100 * mean(loc$fallback)
      )
    }
    cat(sprintf("  [%s] repeat %d/%d done (%.1f min)\n", design, r, length(reps),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }

  readr::write_csv(dplyr::bind_rows(out_agc), f_agc)
  readr::write_csv(dplyr::bind_rows(bal), f_bal)
  cat(sprintf("  [%s] saved %s, %s\n", design, basename(f_agc), basename(f_bal)))
}

for (d in CFG$designs) validate_design(d)

cat("\nP4_03_validation.R COMPLETE\n")
