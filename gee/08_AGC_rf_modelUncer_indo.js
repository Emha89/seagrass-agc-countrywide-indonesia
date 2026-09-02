/************************************************************
 * 08_AGC_rf_modelUncer_indo.js
 * National total AGC, Indonesia -- DEFF/ICC uncertainty version.
 * Internally labelled "Script 09b" in the original file; renamed here
 * to match this repository's numbering, since its content (national
 * total AGC with the adopted uncertainty method) is the natural
 * successor to 07_AGC_rf_modelDev_indo.js.
 *
 * This is the script that actually produces the paper's headline
 * national total AGC figures with confidence intervals -- everything
 * upstream (the R-side DEFF/ICC derivation, the per-year Monte Carlo
 * RMSE, the 60 gC/m2 exclusion threshold) converges here.
 *
 * CHANGE FROM THE PRIOR (L=30-based) VERSION: the spatial
 * autocorrelation correction (L = 30 m, a fixed inflation factor of
 * sqrt(pi*L^2/SCALE^2)) is replaced with a cluster survey design
 * effect (DEFF/ICC), which better matches the training data's actual
 * structure -- discrete field campaigns rather than a spatially
 * continuous survey (Kish, 1965; Griffith, 2005).
 *
 * N_EFFECTIVE below is a fixed constant (201.3), matching the pooled
 * value derived in recompute_N_chapter5_DEFF.R (33 DBSCAN clusters,
 * ICC = 0.188, DEFF = 21.1, N = 4,250).
 *
 * SDdata (5.5046) matches the value calc_SDdata_AGC.R computes and
 * instructs you to copy into this script.
 *
 * The 60 gC/m2 upper AGC_MAX bound below is the exclusion threshold
 * whose validity is checked in F3_rf_evalModel_AGC.R's BLOCK 5
 * (prediction reliability declines above this value due to sparse
 * high-density training samples).
 *
 * The mask asset name referenced below
 * (indo_Seagrass_MaxProb_0_7_mask) matches exactly what
 * 02b_seagrassExtent_indo.js exports as its recommended
 * PERSISTENCE_MASK, given MAX_PROB_THRESHOLD=0.7 there.
 *
 * Seagrass extent mask is held STATIC across all years (MaxProb > 0.7
 * pooled across 2017-2024) so that year-to-year AGC changes reflect
 * canopy/composition dynamics within a fixed extent, not changes in
 * the extent itself.
 *
 * Equations: Section 2.5.1 manuscript (revised)
 *   Eq 3: C_total  = AGC_sum x Apx / 1e6  [ton C]
 *   Eq 4: SDtotal  = sqrt(RMSE_yr^2 + SDdata^2)
 *   Eq 7: SDC_adj  = SDtotal x Apx_total / sqrt(N_EFFECTIVE) / 1e6
 *   Eq 8: CI       = C_total +/- 1.96 x SDC_adj
 *
 * CI available: 2017-2023 (RMSE_all from F3_rf_evalModel_AGC.R)
 * 2024: total AGC only -- no field validation data
 ************************************************************/

// ==========================================================
// 1. CONFIGURATION -- edit here only
// ==========================================================
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var TARGET_YEAR    = 2023;   // 2017-2023 (with CI) | 2024 (AGC only)
// Prob threshold removed -- using static MaxProb_0_7_mask asset
var AGC_MIN        = 1;      // gC/m2
var AGC_MAX        = 60;     // gC/m2 -- see F3_rf_evalModel_AGC.R BLOCK 5
var REDUCE_SCALE   = 100;    // m

// ==========================================================
// 2. UNCERTAINTY PARAMETERS
// ==========================================================
var SDdata      = 5.5046;   // from calc_SDdata_AGC.R
var N_EFFECTIVE = 201.3;    // from recompute_N_chapter5_DEFF.R (DEFF/ICC),
                             // fixed constant, pooled across all years, 33 clusters
var SCALE       = 10;
// var L = 30;              // ORIGINAL, no longer used -- kept for reference
// var A_corr = Math.PI * L * L;  // ORIGINAL, no longer used

// RMSE_all from per-year MC evaluation (F3_rf_evalModel_AGC.R)
// 2024 excluded -- no field validation data
var YEARLY_STATS = {
  2017: { RMSE: 8.724  },
  2018: { RMSE: 14.693 },
  2019: { RMSE: 11.451 },
  2020: { RMSE: 9.275  },
  2021: { RMSE: 9.363  },
  2022: { RMSE: 12.876 },
  2023: { RMSE: 10.785 }
};
var CI_YEARS = [2017, 2018, 2019, 2020, 2021, 2022, 2023];

// ==========================================================
// 3. ASSET PATHS
// ==========================================================
var AGC_PREFIX  = ASSET_ROOT + '/chapter_3/output/07_AGC/indo_RF_AGCsimulated_reduceVar_limitCI2_';
// PROB_PREFIX removed -- static mask used instead
var ACA_ASSET   = ASSET_ROOT + '/chapter_2/InputArea/acav2_indoarch_stack';
var DEPTH_ASSET = ASSET_ROOT + '/chapter_2/InputArea/raster_INA_ACA_Bathymetry';

// ==========================================================
// 4. INDONESIA AOI (ZEE boundary)
// ==========================================================
var INDO_AOI = ee.FeatureCollection(
  ASSET_ROOT + '/packard_fieldData/ecoregion_indo_diss'
).geometry();

// ==========================================================
// 5. BASE LAYERS
// ==========================================================
var depthMask = ee.Image(DEPTH_ASSET).unmask(0).gt(0)
                  .and(ee.Image(DEPTH_ASSET).unmask(0).lte(500));

// Static seagrass presence mask: MaxProb > 0.7 across 2017-2024
// Seagrass extent is constant -- AGC changes reflect canopy/composition dynamics
var probMask = ee.Image(
  ASSET_ROOT + '/chapter_3/output/02_PA/indo_Seagrass_MaxProb_0_7_mask'
).select('seagrass_presence');
var geomorphicMask = ee.Image(ACA_ASSET).select('geomorphic').neq(16);

// ==========================================================
// 6. AGC IMAGE -- all masks applied
// ==========================================================
var agcRaw = ee.Image(AGC_PREFIX + TARGET_YEAR).select('AGC_pred');

var agcFiltered = agcRaw
  .updateMask(agcRaw.gte(AGC_MIN).and(agcRaw.lte(AGC_MAX)))
  .updateMask(probMask)
  .updateMask(depthMask)
  .updateMask(geomorphicMask)
  .rename('AGC');

// ==========================================================
// 7. COMPUTE -- single reduceRegion
// ==========================================================
print('=== National AGC - Indonesia ===');
print('Year:', TARGET_YEAR, '| Mask: MaxProb_0_7 (static 2017-2024)');
print('AGC filter:', AGC_MIN, '-', AGC_MAX, 'gC/m2');
print('Scale:', REDUCE_SCALE + 'm | Class 16 excluded');
print('Computing...');

var stats = agcFiltered.reduceRegion({
  reducer:    ee.Reducer.sum()
                .combine(ee.Reducer.count(), null, true)
                .combine(ee.Reducer.mean(),  null, true),
  geometry:   INDO_AOI,
  scale:      REDUCE_SCALE,
  maxPixels:  1e13,
  bestEffort: true,
  tileScale:  8
});

// ==========================================================
// 8. EVALUATE AND PRINT
// ==========================================================
stats.evaluate(function(res) {

  if (!res) {
    print('WARNING: Timeout. Try REDUCE_SCALE = 200.');
    return;
  }

  var AGC_sum  = res['AGC_sum'];
  var N_pixels = res['AGC_count'];
  var AGC_mean = res['AGC_mean'];

  if (!AGC_sum || !N_pixels) {
    print('WARNING: No valid pixels. Check masks.');
    print('Raw result:', res);
    return;
  }

  // Eq 3: total AGC
  var totalC_ton = (AGC_sum * REDUCE_SCALE * REDUCE_SCALE) / 1e6;

  // Area
  var N_10m    = N_pixels * (REDUCE_SCALE / SCALE) * (REDUCE_SCALE / SCALE);
  var area_km2 = (N_10m * SCALE * SCALE) / 1e6;
  var area_m2  = N_10m * SCALE * SCALE;

  print('');
  print('------------------------------');
  print('METHOD A -- Pixel-based');
  print('------------------------------');
  print('Year              :', TARGET_YEAR);
  print('Seagrass area     :', area_km2.toFixed(2), 'km2');
  print('N pixels (10m eq) :', N_10m);
  print('Total AGC         :', totalC_ton.toFixed(2), 'ton C');
  print('Mean AGC (pixel)  :', AGC_mean.toFixed(4), 'gC/m2');

  // CI -- only for 2017-2023
  var hasRMSE = CI_YEARS.indexOf(TARGET_YEAR) !== -1;
  if (!hasRMSE) {
    print('');
    print('WARNING: No CI for ' + TARGET_YEAR + ' -- no field validation data.');
  } else {
    var RMSE_yr = YEARLY_STATS[TARGET_YEAR].RMSE;
    var SDtotal = Math.sqrt(RMSE_yr * RMSE_yr + SDdata * SDdata);  // Eq 4, unchanged

    // ---- CHANGED BLOCK: DEFF/ICC replaces L-based inflation ----
    // ORIGINAL (prior version):
    //   var SDC_ton    = Math.sqrt(N_10m) * SDtotal * (SCALE * SCALE) / 1e6;  // Eq 6
    //   var inflFactor = Math.sqrt(A_corr / (SCALE * SCALE));                 // Eq 7
    //   var SDC_adj    = SDC_ton * inflFactor;
    //
    // NEW: SDC_adj = SDtotal * area_m2 / sqrt(N_EFFECTIVE) / 1e6
    // (area_m2 already computed above as N_10m * SCALE * SCALE)
    var SDC_adj = (SDtotal * area_m2 / Math.sqrt(N_EFFECTIVE)) / 1e6;  // Eq 7 (revised)
    // ---- END CHANGED BLOCK ----

    var CI_half  = 1.96 * SDC_adj;                                  // Eq 8, unchanged
    var CI_lower = Math.max(0, totalC_ton - CI_half);
    var CI_upper = totalC_ton + CI_half;

    print('');
    print('RMSE_yr (Eq 4)    :', RMSE_yr, 'gC/m2');
    print('SDdata  (Eq 4)    :', SDdata, 'gC/m2');
    print('SDtotal (Eq 4)    :', SDtotal.toFixed(4), 'gC/m2');
    print('N_EFFECTIVE       :', N_EFFECTIVE, '(DEFF/ICC, pooled)');
    print('SDC_adj (Eq 7)    :', SDC_adj.toFixed(2), 'ton C');
    print('CI half (Eq 8)    :', CI_half.toFixed(2), 'ton C');
    print('CI lower          :', CI_lower.toFixed(2), 'ton C');
    print('CI upper          :', CI_upper.toFixed(2), 'ton C');
    print('error_lower       :', CI_half.toFixed(2), 'ton C');
    print('error_upper       :', CI_half.toFixed(2), 'ton C');
  }

  print('');
  print('-- Method B input (use in Excel) --');
  print('area_km2          :', area_km2.toFixed(2), 'km2');
  print('C_B = AGC_mean_field x area_km2  [ton C]');
});
