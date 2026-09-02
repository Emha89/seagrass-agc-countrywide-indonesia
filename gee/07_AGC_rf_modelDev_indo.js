/************************************************************
 * 07_AGC_rf_modelDev_indo.js
 * Above-ground carbon / AGC -- full Indonesia. Final stage of the
 * GEE pipeline. Confirmed final version (the third duplicate of this
 * script; earlier duplicates existed but were superseded).
 *
 * Adopted 1:1 from the proven chapter_2 (Study 1-era) script
 * (07_AGC_rf_modelDev_04022026).
 * Changes:
 *   [1] AOI: ecoregion_indo_diss (full Indonesia)
 *   [2] Training: indo_trainingAGC_complete_YYYY
 *         AGC_simula stored as String -- cleaned via cleanAGC()
 *         (Note: property name truncated from "AGC_simulated" to
 *         "AGC_simula" on upload, same pattern as "carbon_ind" noted
 *         in 06_cIndex_rf_modelDev_indo.js.)
 *   [3] GSE bands: GSE_A00..GSE_A63
 *   [4] Output: chapter_3/output/AGC/
 *   [5] RF params from R tuning (gee_best_hyperparameters.csv):
 *         AGC | ntree=500 | mtry=8 | nodesize=7 | bagFraction=0.8
 *   [6] No .reproject(), depth mask: unmask(0) pattern
 *   [7] No persistence mask in display -- masking in Apps
 *   [8] CONFIRMED PREDICTOR SET (matches F1_rf_model_AGC.R exactly):
 *         GSE bands + tSPC_pred + AGB_pred + cIndex_pred only.
 *         depth, PA_prob, and the morphology probability bands are
 *         present in the code but commented out -- this GEE model uses
 *         the exact same leaner predictor set as the R-side final AGC
 *         model, not the fuller set used by earlier proxy stages.
 *         distToLand REMOVED -- incomplete coverage causes sample loss
 *   [9] All proxy asset prefixes updated to chapter_3 output paths
 *  [10] enrichTraining: all proxy bands extracted from images (not from training columns)
 *  [CI] CI-width filter applied in R (prepare_agc_training withFilter, p75)
 *         -- GEE filter not needed, training asset already reflects that
 *         same filtered subset used to train F1_rf_model_AGC.R
 *  [11] K-fold CV REMOVED -- GEE memory limit exceeded with full Indonesia
 *       dataset (same reason as 04_COVER_rf_modelDev_indo.js). CV metrics
 *       are evaluated via R instead -- specifically F1_rf_model_AGC.R's
 *       Monte Carlo, which for AGC (like AGB) uses 10 bootstrap-expanded
 *       iterations, not 100 (a "100 iterations" note elsewhere in this
 *       script looks like it was carried over from the SPC/cIndex
 *       pattern and doesn't match AGC's actual iteration count).
 *  [12] Map.addLayer training: limit(2000), shown=false -- prevents Too Many Requests
 *  [13] TWO-STEP WORKFLOW to decouple sampleRegions() from export computation graph:
 *         STEP A -- run once: enrichTraining() -> export table asset
 *         STEP B -- run after STEP A completes: load table asset -> train RF -> export images
 *         Toggle RUN_STEP_A below to switch between steps
 *
 * DEPTH SIGN CONVENTION: the depth band is computed here (for the depth
 * mask) using the same unmask(1).multiply(-1) pattern as every other
 * apply-stage script in this repository, even though depth itself isn't
 * an active AGC predictor. See the note in 02_PROB_rf_modelDev_indo.js.
 ************************************************************/

// ==========================================================
// CONFIGURATION
// ==========================================================
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var RESPONSE    = 'AGC_simula';
var SCALE       = 10;
var SEED        = 42;

var TRAIN_YEARS = [2017, 2018, 2019, 2020, 2021, 2022, 2023];
var APPLY_YEARS = ee.List.sequence(2017, 2024);

// [9] All paths updated to chapter_3
var TRAIN_PREFIX  = ASSET_ROOT + '/chapter_3/TrainingPoint/indo_trainingAGC_complete_';
var PA_PREFIX     = ASSET_ROOT + '/chapter_3/output/02_PA/indo_RF_probability2_';
var MORPH_PREFIX  = ASSET_ROOT + '/chapter_3/output/03_MORPH3/indo_MORPH3_probs2_';
var SPC_PREFIX    = ASSET_ROOT + '/chapter_3/output/04_tSPC/indo_RF_tSPC2_';
var AGB_PREFIX    = ASSET_ROOT + '/chapter_3/output/05_AGB/indo_RF_AGB2_';
var CINDEX_PREFIX = ASSET_ROOT + '/chapter_3/output/06_CIndex/indo_RF_cIndex2_';
var EXPORT_PREFIX = ASSET_ROOT + '/chapter_3/output/07_AGC/indo_RF_AGCsimulated_reduceVar_limitCI2_';

// [13] TWO-STEP TOGGLE
// true  -> STEP A: export enriched training table (run once, wait for task to complete)
// false -> STEP B: load table asset, train RF, export images
var RUN_STEP_A = false;

// Asset path for enriched training table (output of STEP A, input of STEP B)
var ENRICHED_TRAINING_ASSET = ASSET_ROOT + '/chapter_3/models/training_AGC_enriched';

// [5] RF params from R tuning (gee_best_hyperparameters.csv)
// AGC | Regressor | ntree=500 | mtry=8 | nodesize=7 | bagFraction=0.8
// Best CV RMSE: 14.697
var RF_PARAMS = {
  numberOfTrees:     500,
  variablesPerSplit:   8,
  minLeafPopulation:   7,
  bagFraction:        0.8,
  seed:               SEED
};

// distToLand REMOVED -- incomplete coverage causes sample loss

// Morph probability band names (as in Script 03 image output)
var MORPH_BANDS = [
  'P_mixed_long',
  'P_mixed_short_plus_mono_short',
  'P_mono_Ea'
];

var EMB_COLL = 'GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL';


// ==========================================================
// AOI -- [1] Full Indonesia
// ==========================================================
var AOI_FC = ee.FeatureCollection(
  ASSET_ROOT + '/packard_fieldData/ecoregion_indo_diss'
);
var AOI = AOI_FC.geometry();
Map.addLayer(AOI, {color: 'yellow'}, 'AOI (Full Indonesia)', false);
Map.centerObject(AOI, 5);


// ==========================================================
// STATIC LAYERS -- [6] no .reproject(), unmask(0) pattern
// distToLand REMOVED -- incomplete coverage causes sample loss
// ==========================================================
var depthRaw       = ee.Image(ASSET_ROOT + '/chapter_2/InputArea/raster_INA_ACA_Bathymetry');
var depthRawFilled = depthRaw.unmask(0);
var depthMask      = depthRawFilled.gt(0).and(depthRawFilled.lte(500));
var depth          = depthRaw.unmask(1).multiply(-1).rename('depth');


// ==========================================================
// PREDICTOR BANDS -- [3][8] GSE_A00..GSE_A63 + proxies only
// (depth, PA_prob, and morphology bands intentionally commented out --
// see [8] above)
// ==========================================================
var GSE_BANDS = ee.List.sequence(0, 63).map(function(i) {
  return ee.String('GSE_A').cat(ee.Number(i).format('%02d'));
});

var FEATURE_BANDS = GSE_BANDS
  //.add('depth')
  //.add('PA_prob')
  //.cat(MORPH_BANDS)
  .add('tSPC_pred')
  .add('AGB_pred')
  .add('cIndex_pred');

print('AGC predictor bands:', FEATURE_BANDS);


// ==========================================================
// PROXY IMAGE LOADERS -- [9] all from chapter_3 assets
// ==========================================================
function loadPA(year)     { return ee.Image(PA_PREFIX     + year).rename('PA_prob'); }
function loadMorph(year)  { return ee.Image(MORPH_PREFIX  + year).select(MORPH_BANDS); }
function loadSPC(year)    { return ee.Image(SPC_PREFIX    + year).rename('tSPC_pred'); }
function loadAGB(year)    { return ee.Image(AGB_PREFIX    + year).rename('AGB_pred'); }
function loadCIndex(year) { return ee.Image(CINDEX_PREFIX + year).rename('cIndex_pred'); }


// ==========================================================
// CLEAN RESPONSE LABEL (String 'NA' -> -9999 sentinel)
// AGC_simula may contain String 'NA' -- parse safely
// ==========================================================
function cleanAGC(f) {
  var val = f.get(RESPONSE);
  var safeVal = ee.Algorithms.If(
    ee.Algorithms.IsEqual(val, 'NA'), -9999, val
  );
  return f.set(RESPONSE, ee.Number.parse(ee.String(safeVal)));
}


// ==========================================================
// ENRICH TRAINING -- [10] extract all proxy bands from images
// sampleRegions() runs here -- decoupled from export graph via STEP A
// ==========================================================
function enrichTraining(year) {
  var pts = ee.FeatureCollection(TRAIN_PREFIX + year);

  var proxyStack = loadPA(year)
    .addBands(loadMorph(year))
    .addBands(loadSPC(year))
    .addBands(loadAGB(year))
    .addBands(loadCIndex(year));

  var proxyColsToExclude = ee.List([
    'PA_prob', 'PA',
    'P_mixed_long', 'P_mixed_short_plus_mono_short', 'P_mono_Ea',
    'P_mixed_lo', 'P_mixed_sh',
    'tSPC_pred', 'tSPC', 'Th_SPC',
    'AGB_pred', 'AGB_simulated',
    'cIndex_pred', 'carbon_ind', 'carbon_i'
  ]);

  var propsToKeep = pts.first().propertyNames().removeAll(proxyColsToExclude);

  var sampled = proxyStack.sampleRegions({
    collection:  pts,
    properties:  propsToKeep,
    scale:       SCALE,
    tileScale:   4,
    geometries:  true
  });

  return sampled
    .map(cleanAGC)
    .filter(ee.Filter.neq(RESPONSE, -9999))
    .filter(ee.Filter.gt(RESPONSE, 0))
    .filter(ee.Filter.lte(RESPONSE, 2000));
}


// ==========================================================
// STACK BUILDER -- [6] no .reproject(), correct depth mask
// [9] All proxy layers from chapter_3 output assets
// distToLand REMOVED -- incomplete coverage causes sample loss
// ==========================================================
function buildStack(year) {
  var emb = ee.ImageCollection(EMB_COLL)
    .filterDate(
      ee.Date.fromYMD(year, 1, 1),
      ee.Date.fromYMD(ee.Number(year).add(1), 1, 1)
    )
    .mosaic()
    .rename(ee.List.sequence(0, 63).map(function(i) {
      return ee.String('GSE_A').cat(ee.Number(i).format('%02d'));
    }));

  return emb
    .addBands(depth)
    .addBands(loadPA(year))
    .addBands(loadMorph(year))
    .addBands(loadSPC(year))
    .addBands(loadAGB(year))
    .addBands(loadCIndex(year))
    .updateMask(depthMask)
    .clip(AOI);
}


// ==========================================================
// [13] STEP A -- export enriched training table
// Run once. Wait for task to complete before running STEP B.
// Set RUN_STEP_A = false after task succeeds.
// ==========================================================
if (RUN_STEP_A) {

  print('>>> STEP A: Building enriched training table...');

  var training_raw = ee.FeatureCollection(
    TRAIN_YEARS.map(enrichTraining)
  ).flatten();

  print('Enriched training size (STEP A):', training_raw.size());
  print('Per-year distribution:', training_raw.aggregate_histogram('year'));

  // [12] limit(2000) for display only -- training stays full for export
  Map.addLayer(
    training_raw.limit(2000),
    {color: '08306b'},
    'AGC training points (sample 2000)',
    false
  );

  Export.table.toAsset({
    collection:  training_raw,
    description: 'export_AGC_enriched_training',
    assetId:     ENRICHED_TRAINING_ASSET
  });

  print('STEP A export queued -- wait for task to complete, then set RUN_STEP_A = false.');

// ==========================================================
// [13] STEP B -- load enriched table, train RF, export images
// Run only after STEP A task completes successfully.
// ==========================================================
} else {

  print('>>> STEP B: Loading enriched training table...');

  var training = ee.FeatureCollection(ENRICHED_TRAINING_ASSET);

  print('Training size (STEP B):', training.size());
  print('Per-year distribution:', training.aggregate_histogram('year'));

  // [12] limit(2000) for display only -- training stays full for the model
  Map.addLayer(
    training.limit(2000),
    {color: '08306b'},
    'AGC training points (sample 2000)',
    false
  );

  // ==========================================================
  // K-FOLD CROSS VALIDATION -- REMOVED [11]
  // CV metrics (MAE, RMSE, R2) are computed via R instead -- GEE memory
  // limit exceeded with the full Indonesia dataset.
  //
  // NOTE: authoritative CV is done in R (F1_rf_model_AGC.R):
  //   - Hyperparameters tuned via grid search (ranger)
  //   - Response: AGC_simulated = runif(AGC_low, AGC_up) per iteration
  //   - Training filtered: AGC_CIwidt <= p75 (withFilter)
  //   - Monte Carlo bootstrap (n=10, not 100 -- see [11] above) for CI
  //     on CV metrics
  // ==========================================================

  // ==========================================================
  // FINAL MODEL -- trained on the entire training dataset
  // ==========================================================
  var rf = ee.Classifier.smileRandomForest(RF_PARAMS)
    .setOutputMode('REGRESSION')
    .train({ features: training, classProperty: RESPONSE, inputProperties: FEATURE_BANDS });

  print('Final RF AGC model trained');

  // ==========================================================
  // VARIABLE IMPORTANCE
  // ==========================================================
  var importance = ee.Dictionary(rf.explain().get('importance'));
  var impFC = ee.FeatureCollection(
    importance.keys().map(function(k) {
      k = ee.String(k);
      return ee.Feature(null, { feature: k, importance: ee.Number(importance.get(k)) });
    })
  ).sort('importance', false);

  print('Feature importance (Top 15):', impFC.limit(15));
  print(ui.Chart.feature.byFeature({
    features: impFC.limit(72), xProperty: 'feature', yProperties: ['importance']
  }).setChartType('ColumnChart').setOptions({
    title: 'Feature Importance -- AGC',
    hAxis: {title: 'Feature'}, vAxis: {title: 'Importance'},
    legend: {position: 'none'}, colors: ['#08306b']
  }));

  // ==========================================================
  // MODEL DIAGNOSTICS EXPORT -> chapter_3/models
  // CV metrics are computed in R -- not stored here
  // ==========================================================
  Export.table.toAsset({
    collection:  impFC.limit(72),
    description: 'export_model_importance_AGC',
    assetId:     ASSET_ROOT + '/chapter_3/models/importance_AGC'
  });
  Export.table.toAsset({
    collection:  ee.FeatureCollection([ee.Feature(null, {
      model:   'AGC',
      n_train: training.size()
      // CV metrics (MAE, RMSE, R2) computed via R (see F1_rf_model_AGC.R)
    })]),
    description: 'export_model_cv_metrics_AGC',
    assetId:     ASSET_ROOT + '/chapter_3/models/cv_metrics_AGC'
  });
  print('Model diagnostics export queued.');

  // ==========================================================
  // APPLY PER YEAR & EXPORT -- model applied across 2017-2024
  // [7] no persistence mask
  // ==========================================================
  APPLY_YEARS.getInfo().forEach(function(y) {
    var stack = buildStack(y).select(FEATURE_BANDS);
    var pred  = stack.classify(rf).rename('AGC_pred');

    Map.addLayer(pred,
      {min: 0, max: 100, palette: ['ffffff', 'c6dbef', '6baed6', '08306b']},
      'AGC ' + y, false
    );

    Export.image.toAsset({
      image:       pred.clip(AOI),
      description: 'RF_indo_AGC_' + y,
      assetId:     EXPORT_PREFIX + y,
      region:      AOI,
      scale:       SCALE,
      crs:         'EPSG:4326',
      maxPixels:   1e13
    });

    print('Export queued AGC:', y);
  });

}
