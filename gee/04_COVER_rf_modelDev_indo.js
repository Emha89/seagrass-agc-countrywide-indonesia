/************************************************************
 * 04_COVER_rf_modelDev_indo.js
 * Seagrass percent cover / tSPC -- full Indonesia.
 * This is the confirmed final version (the third file iteration of
 * this script; earlier duplicates existed but were superseded).
 *
 * Adopted 1:1 from the proven chapter_2 (Study 1-era) script
 * (04_COVER_rf_modelDev_04022026).
 * Original changes from chapter_2:
 *   [1]  AOI       : ecoregion_indo_diss (full Indonesia)
 *   [2]  Training  : indo_training_complete_YYYY
 *   [3]  GSE bands : GSE_A00..GSE_A63
 *   [4]  Output    : chapter_3/output/tSPC/
 *   [5]  RF params : ntree=900, mtry=10, nodesize=9, bagFraction=0.8
 *   [6]  No .reproject(), correct depth mask
 *   [7]  No persistence mask in display
 *   [8]  enrichTraining() extracts morph probs from image at point locations
 *   [9]  tSPC stored as String in asset -- parsed safely in enrichTraining
 *   [10] K-fold CV REMOVED -- GEE memory limit exceeded with full Indonesia
 *        dataset. CV metrics (MAE, RMSE, R2) evaluated via R (100 Monte
 *        Carlo iterations) instead. CONFIRMED against the thesis: Section
 *        3.1.3 reports the SPC model's mean RMSE (16.97%, 95% CI
 *        16.04-17.80) and R2 (0.57) from the R-side 100-iteration Monte
 *        Carlo evaluation on 5,478 samples -- this GEE script's role is
 *        map production, not the reported accuracy figures.
 *   [11] Map.addLayer training: limit(2000), shown=false -- prevents Too
 *        Many Requests
 *
 * Performance fixes over the prior iteration (export results are
 * identical, this only changes how they're computed):
 *   [FIX-1] FEATURE_BANDS converted to client-side JS array
 *           -- eliminates per-tile server-side ee.List evaluation overhead
 *   [FIX-2] Training data pre-materialized as a table asset (2-step workflow)
 *           -- breaks heavy sampleRegions() dependency from export computation graph
 *           -- each export tile no longer re-evaluates sampleRegions over morph images
 *   [FIX-3] buildStack() year arithmetic uses native JS (+1) instead of ee.Number.add()
 *           -- year is a JS number in the forEach loop, no need for server-side math
 *   [FIX-4] APPLY_YEARS as client-side JS array with native forEach()
 *           -- avoids unnecessary .getInfo() round-trip to GEE servers
 *
 * TWO-STEP WORKFLOW:
 *   STEP A -- Set BUILD_TRAINING_ASSET = true
 *            Run script once -> submit task 'export_training_tSPC_enriched'
 *            Wait until the table export task completes in GEE Tasks panel
 *   STEP B -- Set BUILD_TRAINING_ASSET = false
 *            Run script -> loads training from asset -> trains model -> submits image exports
 *            Export tasks are now much faster (no sampleRegions in computation graph)
 *
 * DEPTH SIGN CONVENTION: matches 02_PROB_rf_modelDev_indo.js and
 * 03_MORPHO_rf_modelDev_indo.js (raw asset treated as POSITIVE, negated
 * for the predictor band). See the note in 02_PROB_rf_modelDev_indo.js
 * for the open question about this versus 01_trainingData_prep_indo.js's
 * opposite assumption.
 ************************************************************/


// ==========================================================
// CONFIGURATION
// ==========================================================
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var RESPONSE = 'tSPC';
var SCALE    = 10;
var SEED     = 42;

var TRAIN_YEARS = [2017, 2018, 2019, 2020, 2021, 2022, 2023];

// [FIX-4] Client-side JS array -- no .getInfo() round-trip needed
var APPLY_YEARS = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];

var TRAIN_PREFIX  = ASSET_ROOT + '/chapter_3/TrainingPoint/indo_trainingAGC_complete_';
var MORPH_PREFIX  = ASSET_ROOT + '/chapter_3/output/03_MORPH3/indo_MORPH3_probs2_';
var EXPORT_PREFIX = ASSET_ROOT + '/chapter_3/output/04_tSPC/indo_RF_tSPC2_';

// [FIX-2] Asset path for pre-materialized enriched training data
// STEP A exports enrichTraining() output to this asset (run once)
// STEP B loads this asset for model training (no sampleRegions in graph)
var TRAINING_ASSET = ASSET_ROOT + '/chapter_3/models/training_tSPC_enriched';

// Workflow toggle:
//   true  = STEP A: run enrichTraining() and export table asset (run once)
//   false = STEP B: load training asset, train model, export prediction images
var BUILD_TRAINING_ASSET = false; // <-- set true for STEP A, then false for STEP B

// RF hyperparameters from R tuning -- SPC: ntree=900, mtry=10, nodesize=9, bagFraction=0.8
// Best CV accuracy: 13.7191
var RF_PARAMS = {
  numberOfTrees:     900,
  variablesPerSplit: 10,
  minLeafPopulation:  9,
  bagFraction:       0.8,
  seed:              SEED
};

// Morphology probability bands -- full names as in Script 03 image output
var MORPH_BANDS = [
  'P_mixed_long',
  'P_mixed_short_plus_mono_short',
  'P_mono_Ea'
];
// distToLand REMOVED -- incomplete coverage causes training sample loss

var EMB_COLL = 'GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL';

// [FIX-1] Build FEATURE_BANDS as a client-side JS array
// Previously ee.List (server-side), which added evaluation overhead per tile
// Content is identical: 64 GSE bands + depth + 3 morph probability bands = 68 features
var FEATURE_BANDS = (function() {
  var bands = [];
  for (var i = 0; i < 64; i++) {
    bands.push('GSE_A' + (i < 10 ? '0' + i : '' + i));
  }
  bands.push('depth');
  bands = bands.concat(MORPH_BANDS);
  return bands; // plain JS array -- zero server-side overhead
})();

print('Total predictor bands:', FEATURE_BANDS.length); // expected: 68
print('tSPC predictor bands:', FEATURE_BANDS);


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
// STATIC LAYERS -- [6] no .reproject()
// depth raster stores POSITIVE values (2-2441 m)
// depthMask restricts processing to valid coastal zone (0-500 m)
// distToLand REMOVED -- incomplete coverage causes training sample loss
// ==========================================================
var depthRaw       = ee.Image(ASSET_ROOT + '/chapter_2/InputArea/raster_INA_ACA_Bathymetry');
var depthRawFilled = depthRaw.unmask(0);
var depthMask      = depthRawFilled.gt(0).and(depthRawFilled.lte(500));
var depth          = depthRaw.unmask(1).multiply(-1).rename('depth');

// QA: check bathymetry coverage and depth mask
Map.addLayer(
  depthRaw,
  {min: 0, max: 500, palette: ['#08306b','#2171b5','#6baed6','#c6dbef','#ffffff']},
  'Bathymetry (raw, 0-500m)', false
);
Map.addLayer(
  depthMask.selfMask(),
  {palette: ['#41ab5d']},
  'Depth mask (0-500m valid)', false
);


// ==========================================================
// STACK BUILDER
// Assembles per-year predictor stack: 64 GSE bands + depth + 3 morph bands
// [6]     No .reproject() -- native scale preserved
// [FIX-3] Year arithmetic uses native JS (+1); year is a JS number in forEach loop
// ==========================================================
function buildStack(year) {
  // Build GSE band names as client-side array for rename()
  var gseNames = [];
  for (var i = 0; i < 64; i++) {
    gseNames.push('GSE_A' + (i < 10 ? '0' + i : '' + i));
  }

  // Load annual GSE embedding mosaic for the given year
  var emb = ee.ImageCollection(EMB_COLL)
    .filterDate(
      ee.Date.fromYMD(year, 1, 1),
      ee.Date.fromYMD(year + 1, 1, 1)   // [FIX-3] native JS arithmetic, year is a JS number
    )
    .mosaic()
    .rename(gseNames);

  // Stack: GSE + depth + morph probs, masked to valid depth zone, clipped to AOI
  return emb
    .addBands(depth)
    .addBands(ee.Image(MORPH_PREFIX + year).select(MORPH_BANDS))
    .updateMask(depthMask)
    .clip(AOI);
}


// ==========================================================
// STEP A -- BUILD & EXPORT ENRICHED TRAINING ASSET
// Run ONCE with BUILD_TRAINING_ASSET = true
// This isolates the heavy sampleRegions() computation into a one-time table export
// Once the asset exists, STEP B loads it directly -- no sampleRegions in export graph
// ==========================================================
if (BUILD_TRAINING_ASSET) {

  // Enrich one year of training points with morph probability values
  function enrichTraining(year) {
    var pts   = ee.FeatureCollection(TRAIN_PREFIX + year);
    var morph = ee.Image(MORPH_PREFIX + year).select(MORPH_BANDS);

    // Keep all original properties except old morph columns (replaced by sampled values)
    var propsToKeep = pts.first().propertyNames()
      .removeAll(MORPH_BANDS)
      .removeAll(['P_mixed_lo', 'P_mixed_sh']);

    // Extract morph prob values at each training point location
    var sampled = morph.sampleRegions({
      collection:  pts,
      properties:  propsToKeep,
      scale:       SCALE,
      tileScale:   4,
      geometries:  true
    });

    // Filter valid tSPC values and parse String -> Number
    // tSPC is stored as String in the GEE asset -- must parse before use as regression target [9]
    return sampled
      .filter(ee.Filter.notNull(['tSPC']))
      .filter(ee.Filter.neq('tSPC', 'NA'))
      .filter(ee.Filter.neq('tSPC', 'null'))
      .map(function(f) {
        return f.set('tSPC_num', ee.Number.parse(ee.String(f.get('tSPC'))));
      })
      .filter(ee.Filter.notNull(['tSPC_num']))
      .filter(ee.Filter.gt('tSPC_num', 0))
      .filter(ee.Filter.lte('tSPC_num', 100));
  }

  // Merge enriched training points across all training years
  var trainingRaw = ee.FeatureCollection(
    TRAIN_YEARS.map(enrichTraining)
  ).flatten();

  print('Training sample size (STEP A -- pre-export check):', trainingRaw.size());
  print('Per-year distribution:', trainingRaw.aggregate_histogram('year'));

  // Export enriched training data as a reusable table asset
  // After this task completes, set BUILD_TRAINING_ASSET = false and run STEP B
  Export.table.toAsset({
    collection:  trainingRaw,
    description: 'export_training_tSPC_enriched',
    assetId:     TRAINING_ASSET
  });

  print('=== STEP A: Training asset export queued ===');
  print('Wait for task to complete in Tasks panel, then set BUILD_TRAINING_ASSET = false and re-run.');

} else {

// ==========================================================
// STEP B -- LOAD PRE-BUILT TRAINING -> TRAIN MODEL -> APPLY & EXPORT
// Requires TRAINING_ASSET to exist (complete STEP A first)
// Export computation graph is now lightweight:
//   load asset table -> train RF -> classify pixels -> export image
// sampleRegions() is NOT part of this graph
// ==========================================================

  // [FIX-2] Load pre-materialized training data directly from asset
  var training = ee.FeatureCollection(TRAINING_ASSET);

  print('Training loaded from asset:', training.size());
  print('Per-year distribution:', training.aggregate_histogram('year'));

  // Verify asset integrity before proceeding
  print('tSPC_num stats (verify range 0-100):', training.aggregate_stats('tSPC_num'));
  print('Sample properties (verify morph bands present):', training.first().propertyNames());

  // [11] Display sample only (limit 2000) -- prevents Too Many Requests error in Code Editor
  Map.addLayer(training.limit(2000), {color: 'ff6600'}, 'tSPC training points (sample 2000)', false);

  var RESPONSE_NUM = 'tSPC_num';

  // ==========================================================
  // TRAIN FINAL RF MODEL
  // K-fold CV REMOVED [10] -- GEE memory limit exceeded with full
  // Indonesia dataset. CV metrics (MAE, RMSE, R2) are evaluated via R
  // (100 Monte Carlo iterations) instead -- see D1_rf_model_SPC.R.
  // ==========================================================
  var rf = ee.Classifier.smileRandomForest(RF_PARAMS)
    .setOutputMode('REGRESSION')
    .train({
      features:        training,
      classProperty:   RESPONSE_NUM,
      inputProperties: FEATURE_BANDS  // [FIX-1] client-side array -- no overhead
    });

  print('tSPC model trained');

  // ==========================================================
  // VARIABLE IMPORTANCE
  // ==========================================================
  var importance = ee.Dictionary(rf.explain().get('importance'));
  var impFC = ee.FeatureCollection(
    importance.keys().map(function(k) {
      k = ee.String(k);
      return ee.Feature(null, {
        feature:    k,
        importance: ee.Number(importance.get(k))
      });
    })
  ).sort('importance', false);

  print('Feature importance (Top 10):', impFC.limit(10));
  print(
    ui.Chart.feature.byFeature({
      features:    impFC.limit(66),
      xProperty:   'feature',
      yProperties: ['importance']
    })
    .setChartType('ColumnChart')
    .setOptions({
      title:  'Feature Importance -- tSPC',
      hAxis:  {title: 'Feature'},
      vAxis:  {title: 'Importance'},
      legend: {position: 'none'},
      colors: ['#3d85c8']
    })
  );

  // ==========================================================
  // MODEL DIAGNOSTICS EXPORT -> chapter_3/models
  // CV metrics are computed in R -- not stored here [10]
  // ==========================================================
  Export.table.toAsset({
    collection:  impFC.limit(66),
    description: 'export_model_importance_tSPC',
    assetId:     ASSET_ROOT + '/chapter_3/models/importance_tSPC'
  });
  Export.table.toAsset({
    collection: ee.FeatureCollection([ee.Feature(null, {
      model:   'tSPC',
      n_train: training.size()
      // CV metrics (MAE, RMSE, R2) computed via R (100 Monte Carlo iterations) [10]
    })]),
    description: 'export_model_cv_metrics_tSPC',
    assetId:     ASSET_ROOT + '/chapter_3/models/cv_metrics_tSPC'
  });
  print('Model diagnostics export queued.');

  // ==========================================================
  // APPLY MODEL PER YEAR & EXPORT IMAGE ASSETS (2017-2024)
  // [FIX-4] Native JS forEach -- no .getInfo() round-trip to GEE servers
  // [FIX-1] stack.select() uses client-side array -- no server-side ee.List evaluation
  // [7]     No persistence mask applied to output
  // ==========================================================
  APPLY_YEARS.forEach(function(y) {
    // Build predictor stack for this year and select features in model order
    var stack = buildStack(y).select(FEATURE_BANDS); // [FIX-1] client-side select

    // Apply trained RF model to produce per-pixel tSPC prediction (0-100%)
    var pred = stack.classify(rf).rename('tSPC_pred');

    // Display preview -- not shown by default to avoid Too Many Requests in Code Editor
    Map.addLayer(
      pred,
      {min: 0, max: 100, palette: ['ffffff','d9f0a3','78c679','238443']},
      'tSPC ' + y, false
    );

    // Export prediction image to asset
    Export.image.toAsset({
      image:       pred.clip(AOI),
      description: 'RF_indo_tSPC_' + y,
      assetId:     EXPORT_PREFIX + y,
      region:      AOI,
      scale:       SCALE,
      crs:         'EPSG:4326',
      maxPixels:   1e13
    });

    print('Export queued -- tSPC:', y);
  });

} // end STEP B
