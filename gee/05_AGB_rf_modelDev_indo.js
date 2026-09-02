/************************************************************
 * 05_AGB_rf_modelDev_indo.js
 * Above-ground biomass / AGB -- full Indonesia.
 *
 * Adopted 1:1 from the proven chapter_2 (Study 1-era) script
 * (05_AGB_rf_modelDev_04022026).
 * Changes:
 *   [1] AOI: ecoregion_indo_diss (full Indonesia)
 *   [2] Training: indo_training_complete_YYYY
 *         AGB_pred stored as String -- cleaned via cleanAGB()
 *   [3] GSE bands: GSE_A00..GSE_A63
 *   [4] Output: chapter_3/output/AGB/
 *   [5] RF params from R tuning (gee_best_hyperparameters.csv):
 *         AGB | ntree=500 | mtry=10 | nodesize=7 | bagFraction=0.7
 *         (top-of-file note above once read ntree=300/mtry=5 -- the
 *         values actually used are in the inline comment and RF_PARAMS below)
 *   [6] No .reproject(), correct depth mask (raw raster = positive 2-2441)
 *         depthRaw.unmask(0) to preserve coastal areas with no depth data
 *   [7] No persistence mask in display -- masking in Apps
 *   [8] Predictors: GSE + depth only
 *         distToLand REMOVED -- incomplete coverage causes sample loss
 *
 * K-FOLD CV: retained here (unlike 04_COVER_rf_modelDev_indo.js, where
 * it was removed for a GEE memory limit). As with cIndex, this looks
 * like an internal/diagnostic check rather than the reported figure --
 * the thesis (Table S17) reports the AGB model's evaluation from
 * E1_rf_model_AGB.R's R-side Monte Carlo, which for AGB specifically
 * uses 10 bootstrap-expanded iterations (not the 100 used for PA,
 * morphology, SPC, and carbon index).
 *
 * DEPTH SIGN CONVENTION: matches every other apply-stage script in this
 * repository (raw asset treated as POSITIVE, negated for the predictor
 * band) -- now five independent deployment scripts agreeing with each
 * other against 01_trainingData_prep_indo.js's opposite assumption. See
 * the note in 02_PROB_rf_modelDev_indo.js.
 ************************************************************/

// ==========================================================
// CONFIGURATION
// ==========================================================
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var RESPONSE    = 'AGB_pred';
var K_FOLDS     = 5;
var SCALE       = 10;
var SEED        = 42;

var TRAIN_YEARS = [2017, 2018, 2019, 2020, 2021, 2022, 2023];
var APPLY_YEARS = ee.List.sequence(2017, 2024);

var TRAIN_PREFIX  = ASSET_ROOT + '/chapter_3/TrainingPoint/indo_trainingAGC_complete_';
var EXPORT_PREFIX = ASSET_ROOT + '/chapter_3/output/05_AGB/indo_RF_AGB2_';

// [5] RF params from R tuning (gee_best_hyperparameters.csv)
// AGB | Regressor | ntree=500 | mtry=10 | nodesize=7 | bagFraction=0.7
// Best CV RMSE: 73.3127
var RF_PARAMS = {
  numberOfTrees:     500,
  variablesPerSplit:   10,
  minLeafPopulation:   7,
  bagFraction:        0.7,
  seed:               SEED
};
// distToLand REMOVED -- incomplete coverage causes sample loss

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
// STATIC LAYERS -- [6] no .reproject()
// depthRaw: POSITIVE values (2-2441 m); unmask(0) so null areas
// become 0 and are excluded by gt(0) rather than masking everything out
// distToLand REMOVED -- incomplete coverage causes sample loss
// ==========================================================
var depthRaw       = ee.Image(ASSET_ROOT + '/chapter_2/InputArea/raster_INA_ACA_Bathymetry');
var depthRawFilled = depthRaw.unmask(0);
var depthMask      = depthRawFilled.gt(0).and(depthRawFilled.lte(500));
var depth          = depthRaw.unmask(1).multiply(-1).rename('depth');


// ==========================================================
// PREDICTOR BANDS -- [3] [8]
// ==========================================================
var GSE_BANDS = ee.List.sequence(0, 63).map(function(i) {
  return ee.String('GSE_A').cat(ee.Number(i).format('%02d'));
});

var STATIC_BANDS = ee.List(['depth']);

var FEATURE_BANDS = GSE_BANDS.cat(STATIC_BANDS);
print('AGB predictor bands:', FEATURE_BANDS);


// ==========================================================
// TRAINING DATA PREP
// [2] AGB_pred stored as String -- 'NA' replaced with -9999 sentinel
// ==========================================================
function cleanAGB(f) {
  var val = f.get(RESPONSE);
  var safeVal = ee.Algorithms.If(
    ee.Algorithms.IsEqual(val, 'NA'), -9999, val
  );
  return f.set(RESPONSE, ee.Number.parse(ee.String(safeVal)));
}

function enrichTraining(year) {
  var pts = ee.FeatureCollection(TRAIN_PREFIX + year);
  return pts.map(cleanAGB).filter(ee.Filter.neq(RESPONSE, -9999));
}

var training = ee.FeatureCollection(
  TRAIN_YEARS.map(enrichTraining)
).flatten()
 .filter(ee.Filter.gt(RESPONSE, 0));

print('Training size:', training.size());
print('Per-year distribution:', training.aggregate_histogram('year'));

Map.addLayer(training, {color: 'ff6600'}, 'AGB training points', true);


// ==========================================================
// K-FOLD CROSS VALIDATION
// ==========================================================
var trainingWithFold = training.randomColumn('rand', SEED).map(function(f) {
  return f.set('fold', ee.Number(f.get('rand')).multiply(K_FOLDS).int());
});
print('Fold distribution:', trainingWithFold.aggregate_histogram('fold'));

var cv = ee.FeatureCollection(
  ee.List.sequence(0, K_FOLDS - 1).map(function(k) {
    k = ee.Number(k);
    var trainFold = trainingWithFold.filter(ee.Filter.neq('fold', k));
    var validFold = trainingWithFold.filter(ee.Filter.eq('fold', k));

    var model = ee.Classifier.smileRandomForest(RF_PARAMS)
      .setOutputMode('REGRESSION')
      .train({ features: trainFold, classProperty: RESPONSE, inputProperties: FEATURE_BANDS });

    return validFold.classify(model, 'pred').map(function(f) {
      var y    = ee.Number(f.get(RESPONSE));
      var yhat = ee.Number(f.get('pred'));
      var diff = y.subtract(yhat);
      return f.set({ fold: k, abs_res: diff.abs(), sq_res: diff.pow(2) });
    });
  })
).flatten();

var mae   = cv.aggregate_mean('abs_res');
var rmse  = ee.Number(cv.aggregate_mean('sq_res')).sqrt();
var yMean = ee.Number(cv.aggregate_mean(RESPONSE));
var cv2   = cv.map(function(f) {
  return f.set('sq_tot', ee.Number(f.get(RESPONSE)).subtract(yMean).pow(2));
});
var r2 = ee.Number(1).subtract(
  ee.Number(cv2.aggregate_sum('sq_res')).divide(ee.Number(cv2.aggregate_sum('sq_tot')))
);
print('CV MAE:', mae, '| CV RMSE:', rmse, '| CV R2:', r2);


// ==========================================================
// FINAL MODEL
// ==========================================================
var rf = ee.Classifier.smileRandomForest(RF_PARAMS)
  .setOutputMode('REGRESSION')
  .train({ features: training, classProperty: RESPONSE, inputProperties: FEATURE_BANDS });

print('Final RF AGB model trained');


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

print('Feature importance (Top 10):', impFC.limit(10));
print(ui.Chart.feature.byFeature({
  features: impFC.limit(66), xProperty: 'feature', yProperties: ['importance']
}).setChartType('ColumnChart').setOptions({
  title: 'Feature Importance -- AGB',
  hAxis: {title: 'Feature'}, vAxis: {title: 'Importance'},
  legend: {position: 'none'}, colors: ['#e6550d']
}));


// ==========================================================
// MODEL DIAGNOSTICS EXPORT -> chapter_3/models
// ==========================================================
Export.table.toAsset({
  collection:  impFC.limit(66),
  description: 'export_model_importance_AGB',
  assetId:     ASSET_ROOT + '/chapter_3/models/importance_AGB'
});
Export.table.toAsset({
  collection:  ee.FeatureCollection([ee.Feature(null, {
    model: 'AGB', cv_mae: mae, cv_rmse: rmse, cv_r2: r2, n_train: training.size()
  })]),
  description: 'export_model_cv_metrics_AGB',
  assetId:     ASSET_ROOT + '/chapter_3/models/cv_metrics_AGB'
});
print('Model diagnostics export queued.');


// ==========================================================
// STACK BUILDER -- [6] no .reproject(), correct depth mask
// [8] GSE + depth only (distToLand REMOVED)
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
    .updateMask(depthMask)
    .clip(AOI);
}


// ==========================================================
// APPLY PER YEAR & EXPORT -- [7] no persistence mask
// ==========================================================
APPLY_YEARS.getInfo().forEach(function(y) {
  var stack = buildStack(y).select(FEATURE_BANDS);
  var pred  = stack.classify(rf).rename('AGB_pred');

  Map.addLayer(pred,
    {min: 0, max: 300, palette: ['ffffff', 'fee6ce', 'fdae6b', 'e6550d']},
    'AGB ' + y, false
  );

  Export.image.toAsset({
    image:       pred.clip(AOI),
    description: 'RF_indo_AGB_' + y,
    assetId:     EXPORT_PREFIX + y,
    region:      AOI,
    scale:       SCALE,
    crs:         'EPSG:4326',
    maxPixels:   1e13
  });

  print('Export queued AGB:', y);
});
