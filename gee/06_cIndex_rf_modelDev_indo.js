/************************************************************
 * 06_cIndex_rf_modelDev_indo.js
 * Carbon index / cIndex -- full Indonesia.
 *
 * Adopted 1:1 from the proven chapter_2 (Study 1-era) script
 * (06_cIndex_rf_modelDev_04022026).
 * Changes:
 *   [1] AOI: ecoregion_indo_diss (full Indonesia)
 *   [2] Training: indo_training_complete_YYYY
 *         carbon_ind stored as String -- cleaned via cleanCol()
 *         (Note: the asset property is named "carbon_ind", not
 *         "carbon_index" -- likely truncated on upload, same pattern
 *         as the "P_mixed_lo"/"P_mixed_sh" truncated names noted in
 *         04_COVER_rf_modelDev_indo.js.)
 *   [3] GSE bands: GSE_A00..GSE_A63
 *   [4] Output: chapter_3/output/cIndex/
 *   [5] RF params from R tuning (gee_best_hyperparameters.csv):
 *         cIndex | ntree=500 | mtry=8 | nodesize=9 | bagFraction=0.6
 *         (top-of-file note above once read ntree=300/mtry=4 -- the
 *         values actually used are in the inline comment and RF_PARAMS below)
 *   [6] No .reproject(), depth mask: unmask(0) pattern
 *   [7] No persistence mask in display -- masking in Apps
 *   [8] Predictors: GSE + depth only
 *         distToLand REMOVED -- incomplete coverage causes sample loss
 *
 * K-FOLD CV: unlike 04_COVER_rf_modelDev_indo.js (where K-fold was
 * removed due to a GEE memory limit), this script keeps its own 5-fold
 * CV and exports cv_mae/cv_rmse/cv_r2. This looks like an internal/
 * diagnostic check rather than the reported figure, though -- the
 * thesis (Table S22) reports the carbon index model's Monte Carlo
 * evaluation from C1_rf_model_cIndex.R (100 R-side resampling
 * iterations) as the primary result, consistent with how the R-side
 * and GEE-side models throughout this repository are trained
 * independently and evaluated separately.
 *
 * DEPTH SIGN CONVENTION: matches 02_PROB_rf_modelDev_indo.js,
 * 03_MORPHO_rf_modelDev_indo.js, and 04_COVER_rf_modelDev_indo.js (raw
 * asset treated as POSITIVE, negated for the predictor band) -- now
 * four independent deployment scripts agreeing with each other against
 * 01_trainingData_prep_indo.js's opposite assumption. See the note in
 * 02_PROB_rf_modelDev_indo.js.
 ************************************************************/

// ==========================================================
// CONFIGURATION
// ==========================================================
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var RESPONSE    = 'carbon_ind';
var K_FOLDS     = 5;
var SCALE       = 10;
var SEED        = 42;

var TRAIN_YEARS = [2017, 2018, 2019, 2020, 2021, 2022, 2023];
var APPLY_YEARS = ee.List.sequence(2017, 2024);

var TRAIN_PREFIX  = ASSET_ROOT + '/chapter_3/TrainingPoint/indo_trainingAGC_complete_';
var EXPORT_PREFIX = ASSET_ROOT + '/chapter_3/output/06_CIndex/indo_RF_cIndex2_';

// [5] RF params from R tuning (gee_best_hyperparameters.csv)
// cIndex | Regressor | ntree=500 | mtry=8 | nodesize=9 | bagFraction=0.6
// Best CV RMSE: 0.1055
var RF_PARAMS = {
  numberOfTrees:     500,
  variablesPerSplit:   8,
  minLeafPopulation:   9,
  bagFraction:        0.6,
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
// STATIC LAYERS -- [6] no .reproject(), unmask(0) pattern
// distToLand REMOVED -- incomplete coverage causes sample loss
// ==========================================================
var depthRaw       = ee.Image(ASSET_ROOT + '/chapter_2/InputArea/raster_INA_ACA_Bathymetry');
var depthRawFilled = depthRaw.unmask(0);
var depthMask      = depthRawFilled.gt(0).and(depthRawFilled.lte(500));
var depth          = depthRaw.unmask(1).multiply(-1).rename('depth');


// ==========================================================
// PREDICTOR BANDS -- [3] GSE_A00..GSE_A63 + depth
// ==========================================================
var GSE_BANDS = ee.List.sequence(0, 63).map(function(i) {
  return ee.String('GSE_A').cat(ee.Number(i).format('%02d'));
});

var STATIC_BANDS = ee.List(['depth']);

var FEATURE_BANDS = GSE_BANDS.cat(STATIC_BANDS);
print('cIndex predictor bands:', FEATURE_BANDS);


// ==========================================================
// TRAINING DATA PREP
// [2] carbon_ind stored as String -- 'NA' replaced with -9999 sentinel
// ==========================================================
function cleanCol(f) {
  var val = f.get(RESPONSE);
  var safeVal = ee.Algorithms.If(
    ee.Algorithms.IsEqual(val, 'NA'), -9999, val
  );
  return f.set(RESPONSE, ee.Number.parse(ee.String(safeVal)));
}

function enrichTraining(year) {
  var pts = ee.FeatureCollection(TRAIN_PREFIX + year);
  return pts.map(cleanCol).filter(ee.Filter.neq(RESPONSE, -9999));
}

var training = ee.FeatureCollection(
  TRAIN_YEARS.map(enrichTraining)
).flatten()
 .filter(ee.Filter.gt(RESPONSE, 0))
 .filter(ee.Filter.lte(RESPONSE, 1));

print('Training size:', training.size());
print('Per-year distribution:', training.aggregate_histogram('year'));

Map.addLayer(training, {color: 'cb181d'}, 'cIndex training points', true);


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

print('cIndex model trained');


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
  title: 'Feature Importance -- cIndex',
  hAxis: {title: 'Feature'}, vAxis: {title: 'Importance'},
  legend: {position: 'none'}, colors: ['#cb181d']
}));


// ==========================================================
// MODEL DIAGNOSTICS EXPORT -> chapter_3/models
// ==========================================================
Export.table.toAsset({
  collection:  impFC.limit(66),
  description: 'export_model_importance_cIndex',
  assetId:     ASSET_ROOT + '/chapter_3/models/importance_cIndex'
});
Export.table.toAsset({
  collection:  ee.FeatureCollection([ee.Feature(null, {
    model: 'cIndex', cv_mae: mae, cv_rmse: rmse, cv_r2: r2, n_train: training.size()
  })]),
  description: 'export_model_cv_metrics_cIndex',
  assetId:     ASSET_ROOT + '/chapter_3/models/cv_metrics_cIndex'
});
print('Model diagnostics export queued.');


// ==========================================================
// STACK BUILDER -- [6] no .reproject(), correct depth mask
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
    .updateMask(depthMask)
    .clip(AOI);
}


// ==========================================================
// APPLY PER YEAR & EXPORT -- [7] no persistence mask
// ==========================================================
APPLY_YEARS.getInfo().forEach(function(y) {
  var stack = buildStack(y).select(FEATURE_BANDS);
  var pred  = stack.classify(rf).rename('cIndex_pred');

  Map.addLayer(pred,
    {min: 0, max: 0.3, palette: ['ffffff', 'fcae91', 'fb6a4a', 'cb181d']},
    'cIndex ' + y, false
  );

  Export.image.toAsset({
    image:       pred.clip(AOI),
    description: 'RF_indo_cIndex_' + y,
    assetId:     EXPORT_PREFIX + y,
    region:      AOI,
    scale:       SCALE,
    crs:         'EPSG:4326',
    maxPixels:   1e13
  });

  print('Export queued cIndex:', y);
});
