/************************************************************
 * 02_PROB_rf_modelDev_indo.js
 * PA probability model -- full Indonesia.
 *
 * Adopted 1:1 from the proven chapter_2 (Study 1-era) script.
 * Changes from original:
 *   [1] AOI: R2_case_study -> ecoregion_indo_diss (full Indonesia)
 *   [2] Training: training_embedding_depth3_ -> indo_trainingAGC_complete_
 *       (Script 01's export prefix)
 *   [3] GSE bands: A00..A63 -> GSE_A00..GSE_A63 (match Script 01 output)
 *   [4] Output: chapter_2/output -> chapter_3/output/PA/
 *   [5] RF params: confirmed from R tuning (gee_best_hyperparameters.csv)
 *
 * OPEN QUESTION -- DEPTH SIGN CONVENTION: this script's comments state
 * the raw bathymetry raster holds POSITIVE values (2-2441 m), and the
 * mask below is built directly from that assumption
 * (depthRaw.gt(0).and(depthRaw.lte(500))). The predictor band is then
 * built as depthRaw.multiply(-1) -- i.e. NEGATIVE depth values.
 * 01_trainingData_prep_indo.js uses the exact same source asset but
 * states the opposite: "bathymetry stored as negative, converted to
 * positive depth" via the same multiply(-1) step, producing POSITIVE
 * depth values for the training data. If both comments are accurate
 * about their own script, the "depth" predictor would carry opposite
 * signs between training (Script 01, positive) and deployment (this
 * script, negative) -- worth verifying directly (e.g. print/inspect
 * actual depth values from both) before relying on any output that
 * used this stack. Left unchanged here since this is what actually
 * produced the existing results; flagging rather than guessing at a fix.
 *
 * Mask: depthRaw.gt(0).and(depthRaw.lte(500))
 ************************************************************/

// ==========================================================
// CONFIGURATION
// ==========================================================
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var RESPONSE = 'PA';
var K_FOLDS  = 5;
var SCALE    = 10;
var SEED     = 42;

var START_YEAR = 2017;
var END_YEAR   = 2024;

var EXPORT_BASE        = ASSET_ROOT + '/chapter_3/output';
var TRAIN_ASSET_PREFIX = ASSET_ROOT + '/chapter_3/TrainingPoint/indo_trainingAGC_complete_';

var PROB_THRESHOLD = 0.7;
var APPLY_THRESHOLD = {
  2017: 0.5, 2018: 0.5, 2019: 0.5, 2020: 0.5,
  2021: 0.5, 2022: 0.5, 2023: 0.5, 2024: 0.5
};

var STATIC_LAYER_CONFIG = {
  'depth':    { include: true  },
  'slope':    { include: false },
  'rugosity': { include: false },
  'wave':     { include: false }
};
// distToLand REMOVED -- incomplete coverage causes sample loss (see
// 01_trainingData_prep_indo.js)

// [5] RF hyperparameters from R tuning (gee_best_hyperparameters.csv)
// PA | Classifier | ntree=500 | mtry=10 | nodesize=3 | bagFraction=0.8
// Best CV accuracy: 0.9317
// (The top-of-file summary above once read ntree=700/mtry=9 -- that
// looks like a stale note; the values actually used are the ones below.)
var RF_PARAMS = {
  numberOfTrees:     500,
  variablesPerSplit: 10,
  minLeafPopulation: 3,
  bagFraction:       0.8,
  seed:              SEED
};

var APPLY_YEARS = ee.List.sequence(START_YEAR, END_YEAR);

// ==========================================================
// PREDICTOR BANDS
// [3] Prefix changed from 'A' to 'GSE_A' to match Script 01 output
// ==========================================================
var GSE_BANDS = ee.List.sequence(0, 63).map(function(i) {
  return ee.String('GSE_A').cat(ee.Number(i).format('%02d'));
});

var STATIC_BANDS = ee.List([]);
if (STATIC_LAYER_CONFIG['depth'].include)    STATIC_BANDS = STATIC_BANDS.add('depth');
if (STATIC_LAYER_CONFIG['slope'].include)    STATIC_BANDS = STATIC_BANDS.add('slope');
if (STATIC_LAYER_CONFIG['rugosity'].include) STATIC_BANDS = STATIC_BANDS.add('rugosity');
if (STATIC_LAYER_CONFIG['wave'].include)     STATIC_BANDS = STATIC_BANDS.cat(['wElevation','wHeight','wPeriod']);

var FEATURE_BANDS = GSE_BANDS.cat(STATIC_BANDS);
print('Predictors:', FEATURE_BANDS);

// ==========================================================
// AOI
// [1] Full Indonesia ecoregion boundary
// ==========================================================
var AOI_FC   = ee.FeatureCollection(
  ASSET_ROOT + '/packard_fieldData/ecoregion_indo_diss'
);
var AOI_GEOM = AOI_FC.geometry();

// ==========================================================
// LOAD TRAINING DATA
// ==========================================================
var training_all = ee.FeatureCollection([]);
var available_years = [2018, 2019, 2020, 2021, 2022, 2023];

available_years.forEach(function(y) {
  training_all = training_all.merge(ee.FeatureCollection(TRAIN_ASSET_PREFIX + y));
});
print('Raw training size:', training_all.size());

function notNullFor(list) {
  return ee.Filter.and.apply(null, ee.List(list).map(function(b) {
    return ee.Filter.notNull([b]);
  }));
}
training_all = training_all.filter(notNullFor(FEATURE_BANDS.add(RESPONSE)));

// Cast PA from Boolean to Integer (0/1)
training_all = training_all.map(function(f) {
  var paVal = ee.Algorithms.If(f.get(RESPONSE), 1, 0);
  return f.set(RESPONSE, ee.Number(paVal));
});

print('Filtered training size:', training_all.size());
print('PA distribution:', training_all.aggregate_histogram(RESPONSE));

// ==========================================================
// TRAINING SAMPLE DISTRIBUTION -- QA BEFORE EXPORT
// ==========================================================
var train_present = training_all.filter(ee.Filter.eq(RESPONSE, 1));
var train_absent  = training_all.filter(ee.Filter.eq(RESPONSE, 0));

print('--- Training sample QA ---');
print('PA = 1 (present):', train_present.size());
print('PA = 0 (absent):', train_absent.size());
print('Per-year distribution:', training_all.aggregate_histogram('year'));

Map.addLayer(train_absent,  {color: 'ff4444'}, 'Training -- PA=0 (absent)',  true);
Map.addLayer(train_present, {color: '0066ff'}, 'Training -- PA=1 (present)', true);
Map.addLayer(AOI_GEOM, {color: 'red'}, 'AOI (Full Indonesia)', false);
Map.centerObject(AOI_GEOM, 5);

// ==========================================================
// K-FOLD VALIDATION
// ==========================================================
var withFold = training_all.randomColumn('rand', SEED).map(function(f) {
  return f.set('fold', ee.Number(f.get('rand')).multiply(K_FOLDS).int());
});
print('Fold distribution:', withFold.aggregate_histogram('fold'));

var foldIndices   = ee.List.sequence(0, K_FOLDS - 1);
var cvPredictions = foldIndices.map(function(k) {
  k = ee.Number(k);
  var train = withFold.filter(ee.Filter.neq('fold', k));
  var valid = withFold.filter(ee.Filter.eq('fold', k));

  var rf = ee.Classifier.smileRandomForest(RF_PARAMS).train({
    features: train, classProperty: RESPONSE, inputProperties: FEATURE_BANDS
  });
  return valid.classify(rf.setOutputMode('PROBABILITY'), 'prob').map(function(f) {
    return f.set('pred', ee.Number(f.get('prob')).gte(PROB_THRESHOLD).int());
  });
});

var cvResult = ee.FeatureCollection(cvPredictions).flatten();
var cm       = cvResult.errorMatrix(RESPONSE, 'pred');
print('Confusion matrix:', cm);
print('Accuracy:', cm.accuracy());
print('Kappa:', cm.kappa());

var arr  = ee.Array(cm.array());
var TP   = arr.get([1,1]); var FP = arr.get([0,1]); var FN = arr.get([1,0]);
var precision = TP.divide(TP.add(FP));
var recall    = TP.divide(TP.add(FN));
var f1        = precision.multiply(recall).multiply(2).divide(precision.add(recall));
print('Precision:', precision); print('Recall:', recall); print('F1:', f1);

// ==========================================================
// FINAL MODEL
// ==========================================================
var finalRF = ee.Classifier.smileRandomForest(RF_PARAMS).train({
  features: training_all, classProperty: RESPONSE, inputProperties: FEATURE_BANDS
});
print('Final model trained');

// ==========================================================
// VARIABLE IMPORTANCE
// ==========================================================
var importanceDict = ee.Dictionary(finalRF.explain().get('importance'));
var sortedFC       = ee.FeatureCollection(
  importanceDict.keys().map(function(k) {
    k = ee.String(k);
    return ee.Feature(null, {feature: k, importance: ee.Number(importanceDict.get(k))});
  })
).sort('importance', false);

print('Feature Importance (Top 10):', sortedFC.limit(10));
print(ui.Chart.feature.byFeature({
  features: sortedFC.limit(66), xProperty: 'feature', yProperties: ['importance']
}).setChartType('ColumnChart').setOptions({
  title: 'Feature Importance -- PA', hAxis: {title: 'Feature'},
  vAxis: {title: 'Importance'}, legend: {position: 'none'}, colors: ['#444']
}));

// ==========================================================
// MODEL DIAGNOSTICS EXPORT -> chapter_3/models
// ==========================================================
Export.table.toAsset({
  collection:  sortedFC.limit(66),
  description: 'export_model_importance_PA',
  assetId:     ASSET_ROOT + '/chapter_3/models/importance_PA'
});
Export.table.toAsset({
  collection:  ee.FeatureCollection([ee.Feature(null, {
    model: 'PA', cv_acc: cm.accuracy(), cv_kappa: cm.kappa(),
    cv_f1: f1,   n_train: training_all.size()
  })]),
  description: 'export_model_cv_metrics_PA',
  assetId:     ASSET_ROOT + '/chapter_3/models/cv_metrics_PA'
});
print('Model diagnostics export queued.');

// ==========================================================
// STACK BUILDER -- no .reproject(), identical to original
// [3] GSE mosaic renamed GSE_A00..GSE_A63 after .mosaic()
// See the depth sign-convention note at the top of this file.
// distToLand REMOVED -- incomplete coverage causes sample loss
// ==========================================================
function buildStackForYear(y) {
  y = ee.Number(y);

  var emb = ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL')
    .filterDate(ee.Date.fromYMD(y, 1, 1), ee.Date.fromYMD(y.add(1), 1, 1))
    .mosaic()
    .rename(ee.List.sequence(0, 63).map(function(i) {
      return ee.String('GSE_A').cat(ee.Number(i).format('%02d'));
    }));

  var depthRaw = ee.Image(ASSET_ROOT + '/chapter_2/InputArea/raster_INA_ACA_Bathymetry');
  var mask     = depthRaw.gt(0).and(depthRaw.lte(500));
  var depth    = depthRaw.multiply(-1).rename('depth');

  var slope = ee.Terrain.slope(depth).rename('slope');
  var rugosity = slope.reduceNeighborhood({
    reducer: ee.Reducer.stdDev(), kernel: ee.Kernel.circle(9, 'meters')
  }).unmask(0).rename('rugosity');

  var wave = ee.Image(ASSET_ROOT + '/chapter_2/InputArea/era5_gebco_multiband_2024_450m_indo')
    .select([0, 1, 2], ['wElevation', 'wHeight', 'wPeriod']);

  var stack = emb;
  if (STATIC_LAYER_CONFIG['depth'].include)    stack = stack.addBands(depth);
  if (STATIC_LAYER_CONFIG['slope'].include)    stack = stack.addBands(slope);
  if (STATIC_LAYER_CONFIG['rugosity'].include) stack = stack.addBands(rugosity);
  if (STATIC_LAYER_CONFIG['wave'].include)     stack = stack.addBands(wave);

  return stack.updateMask(mask).clip(AOI_GEOM);
}

// ==========================================================
// APPLY MODEL BY YEAR & EXPORT
// ==========================================================
APPLY_YEARS.evaluate(function(years) {
  years.forEach(function(y) {
    var yearNum = parseInt(y);
    var stack   = buildStackForYear(yearNum).select(FEATURE_BANDS);
    var prob    = stack.classify(finalRF.setOutputMode('PROBABILITY')).rename('prob');
    var thr     = APPLY_THRESHOLD[yearNum] || PROB_THRESHOLD;
    var mask    = prob.gte(thr).rename('seagrass').selfMask();

    Map.addLayer(prob, {min: 0, max: 1, palette: ['#d73027', '#fee08b', '#1a9850']},
      'Prob ' + yearNum, false);
    Map.addLayer(mask, {palette: ['#1a9850']}, 'Mask ' + yearNum, false);

    Export.image.toAsset({
      image: prob.clip(AOI_GEOM), description: 'export_indo_RF_prob_' + yearNum,
      assetId: EXPORT_BASE + '/02_PA/indo_RF_probability2_' + yearNum,
      region: AOI_GEOM, scale: SCALE, crs: 'EPSG:4326', maxPixels: 1e13
    });
    Export.image.toAsset({
      image: mask.clip(AOI_GEOM), description: 'export_indo_RF_mask_' + yearNum,
      assetId: EXPORT_BASE + '/02_PA/indo_RF_seagrass_mask_' + yearNum,
      region: AOI_GEOM, scale: SCALE, crs: 'EPSG:4326', maxPixels: 1e13
    });

    print('Export queued:', yearNum);
  });

  Map.centerObject(AOI_GEOM, 6);
});
