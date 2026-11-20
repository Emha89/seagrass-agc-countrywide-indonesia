/************************************************************
 * P4_G09_pa_cv.js
 * Paper 4 (revised design): cross-validation accuracy of the seagrass
 * presence (PA) model that made the annual probability layers.
 *
 * WHY
 *   The thesis script 02_PROB (02 -- PA PROBABILITY MODEL, Full Indonesia)
 *   computed a 5-fold cross-validation, but only printed it and queued an
 *   export task (cv_metrics_PA) that was not run. This script repeats the
 *   same cross-validation exactly, without touching the thesis assets:
 *     - same training assets (indo_trainingAGC_complete_2018 ... _2023)
 *     - same filter (all predictors and PA not null) and PA cast to 0/1
 *     - same predictors (64 GSE bands + depth) and the same RF settings
 *       (500 trees, 10 variables per split, min leaf 3, bag fraction 0.8,
 *        seed 42; thesis Table 5.3.2)
 *     - same fold assignment (randomColumn seed 42, 5 folds)
 *   The cross-validation uses the predictor values stored in the training
 *   points, not the GSE images, so it does not depend on the GSE version.
 *   The row at threshold 0.7 should equal the Accuracy, Kappa and F1 that
 *   script 02 printed. Rows at 0.5 and 0.6 are added for reference.
 *   Folds are random, not spatial, so the values describe the model fit
 *   to points like the training points, as in the thesis.
 *
 * EXPORTS (Google Drive, folder DRIVE_FOLDER; new names, nothing overwritten)
 *   P4_gee_pa_cv.csv         threshold, n, n_pres, n_abs, TP, FP, FN, TN,
 *                            accuracy, kappa, precision, recall, f1
 *   P4_gee_pa_cv_counts.csv  year_pa, PA01, n  (training points used)
 *
 * NEEDS: the PA training assets of the thesis (read only)
 * RUN:   paste into the Code Editor, Run, then run both tasks in the Tasks
 *        tab. The console also prints the metrics.
 ************************************************************/

// ==========================================================
// 1. CONFIGURATION -- edit here only
// ==========================================================
var ASSET_ROOT    = 'projects/YOUR-GEE-PROJECT/assets';   // your Earth Engine asset root
print('P4 script set 2026-11-08 (P4_G09)');

var DRIVE_FOLDER    = 'P4v2_GEE';
var PA_TRAIN_PREFIX = ASSET_ROOT + '/chapter_3/TrainingPoint/indo_trainingAGC_complete_';
var TRAIN_YEARS     = [2018, 2019, 2020, 2021, 2022, 2023];   // as script 02 (2017 not used)

var RESPONSE   = 'PA';
var K_FOLDS    = 5;
var SEED       = 42;
var THRESHOLDS = [0.5, 0.6, 0.7];                             // 0.7 = script 02

var STATIC_BANDS = ['depth'];                                 // as script 02 (distToLand removed)
var RF_PARAMS = {
  numberOfTrees:     500,
  variablesPerSplit: 10,
  minLeafPopulation: 3,
  bagFraction:       0.8,
  seed:              SEED
};


// ==========================================================
// 2. TRAINING POINTS (as script 02)
// ==========================================================
var GSE_BANDS = [];
for (var i = 0; i < 64; i++) GSE_BANDS.push('GSE_A' + (i < 10 ? '0' + i : '' + i));
var FEATURE_BANDS = GSE_BANDS.concat(STATIC_BANDS);

function notNullFor(list) {
  return ee.Filter.and.apply(null, list.map(function (b) { return ee.Filter.notNull([b]); }));
}

// Built exactly as in script 02 (merge order, filter, cast), so that the
// fold assignment below gives the same folds.
function castPA(f) {
  return f.set(RESPONSE, ee.Number(ee.Algorithms.If(f.get(RESPONSE), 1, 0)));
}
var training = ee.FeatureCollection([]);
TRAIN_YEARS.forEach(function (y) {
  training = training.merge(ee.FeatureCollection(PA_TRAIN_PREFIX + y));
});
training = training.filter(notNullFor(FEATURE_BANDS.concat([RESPONSE])));
training = training.map(castPA);

print('Training points used:', training.size());
print('PA distribution:', training.aggregate_histogram(RESPONSE));


// ==========================================================
// 3. 5-FOLD CROSS-VALIDATION (fold assignment as script 02)
// ==========================================================
var withFold = training.randomColumn('rand', SEED).map(function (f) {
  return f.set('fold', ee.Number(f.get('rand')).multiply(K_FOLDS).int());
});
print('Fold distribution:', withFold.aggregate_histogram('fold'));

var cvResult = ee.FeatureCollection(ee.List.sequence(0, K_FOLDS - 1).map(function (k) {
  k = ee.Number(k);
  var train = withFold.filter(ee.Filter.neq('fold', k));
  var valid = withFold.filter(ee.Filter.eq('fold', k));
  var rf = ee.Classifier.smileRandomForest(RF_PARAMS).train({
    features: train, classProperty: RESPONSE, inputProperties: FEATURE_BANDS
  });
  return valid.classify(rf.setOutputMode('PROBABILITY'), 'prob');
})).flatten();


// ==========================================================
// 4. METRICS AT EACH THRESHOLD
// ==========================================================
function metricsAt(t) {
  var fc = cvResult.map(function (f) {
    return f.set('pred', ee.Number(f.get('prob')).gte(t).int());
  });
  var cm  = fc.errorMatrix(RESPONSE, 'pred', [0, 1]);
  var arr = ee.Array(cm.array());
  var TN = ee.Number(arr.get([0, 0])), FP = ee.Number(arr.get([0, 1]));
  var FN = ee.Number(arr.get([1, 0])), TP = ee.Number(arr.get([1, 1]));
  var precision = TP.divide(TP.add(FP));
  var recall    = TP.divide(TP.add(FN));
  var f1        = precision.multiply(recall).multiply(2).divide(precision.add(recall));
  return ee.Feature(null, {
    threshold: t, n: TP.add(FP).add(FN).add(TN), n_pres: TP.add(FN), n_abs: TN.add(FP),
    TP: TP, FP: FP, FN: FN, TN: TN,
    accuracy: cm.accuracy(), kappa: cm.kappa(),
    precision: precision, recall: recall, f1: f1
  });
}
var cvTable = ee.FeatureCollection(THRESHOLDS.map(metricsAt));
print('PA cross-validation (row 0.7 = script 02):', cvTable);

// Counts per year from each asset with the same filter (the main collection
// is left unchanged)
var counts = ee.FeatureCollection(TRAIN_YEARS.map(function (y) {
  var d = ee.FeatureCollection(PA_TRAIN_PREFIX + y)
    .filter(notNullFor(FEATURE_BANDS.concat([RESPONSE]))).map(castPA);
  return ee.FeatureCollection([
    ee.Feature(null, {year_pa: y, PA01: 1, n: d.filter(ee.Filter.eq(RESPONSE, 1)).size()}),
    ee.Feature(null, {year_pa: y, PA01: 0, n: d.filter(ee.Filter.eq(RESPONSE, 0)).size()})
  ]);
})).flatten();


// ==========================================================
// 5. EXPORTS
// ==========================================================
Export.table.toDrive({
  collection:  cvTable,
  description: 'P4_gee_pa_cv',
  folder:      DRIVE_FOLDER,
  fileFormat:  'CSV',
  selectors:   ['threshold', 'n', 'n_pres', 'n_abs', 'TP', 'FP', 'FN', 'TN',
                'accuracy', 'kappa', 'precision', 'recall', 'f1']
});
Export.table.toDrive({
  collection:  counts,
  description: 'P4_gee_pa_cv_counts',
  folder:      DRIVE_FOLDER,
  fileFormat:  'CSV',
  selectors:   ['year_pa', 'PA01', 'n']
});
print('Export queued: P4_gee_pa_cv, P4_gee_pa_cv_counts');
