/************************************************************
 * 03_MORPHO_rf_modelDev_indo.js
 * Leaf morphology model / MORPH3 -- full Indonesia.
 *
 * Adopted 1:1 from the proven chapter_2 (Study 1-era) script
 * (03_MORPHO_rf_modelDev_04022026).
 * Changes from original:
 *   [1] AOI: R2_case_study -> ecoregion_indo_diss (full Indonesia)
 *   [2] Training: indo_training_complete_YYYY (unified, morph3 col already present)
 *   [3] GSE bands: A00..A63 -> GSE_A00..GSE_A63 (match Script 01 output)
 *   [4] Output: chapter_3/output/MORPH3/
 *   [5] RF params from R tuning: ntree=300, mtry=7, nodesize=5, bagFraction=0.9
 *       (top-of-file note above once read mtry=10/bagFraction=0.7 -- the
 *       values actually used are in the inline comment and RF_PARAMS below)
 *   [6] No .reproject() on static layers
 *   [7] Depth mask: depthRaw.gt(0).and(depthRaw.lte(500)) -- raw values assumed
 *       POSITIVE here (see the depth sign-convention note below)
 *   [8] No persistence mask -- confirmed not referenced anywhere in this
 *       script; masking is applied later in the Apps viewer, after all
 *       proxy stages are complete
 *   [9] Training filtered to rows with morph3 label (PA=0 excluded)
 *  [10] UI controls preserved: dropdown (class), slider (threshold), year selector
 *
 * DEPTH SIGN CONVENTION (see also 02_PROB_rf_modelDev_indo.js): this
 * script treats the raw bathymetry asset as POSITIVE and negates it for
 * the predictor band (depthRaw.unmask(1).multiply(-1)) -- matching
 * 02_PROB_rf_modelDev_indo.js's own buildStackForYear(). Both of these
 * apply-stage scripts now agree with each other, but both still
 * disagree with 01_trainingData_prep_indo.js, which states the raw
 * asset is NEGATIVE and negates it to get a POSITIVE predictor value.
 * With two independent deployment scripts agreeing with each other
 * against the training-extraction script, this is worth verifying
 * directly (e.g. print/inspect actual depth values from both sides)
 * before relying on any output built from these stacks. Left unchanged
 * here, since this is what actually produced the existing results.
 ************************************************************/

// ==========================================================
// CONFIGURATION
// ==========================================================
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var RESPONSE = 'morph3';
var K_FOLDS  = 5;
var SCALE    = 10;
var SEED     = 42;

var START_YEAR  = 2017;
var END_YEAR    = 2024;
var TRAIN_YEARS = [2017, 2018, 2019, 2020, 2021, 2022, 2023];

var EXPORT_BASE  = ASSET_ROOT + '/chapter_3/output';
var TRAIN_PREFIX = ASSET_ROOT + '/chapter_3/TrainingPoint/indo_trainingAGC_complete_';

var USE_SLOPE_RUGOSITY = false;
var USE_WAVE           = false;
// distToLand REMOVED -- incomplete coverage causes sample loss

// [5] RF hyperparameters from R tuning (gee_best_hyperparameters.csv)
// leafMorpho | Classifier | ntree=300 | mtry=7 | nodesize=5 | bagFraction=0.9
// Best CV accuracy: 0.8866
var RF_PARAMS = {
  numberOfTrees:     300,
  variablesPerSplit: 7,
  minLeafPopulation:  5,
  bagFraction:       0.9,
  seed:              SEED
};

var APPLY_YEARS = ee.List.sequence(START_YEAR, END_YEAR);

// ==========================================================
// [10] VISUALIZATION CONTROLS
// ==========================================================
var DEFAULT_CLASS_TO_VIEW = 'P_mono_Ea';
var DEFAULT_THR           = 0.7;
var SHOW_UI               = true;
var PROB_VIS              = {min: 0, max: 1, palette: ['red', 'yellow', 'green']};


// ==========================================================
// PREDICTOR BANDS -- [3] GSE_A00..GSE_A63
// ==========================================================
var GSE_BANDS = ee.List.sequence(0, 63).map(function(i) {
  return ee.String('GSE_A').cat(ee.Number(i).format('%02d'));
});

var STATIC_BANDS = ['depth'];
if (USE_SLOPE_RUGOSITY) STATIC_BANDS = STATIC_BANDS.concat(['slope', 'rugosity']);
if (USE_WAVE)           STATIC_BANDS = STATIC_BANDS.concat(['wElevation', 'wHeight', 'wPeriod']);

var FEATURE_BANDS = GSE_BANDS.cat(STATIC_BANDS);
print('MORPH3 predictors:', FEATURE_BANDS);


// ==========================================================
// AOI -- [1] Full Indonesia
// ==========================================================
var AOI_FC   = ee.FeatureCollection(
  ASSET_ROOT + '/packard_fieldData/ecoregion_indo_diss'
);
var AOI_GEOM = AOI_FC.geometry();
Map.addLayer(AOI_GEOM, {color: 'yellow'}, 'AOI (Full Indonesia)', false);
Map.centerObject(AOI_GEOM, 5);


// ==========================================================
// LOAD TRAINING DATA
// [2] Unified asset -- morph3 column already present
// [9] Filter to rows with morph3 label (PA=0 have null morph3)
// ==========================================================
var training_all = ee.FeatureCollection([]);
TRAIN_YEARS.forEach(function(y) {
  training_all = training_all.merge(ee.FeatureCollection(TRAIN_PREFIX + y));
});
print('Raw training size (all rows):', training_all.size());

var totalRaw  = training_all.size();
training_all  = training_all.filter(ee.Filter.notNull([RESPONSE]));
print('With morph3 label:', training_all.size());
print('Removed (PA=0, no morph3):', totalRaw.subtract(training_all.size()));

function notNullFor(list) {
  return ee.Filter.and.apply(null, ee.List(list).map(function(b) {
    return ee.Filter.notNull([b]);
  }));
}
var totalBefore = training_all.size();
training_all    = training_all.filter(notNullFor(FEATURE_BANDS.add(RESPONSE)));
print('Rows removed (null predictors):', totalBefore.subtract(training_all.size()));
print('Final MORPH3 training size:', training_all.size());
print('Class distribution:', training_all.aggregate_histogram(RESPONSE));

// Training QA -- visual check
//Map.addLayer(training_all, {color: '00ccff'}, 'MORPH3 training points', true);

// Show only a 2000-point sample, hidden by default
Map.addLayer(
  training_all.limit(2000),
  {color: '00ccff'},
  'MORPH3 training points (sample 2000)',
  false   // hidden by default
);

//print('Per-year distribution:', training_all.aggregate_histogram('year'));


// ==========================================================
// CLASS ENCODING (string -> integer for RF stability)
// ==========================================================
var classList   = ee.List(training_all.aggregate_array(RESPONSE)).distinct().sort();
print('morph3 classes (sorted):', classList);

var idxList     = ee.List.sequence(0, classList.size().subtract(1));
var labelToCode = ee.Dictionary.fromLists(classList, idxList);
var RESPONSE_ID = 'morph3_id';

training_all = training_all.map(function(f) {
  return f.set(RESPONSE_ID, labelToCode.get(f.get(RESPONSE)));
});
print('Encoded morph3 -> morph3_id (0..n-1)');


// ==========================================================
// K-FOLD CROSS VALIDATION
// ==========================================================
var withFold = training_all.randomColumn('rand', SEED).map(function(f) {
  return f.set('fold', ee.Number(f.get('rand')).multiply(K_FOLDS).int());
});
print('Fold distribution:', withFold.aggregate_histogram('fold'));

var cvPredictions = ee.List.sequence(0, K_FOLDS - 1).map(function(k) {
  k = ee.Number(k);
  var train = withFold.filter(ee.Filter.neq('fold', k));
  var valid = withFold.filter(ee.Filter.eq('fold', k));
  var rf    = ee.Classifier.smileRandomForest(RF_PARAMS).train({
    features: train, classProperty: RESPONSE_ID, inputProperties: FEATURE_BANDS
  });
  return valid.classify(rf, 'pred_id');
});

var cvResult = ee.FeatureCollection(cvPredictions).flatten();
var cm       = cvResult.errorMatrix(RESPONSE_ID, 'pred_id');
print('Confusion matrix (K-fold):', cm);
print('Overall accuracy:', cm.accuracy());
print('Kappa:', cm.kappa());
print('Producer accuracy (recall):', cm.producersAccuracy());
print('Consumer accuracy (precision):', cm.consumersAccuracy());


// ==========================================================
// FINAL MODEL
// ==========================================================
var finalRF = ee.Classifier.smileRandomForest(RF_PARAMS).train({
  features: training_all, classProperty: RESPONSE_ID, inputProperties: FEATURE_BANDS
});
print('Final MORPH3 model trained.');


// ==========================================================
// VARIABLE IMPORTANCE
// ==========================================================
var importance = ee.Dictionary(finalRF.explain().get('importance'));
var impFC      = ee.FeatureCollection(
  importance.keys().map(function(k) {
    k = ee.String(k);
    return ee.Feature(null, {feature: k, importance: ee.Number(importance.get(k))});
  })
).sort('importance', false);

print('Feature Importance (Top 10):', impFC.limit(10));
print(ui.Chart.feature.byFeature({
  features: impFC.limit(66), xProperty: 'feature', yProperties: ['importance']
}).setChartType('ColumnChart').setOptions({
  title: 'Feature Importance -- MORPH3',
  hAxis: {title: 'Feature'}, vAxis: {title: 'Importance'},
  legend: {position: 'none'}, colors: ['grey']
}));


// ==========================================================
// MODEL DIAGNOSTICS EXPORT -> chapter_3/models
// ==========================================================
Export.table.toAsset({
  collection:  impFC.limit(66),
  description: 'export_model_importance_MORPH3',
  assetId:     ASSET_ROOT + '/chapter_3/models/importance_MORPH3'
});
Export.table.toAsset({
  collection:  ee.FeatureCollection([ee.Feature(null, {
    model: 'MORPH3', cv_acc: cm.accuracy(), cv_kappa: cm.kappa(),
    n_train: training_all.size()
  })]),
  description: 'export_model_cv_metrics_MORPH3',
  assetId:     ASSET_ROOT + '/chapter_3/models/cv_metrics_MORPH3'
});
print('Model diagnostics export queued.');


// ==========================================================
// STACK BUILDER -- [6] no .reproject(), [7] depth mask
// distToLand REMOVED -- incomplete coverage causes sample loss
// See the depth sign-convention note at the top of this file.
// ==========================================================
function buildStackForYear(y) {
  var emb = ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL')
    .filterDate(ee.Date.fromYMD(y, 1, 1), ee.Date.fromYMD(ee.Number(y).add(1), 1, 1))
    .mosaic()
    .rename(ee.List.sequence(0, 63).map(function(i) {
      return ee.String('GSE_A').cat(ee.Number(i).format('%02d'));
    }));

  var depthRaw       = ee.Image(ASSET_ROOT + '/chapter_2/InputArea/raster_INA_ACA_Bathymetry');
  var depthRawFilled = depthRaw.unmask(0);
  var mask           = depthRawFilled.gt(0).and(depthRawFilled.lte(500));
  var depth          = depthRaw.unmask(1).multiply(-1).rename('depth');

  var stack = emb.addBands(depth);

  if (USE_SLOPE_RUGOSITY) {
    var slope    = ee.Terrain.slope(depth).rename('slope');
    var rugosity = slope.reduceNeighborhood({
      reducer: ee.Reducer.stdDev(), kernel: ee.Kernel.circle(9, 'meters')
    }).unmask(0).rename('rugosity');
    stack = stack.addBands(slope).addBands(rugosity);
  }
  if (USE_WAVE) {
    var wave = ee.Image(ASSET_ROOT + '/chapter_2/InputArea/era5_gebco_multiband_2024_450m_indo')
      .select([0, 1, 2], ['wElevation', 'wHeight', 'wPeriod']);
    stack = stack.addBands(wave);
  }

  return stack.updateMask(mask).clip(AOI_GEOM);
}


// ==========================================================
// [10] UI STATE -- tracks selected class, threshold, year
// ==========================================================
var uiState = {
  selectedBand: DEFAULT_CLASS_TO_VIEW,
  thr:          DEFAULT_THR,
  currentYear:  null,
  probBands:    null
};

function updateMapLayers() {
  if (!uiState.probBands || uiState.currentYear === null) return;

  var band = uiState.selectedBand;
  var thr  = uiState.thr;

  while (Map.layers().length() > 2) {
    Map.layers().remove(Map.layers().get(Map.layers().length() - 1));
  }

  Map.addLayer(
    uiState.probBands.select(band),
    PROB_VIS,
    'Prob ' + band + ' | ' + uiState.currentYear,
    true
  );

  Map.addLayer(
    uiState.probBands.select(band).gte(thr).selfMask(),
    {palette: ['blue']},
    'Mask ' + band + ' >= ' + thr.toFixed(2) + ' | ' + uiState.currentYear,
    false
  );
}


// ==========================================================
// APPLY MODEL + EXPORT + UI CONTROLS
// ==========================================================
APPLY_YEARS.evaluate(function(years) {

  if (SHOW_UI) {
    var bandNames = classList.map(function(c) {
      return ee.String('P_').cat(ee.String(c));
    }).getInfo();

    var selector = ui.Select({
      items:   bandNames,
      value:   DEFAULT_CLASS_TO_VIEW,
      onChange: function(v) {
        uiState.selectedBand = v;
        updateMapLayers();
      }
    });

    var slider = ui.Slider({
      min: 0, max: 1, step: 0.01,
      value: DEFAULT_THR,
      onChange: function(v) {
        uiState.thr = v;
        updateMapLayers();
      }
    });

    var yearSelector = ui.Select({
      items: years.map(function(yy) { return String(yy); }),
      value: String(years[years.length - 1]),
      onChange: function(v) {
        var y = parseInt(v, 10);
        uiState.currentYear = y;
        var stack     = buildStackForYear(y).select(FEATURE_BANDS);
        var probArr   = stack.classify(finalRF.setOutputMode('MULTIPROBABILITY'));
        var probBands = probArr.arrayFlatten([classList])
          .rename(classList.map(function(c) { return ee.String('P_').cat(ee.String(c)); }));
        uiState.probBands = probBands;
        updateMapLayers();
      }
    });

    print('MORPH3 Probability Viewer');
    print('Select class to display:', selector);
    print('Confidence threshold:', slider);
    print('Preview year:', yearSelector);
  }

  years.forEach(function(y) {
    print('Applying MORPH3 for year:', y);

    var stack     = buildStackForYear(y).select(FEATURE_BANDS);
    var probArr   = stack.classify(finalRF.setOutputMode('MULTIPROBABILITY'));
    var probBands = probArr.arrayFlatten([classList])
      .rename(classList.map(function(c) { return ee.String('P_').cat(ee.String(c)); }));
    var predIdx   = probArr.arrayArgmax().arrayGet([0]).toInt().rename('morph3_pred_id');

    if (uiState.currentYear === null) {
      uiState.currentYear = y;
      uiState.probBands   = probBands;
      if (SHOW_UI) updateMapLayers();
    }

    Export.image.toAsset({
      image: probBands.clip(AOI_GEOM), description: 'export_indo_MORPH3_probs2_' + y,
      assetId: EXPORT_BASE + '/03_MORPH3/indo_MORPH3_probs2_' + y,
      region: AOI_GEOM, scale: SCALE, crs: 'EPSG:4326', maxPixels: 1e13
    });
    Export.image.toAsset({
      image: predIdx.clip(AOI_GEOM), description: 'export_indo_MORPH3_pred2_' + y,
      assetId: EXPORT_BASE + '/03_MORPH3/indo_MORPH3_pred2_' + y,
      region: AOI_GEOM, scale: SCALE, crs: 'EPSG:4326', maxPixels: 1e13
    });

    print('Export queued MORPH3:', y);
  });

  Map.centerObject(AOI_GEOM, 5);
});
