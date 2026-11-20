/************************************************************
 * P4_G01_train_apply.js
 * Paper 4 (revised design): GSE-only AGC and cover models trained in GEE on
 * the same samples as the R validation, applied to 2017-2024 inside one
 * fixed national extent.
 *
 * RUN ORDER
 *   R   : P4_06_gee_prepare.R -> P4_gee_training.csv, P4_meow_indo.shp
 *   GEE : create the asset folder OUT_DIR (chapter_3/output/P4v2)
 *   GEE : upload into OUT_DIR: P4_gee_training.csv as a table asset
 *         (Assets > New > CSV; X column = xcoord, Y column = ycoord, EPSG:4326)
 *         and P4_meow_indo.shp (+ .shx .dbf .prj .cpg)
 *   GEE : this script -> extent asset, 8 annual map assets, 2 Drive tables
 *   GEE : P4_G02, P4_G03, P4_G04 once the map assets have finished
 *   (P4_G00_extract_gse.js: run once before the final P4_06 when the GSE
 *    check below shows differing points -- it re-extracts the training
 *    vectors from the current dataset version, 1.1 since 2025-11-17)
 *
 * NEW NAMES ONLY: every asset and Drive file of this workflow carries a P4_
 * prefix in chapter_3/output/P4v2 and the Drive folder P4v2_GEE, so nothing
 * from the thesis version (07_AGC, 02_PA, models, TrainingPoint) is touched.
 * An export to an asset ID that already exists fails rather than
 * overwriting; delete that P4v2 asset first to rerun a task.
 *
 * CHANGES FROM 07_AGC_rf_modelDev_indo.js
 *   [1] Predictors: the 64 GSE bands only. No tSPC/AGB/cIndex proxies, no
 *       depth, no distToLand (both lost samples). Same predictors as
 *       P4_02_tuning.R and P4_03_validation.R.
 *   [2] Response: AGC_pred (point estimate), one record per sample. No
 *       bootstrap expansion, no CI-width filter, no log transform.
 *   [3] The training table carries the GSE values used in R, so GEE and R
 *       fit on identical inputs and STEP A (sampling) is no longer needed.
 *       RUN.gseCheck compares the table values with the GSE images at a
 *       random subset of points (same reduceRegion-per-point extraction as
 *       01_trainingData_prep_indo.js).
 *   [4] A second GSE-only model predicts seagrass cover (tSPC, present
 *       samples only). Its classes (<30, 30-59.9, >=60 %, KepMenLH 200/2004)
 *       are the Tier 2 strata.
 *   [5] RF settings from P4_02_tuning.R (P4_best_params.csv).
 *   [6] Maps are masked by the national extent before export:
 *       MaxProb_0_7 mask AND ACA depth 0-500 (unmask(0) pattern) AND
 *       ACA geomorphic != 16. Same mask as 09b. The AGC 1-60 gC m-2 filter
 *       of 09b is dropped.
 *   [7] Exports are written on the grid of the MaxProb mask (crsTransform),
 *       so the extent is not resampled and all later sums use one grid.
 *   [8] Fitted values at the training points (FeatureCollection.classify)
 *       are exported. Their cover class is the map class used for the
 *       Tier 2 stratum means in P4_07_gee_results.R. These are in-sample
 *       values; accuracy comes from the R cross-validation (P4_03/P4_04).
 *
 * NOTE ON ASSET PATHS: the GEE folders "chapter_2" and "chapter_3" are
 * internal labels from an earlier stage and do not match thesis chapters.
 ************************************************************/

// ==========================================================
// 1. CONFIGURATION -- edit here only
// ==========================================================
var ASSET_ROOT   = 'projects/YOUR-GEE-PROJECT/assets';   // your Earth Engine asset root
print('P4 script set 2026-10-31');

// All new assets (inputs and outputs) live in OUT_DIR with a P4_ prefix; the
// existing chapter_2 / chapter_3 assets are only read, never written.
var OUT_DIR      = ASSET_ROOT + '/chapter_3/output/P4v2';      // create this folder first
var TRAIN_TABLE  = OUT_DIR + '/P4_gee_training';               // from P4_06
var MAP_PREFIX   = OUT_DIR + '/P4_AGC_COV_';                    // + year
var EXTENT_ASSET = OUT_DIR + '/P4_extent_mask';
var DRIVE_FOLDER = 'P4v2_GEE';

var MAXPROB_MASK = ASSET_ROOT + '/chapter_3/output/02_PA/indo_Seagrass_MaxProb_0_7_mask';
var DEPTH_ASSET  = ASSET_ROOT + '/chapter_2/InputArea/raster_INA_ACA_Bathymetry';
var ACA_ASSET    = ASSET_ROOT + '/chapter_2/InputArea/acav2_indoarch_stack';
var AOI_ASSET    = ASSET_ROOT + '/packard_fieldData/ecoregion_indo_diss';
var EMB_COLL     = 'GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL';

var TRAIN_YEARS  = [2017, 2018, 2019, 2020, 2021, 2022, 2023];
var APPLY_YEARS  = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];
var SCALE        = 10;
var SEED         = 42;
var NA_VALUE     = -9999;        // missing value code in P4_gee_training.csv

// Expected sample counts from P4_01_build_training.R (printed as a check)
var N_EXPECTED   = {agc: 6616, cover: 7014};   // after dropping 1 sample with GSE = -9999 (31 Oct)

// RF settings from P4_02_tuning.R (P4_best_params.csv, pixel-grouped CV).
// Final run 2 Nov 2026 (GSE version 1.1, 6,616 AGC / 7,014 cover samples):
// AGC 20/5/0.5, cover 12/5/0.9. The grid is flat (RMSE within ~0.2 gC m-2 for
// AGC), so the winner can change with small data changes: always copy the
// values from the latest file.
// ranger min.node.size -> minLeafPopulation, sample.fraction -> bagFraction
// (both sample without replacement when the fraction is < 1).
var RF_AGC = {numberOfTrees: 500, variablesPerSplit: 20, minLeafPopulation: 5,
              bagFraction: 0.5, seed: SEED};
var RF_COV = {numberOfTrees: 500, variablesPerSplit: 12, minLeafPopulation: 5,
              bagFraction: 0.9, seed: SEED};

// Cover classes (KepMenLH 200/2004), right-open as cut(right = FALSE) in R:
// 1: < 30 %, 2: 30-59.9 %, 3: >= 60 %
var COVER_BREAKS = [30, 60];

var RUN = {
  gseCheck: true,    // compare table GSE values with the GSE images (print only)
  extent:   true,    // export the extent mask asset (used by G03 and G04)
  maps:     true,    // export the 8 annual map assets (AGC_pred, tSPC_pred)
  tables:   true     // export importance and fitted values at the points (Drive)
};
var GSE_CHECK_N = 300;


// ==========================================================
// 2. BAND NAMES, GSE IMAGE, AOI, GRID
// ==========================================================
var GSE_BANDS = [];
for (var i = 0; i < 64; i++) GSE_BANDS.push('GSE_A' + (i < 10 ? '0' + i : String(i)));

// Annual GSE mosaic renamed GSE_A00..GSE_A63 (same as scripts 01-07)
function gseImage(year) {
  var start = ee.Date.fromYMD(year, 1, 1);
  return ee.ImageCollection(EMB_COLL)
    .filterDate(start, start.advance(1, 'year'))
    .mosaic()
    .rename(GSE_BANDS);
}

var AOI = ee.FeatureCollection(AOI_ASSET).geometry();

// Values of a reduceRegion(first) result in key order; null -> NA_VALUE
// (a point outside every GSE tile returns null even after unmask)
function safeValues(dict, keys) {
  return ee.List(ee.Dictionary(dict).values(keys)).map(function (x) {
    return ee.Algorithms.If(ee.Algorithms.IsEqual(x, null), NA_VALUE, x);
  });
}

// Grid of the MaxProb mask: every export uses it (see [7])
var GRID = ee.Image(MAXPROB_MASK).projection().getInfo();
var GRID_CRS = GRID.crs || GRID.wkt;
print('Export grid (from the MaxProb mask):', GRID_CRS, GRID.transform);


// ==========================================================
// 3. NATIONAL EXTENT (fixed for 2017-2024)
// ==========================================================
var depthRaw = ee.Image(DEPTH_ASSET).select(0).unmask(0);
var probOK   = ee.Image(MAXPROB_MASK).select('seagrass_presence').unmask(0).gt(0);
var depthOK  = depthRaw.gt(0).and(depthRaw.lte(500));
var geoOK    = ee.Image(ACA_ASSET).select('geomorphic').unmask(16).neq(16);  // no ACA = excluded, as 09b

var extent = probOK.and(depthOK).and(geoOK)
  .selfMask()
  .rename('extent')
  .toByte();

Map.centerObject(AOI, 5);
Map.addLayer(AOI, {color: 'yellow'}, 'AOI (Indonesia)', false);
Map.addLayer(extent, {palette: ['00aa00']}, 'National extent (MaxProb 0.7, depth, geomorphic)', false);


// ==========================================================
// 4. TRAINING DATA (from P4_06_gee_prepare.R)
// ==========================================================
// Point geometry rebuilt from xcoord/ycoord (EPSG:4326), so the point
// locations do not depend on the X/Y columns chosen when the CSV was uploaded
function withPoint(f) {
  return f.setGeometry(ee.Geometry.Point([ee.Number(f.get('xcoord')), ee.Number(f.get('ycoord'))]));
}

var table = ee.FeatureCollection(TRAIN_TABLE).map(withPoint);

var agcTrain = table
  .filter(ee.Filter.eq('is_agc', 1))
  .filter(ee.Filter.gt('AGC_pred', 0))
  .filter(ee.Filter.neq('GSE_A00', NA_VALUE)); // no GSE fill codes (P4_01 drops them too)

var covTrain = table
  .filter(ee.Filter.eq('is_cov', 1))
  .filter(ee.Filter.gte('tSPC', 0))            // drops the NA code (-9999)
  .filter(ee.Filter.neq('GSE_A00', NA_VALUE));

var raw1 = ee.Feature(ee.FeatureCollection(TRAIN_TABLE).first());
print('First row: geometry as uploaded', raw1.geometry(),
      'xcoord, ycoord', raw1.get('xcoord'), raw1.get('ycoord'));

print('Columns of the training table:', table.first().propertyNames());
// Same numbers as the FINGERPRINT printed by P4_06. Any difference = the asset
// is an older upload: delete it and upload the latest P4_gee_training.csv.
print('FINGERPRINT (compare with P4_06):', ee.Dictionary({
  rows:     table.size(),
  agc:      table.filter(ee.Filter.eq('is_agc', 1)).size(),
  cover:    table.filter(ee.Filter.eq('is_cov', 1)).size(),
  max_tSPC: table.filter(ee.Filter.eq('is_cov', 1)).aggregate_max('tSPC'),
  sum_AGC:  table.filter(ee.Filter.eq('is_agc', 1)).aggregate_sum('AGC_pred'),
  gse_fill_rows: table.filter(ee.Filter.eq('GSE_A00', NA_VALUE)).size()   // expect 0
}));
print('AGC samples (R: ' + N_EXPECTED.agc + ')', agcTrain.size());
print('Cover samples (R: ' + N_EXPECTED.cover + ')', covTrain.size());
print('AGC samples per year:', agcTrain.aggregate_histogram('year'));
print('tSPC summary (%, max must be 100: capped in P4_01):', covTrain.aggregate_stats('tSPC'));
print('AGC_pred summary (gC m-2):', agcTrain.aggregate_stats('AGC_pred'));

Map.addLayer(agcTrain.limit(2000), {color: '08306b'}, 'AGC samples (first 2000)', false);


// ==========================================================
// 5. CHECK: table GSE values vs GSE images (print only)
// Same extraction as 01_trainingData_prep_indo.js: reduceRegion(first)
// at each point, scale 10. max_abs_diff ~ 1e-6 means the table and the
// images agree; values near 1 for the dot product mean the same pixel.
// ==========================================================
if (RUN.gseCheck) {
  var sub = agcTrain.randomColumn('rnd', SEED).sort('rnd').limit(GSE_CHECK_N);

  var checked = ee.FeatureCollection(TRAIN_YEARS.map(function (y) {
    var img = gseImage(y).unmask(NA_VALUE, false);      // false: fill outside the tile footprints too
    return sub.filter(ee.Filter.eq('year', y)).map(function (f) {
      var v = img.reduceRegion({
        reducer: ee.Reducer.first(), geometry: f.geometry(), scale: SCALE, maxPixels: 1e13
      });
      var vals = safeValues(v, GSE_BANDS);
      var a = ee.Array(f.toDictionary(GSE_BANDS).values(GSE_BANDS));
      var b = ee.Array(vals);
      return ee.Feature(null, {
        year: y,
        image_na: ee.Number(vals.get(0)).eq(NA_VALUE),
        max_abs_diff: a.subtract(b).abs().reduce(ee.Reducer.max(), [0]).get([0]),
        dot: a.multiply(b).reduce(ee.Reducer.sum(), [0]).get([0])
      });
    });
  })).flatten();

  var okPts = checked.filter(ee.Filter.eq('image_na', 0));
  print('GSE check: points checked', checked.size(), 'with no GSE in the image', checked.size().subtract(okPts.size()));
  print('GSE check: max |table - image| (expect ~1e-6)', okPts.aggregate_stats('max_abs_diff'));
  print('GSE check: dot product table x image (expect ~1)', okPts.aggregate_stats('dot'));
  var diffPts = okPts.filter(ee.Filter.gt('max_abs_diff', 1e-4));
  print('GSE check: points with max |diff| > 1e-4', diffPts.size());
  print('GSE check: those points per year', diffPts.aggregate_histogram('year'));
  print('GSE check: dot product at those points (min, p05, p50)',
        diffPts.reduceColumns(ee.Reducer.percentile([0, 5, 50]), ['dot']));
  // Reading the check: dot >= ~0.9 at the differing points = a neighbouring
  // 10 m pixel (the table was extracted on another pixel grid); the model uses
  // the table values and the maps use the image, so the maps can be exported.
  // dot well below 0.9 = different data (year or location): stop and check.
}


// ==========================================================
// 6. MODELS (final, all samples)
// ==========================================================
var rfAgc = ee.Classifier.smileRandomForest(RF_AGC)
  .setOutputMode('REGRESSION')
  .train({features: agcTrain, classProperty: 'AGC_pred', inputProperties: GSE_BANDS});

var rfCov = ee.Classifier.smileRandomForest(RF_COV)
  .setOutputMode('REGRESSION')
  .train({features: covTrain, classProperty: 'tSPC', inputProperties: GSE_BANDS});

print('AGC model (explain):', rfAgc.explain());
print('Cover model (explain):', rfCov.explain());

function coverClass(x) {
  x = ee.Number(x);
  return ee.Number(ee.Algorithms.If(x.gte(COVER_BREAKS[1]), 3,
                   ee.Algorithms.If(x.gte(COVER_BREAKS[0]), 2, 1)));
}

function importanceFC(rf, model) {
  var imp = ee.Dictionary(ee.Dictionary(rf.explain()).get('importance'));
  return ee.FeatureCollection(imp.keys().map(function (k) {
    return ee.Feature(null, {model: model, band: k, importance: imp.get(k)});
  }));
}


// ==========================================================
// 7. EXPORT: extent mask asset
// ==========================================================
if (RUN.extent) {
  Export.image.toAsset({
    image:            extent,
    description:      'P4_extent_mask',
    assetId:          EXTENT_ASSET,
    region:           AOI,
    crs:              GRID_CRS,
    crsTransform:     GRID.transform,
    maxPixels:        1e13,
    pyramidingPolicy: {'.default': 'mode'}
  });
  print('Export queued: extent mask');
}


// ==========================================================
// 8. EXPORT: annual maps (AGC_pred gC m-2, tSPC_pred %), 2017-2024
// ==========================================================
if (RUN.maps) {
  APPLY_YEARS.forEach(function (y) {
    var emb = gseImage(y).updateMask(extent);
    var img = emb.classify(rfAgc).rename('AGC_pred')
      .addBands(emb.classify(rfCov).rename('tSPC_pred'))
      .toFloat()
      .set({
        year: y,
        predictors: 'GSE_A00-GSE_A63',
        rf_agc: JSON.stringify(RF_AGC),
        rf_cover: JSON.stringify(RF_COV),
        extent: 'MaxProb_0_7 AND depth 0-500 AND geomorphic != 16'
      });

    Export.image.toAsset({
      image:            img,
      description:      'P4_AGC_COV_' + y,
      assetId:          MAP_PREFIX + y,
      region:           AOI,
      crs:              GRID_CRS,
      crsTransform:     GRID.transform,
      maxPixels:        1e13,
      pyramidingPolicy: {'.default': 'mean'}
    });
  });
  print('Export queued: ' + APPLY_YEARS.length + ' annual map assets');

  var preview = gseImage(2023).updateMask(extent);
  Map.addLayer(preview.classify(rfAgc),
    {min: 0, max: 80, palette: ['ffffcc', 'a1dab4', '41b6c4', '2c7fb8', '253494']},
    'AGC 2023 (preview, gC m-2)', false);
  Map.addLayer(preview.classify(rfCov),
    {min: 0, max: 100, palette: ['f7fcb9', 'addd8e', '31a354']},
    'Cover 2023 (preview, %)', false);
}


// ==========================================================
// 9. EXPORT: importance and fitted values at the points (Drive)
// ==========================================================
if (RUN.tables) {
  Export.table.toDrive({
    collection:  importanceFC(rfAgc, 'AGC').merge(importanceFC(rfCov, 'COVER')),
    description: 'P4_gee_importance',
    folder:      DRIVE_FOLDER,
    fileFormat:  'CSV',
    selectors:   ['model', 'band', 'importance']
  });

  var fitted = table
    .classify(rfAgc, 'AGC_map')
    .classify(rfCov, 'tSPC_map')
    .map(function (f) { return f.set('cov_class_map', coverClass(f.get('tSPC_map'))); });

  Export.table.toDrive({
    collection:  fitted,
    description: 'P4_gee_points_fitted',
    folder:      DRIVE_FOLDER,
    fileFormat:  'CSV',
    selectors:   ['rid', 'gee_id', 'year', 'is_agc', 'is_cov', 'AGC_pred', 'tSPC', 'eco_id',
                  'AGC_map', 'tSPC_map', 'cov_class_map']
  });
  print('Export queued: P4_gee_importance, P4_gee_points_fitted (Drive folder ' + DRIVE_FOLDER + ')');
}

print('==============================================');
print('NEXT: when the map assets finish, run P4_G02_totals.js,');
print('P4_G03_extent_coverage.js and P4_G04_temporal.js.');
print('==============================================');
