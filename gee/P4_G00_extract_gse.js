/************************************************************
 * P4_G00_extract_gse.js
 * Paper 4 (revised design): re-extract the 64 GSE values at every training
 * sample from the CURRENT Satellite Embedding collection, so that the
 * training vectors and the maps come from the same dataset version.
 *
 * WHY (31 Oct 2026): the GSE check in P4_G01 found table vectors that differ
 * from the images at 59 of the 2017 points checked (dot 0.95-0.996) and at
 * 2-3 points in other years. The catalog page states that as of 2025-11-17
 * the layers have DATASET_VERSION 1.1, with a regenerated 2017 layer that
 * adds Sentinel-1 acquisitions. master_raw was extracted from the earlier
 * version.
 *
 * RUN ORDER
 *   GEE : this script (needs P4_gee_training, already uploaded) -> 1 Drive table
 *   R   : copy P4_gee_gse_extract.csv to result/P4v2/gee, set CFG$gse_path,
 *         rerun P4_01 -> P4_02 -> P4_06, re-upload P4_gee_training
 *   GEE : P4_G01 (the GSE check should now show ~0 differing points)
 *
 * Extraction as 01_trainingData_prep_indo.js: reduceRegion(first) per point,
 * scale 10, unmask(-9999, false). Points are rebuilt from xcoord/ycoord.
 * Writes only a Drive table; no asset is created or overwritten.
 ************************************************************/

// ==========================================================
// 1. CONFIGURATION -- edit here only
// ==========================================================
var ASSET_ROOT   = 'projects/YOUR-GEE-PROJECT/assets';   // your Earth Engine asset root
print('P4 script set 2026-10-31');

// All new assets (inputs and outputs) live in OUT_DIR with a P4_ prefix; the
// existing chapter_2 / chapter_3 assets are only read, never written.
var OUT_DIR      = ASSET_ROOT + '/chapter_3/output/P4v2';
var TRAIN_TABLE  = OUT_DIR + '/P4_gee_training';               // from P4_06
var EMB_COLL     = 'GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL';
var DRIVE_FOLDER = 'P4v2_GEE';

var TRAIN_YEARS  = [2017, 2018, 2019, 2020, 2021, 2022, 2023];
var SCALE        = 10;
var NA_VALUE     = -9999;


// ==========================================================
// 2. HELPERS
// ==========================================================
var GSE_BANDS = [];
for (var i = 0; i < 64; i++) GSE_BANDS.push('GSE_A' + (i < 10 ? '0' + i : String(i)));

function gseImage(year) {
  var start = ee.Date.fromYMD(year, 1, 1);
  return ee.ImageCollection(EMB_COLL)
    .filterDate(start, start.advance(1, 'year'))
    .mosaic()
    .rename(GSE_BANDS);
}

// Values of a reduceRegion(first) result in key order; null -> NA_VALUE
function safeValues(dict, keys) {
  return ee.List(ee.Dictionary(dict).values(keys)).map(function (x) {
    return ee.Algorithms.If(ee.Algorithms.IsEqual(x, null), NA_VALUE, x);
  });
}

// Point geometry rebuilt from xcoord/ycoord (the CSV upload left it empty)
function withPoint(f) {
  return f.setGeometry(ee.Geometry.Point([ee.Number(f.get('xcoord')), ee.Number(f.get('ycoord'))]));
}


// ==========================================================
// 3. DATASET VERSION (for Methods)
// ==========================================================
TRAIN_YEARS.forEach(function (y) {
  var col = ee.ImageCollection(EMB_COLL)
    .filterDate(ee.Date.fromYMD(y, 1, 1), ee.Date.fromYMD(y + 1, 1, 1));
  print('GSE ' + y + ': images, DATASET_VERSION',
        col.size(), col.aggregate_array('DATASET_VERSION').distinct());
});


// ==========================================================
// 4. EXTRACTION (one reduceRegion per point, survey year)
// ==========================================================
var table = ee.FeatureCollection(TRAIN_TABLE).map(withPoint);

var extracted = ee.FeatureCollection(TRAIN_YEARS.map(function (y) {
  var img = gseImage(y).unmask(NA_VALUE, false);
  return table.filter(ee.Filter.eq('year', y)).map(function (f) {
    var v = img.reduceRegion({
      reducer: ee.Reducer.first(), geometry: f.geometry(), scale: SCALE, maxPixels: 1e13
    });
    return ee.Feature(null, ee.Dictionary.fromLists(GSE_BANDS, safeValues(v, GSE_BANDS)))
      .set({rid: f.get('rid'), gee_id: f.get('gee_id'), year: f.get('year')});
  });
})).flatten();

print('Samples to extract (R: 7014)', table.size());

Export.table.toDrive({
  collection:  extracted,
  description: 'P4_gee_gse_extract',
  folder:      DRIVE_FOLDER,
  fileFormat:  'CSV',
  selectors:   ['rid', 'gee_id', 'year'].concat(GSE_BANDS)
});
print('Export queued: P4_gee_gse_extract -> result/P4v2/gee, then CFG$gse_path in R');
