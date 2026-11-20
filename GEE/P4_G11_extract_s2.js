/************************************************************
 * P4_G11_extract_s2.js
 * Paper 4, design v3 (14 Nov 2026, revised 14 Nov c): Sentinel-2 predictors
 * at every AGC sample, for the Sentinel-2 benchmark model (validation only,
 * no map).
 *
 * Product level: the best-corrected Sentinel-2 product in Earth Engine, L2A
 * surface reflectance (COPERNICUS/S2_SR_HARMONIZED, Sen2Cor). The benchmark
 * stands for what an operational Sentinel-2 inventory would use, a mission
 * that keeps recording; it is not meant to copy the input level of GSE.
 *
 * Annual mosaic: built the way the GSE annual layer summarises a year
 * (Google Satellite Embedding V1, AlphaEarth Foundations):
 *   - period  : the calendar year of the survey (GSE annual layer = calendar year)
 *   - scenes  : every scene of the year; clouds and shadows are masked per
 *               pixel with Cloud Score+ (cs_cdf >= CS_THRESHOLD) instead of
 *               dropping cloudy scenes, since GSE also draws on all scenes
 *   - passes  : one observation per satellite pass (tiles of the same pass
 *               overlap; AlphaEarth removes these duplicates too)
 *   - bands   : B2, B3, B4, B8, B11 (the Sentinel-2 bands GSE receives)
 *   - summary : percentiles 10, 25, 50, 75 and 90 of each band over the
 *               year (25 predictors, S2_<band>_p<pct>). GSE summarises the
 *               whole annual series; a set of percentiles is the closest
 *               conventional equivalent. p50 is the median composite (the
 *               composite baseline of the AlphaEarth paper); R also fits a
 *               median-only model as a check. p10/p90 instead of p0/p100,
 *               because the extremes keep single residual clouds or glint.
 *   - grid    : 10 m (B11, 20 m, is read at the 10 m sample pixel)
 * Not matched: GSE also uses Landsat, Sentinel-1 and other sources, and may
 * draw on observations just outside the year. No water-column or glint
 * correction is applied (not operational at national scale).
 *
 * L2A gap: Earth Engine warns that its 2017-2018 L2A coverage is not yet
 * global, and 26% of the AGC samples are from 2017-2018. The console prints
 * the L2A and L1C scenes per survey year. Samples without a clear L2A pixel
 * in their survey year leave the Sentinel-2 comparison (has_s2 = FALSE in
 * P4_01, reported in P4_s2_availability.csv); GSE is refitted on the same
 * samples, so the comparison stays equal. L1C is kept as an option for a
 * sensitivity run only (do not mix levels between years).
 *
 * Reads : P4_gee_training (uploaded with P4_06; rid, gee_id, year, xcoord,
 *         ycoord, is_agc). Nothing else is read or written in Assets.
 * Writes: Drive folder DRIVE_FOLDER, one table per survey year:
 *         P4_gee_s2_extract_YYYY.csv  rid, gee_id, year, S2_*, S2_nobs
 *         S2_nobs = clear observations (passes) of B2 in that year at the pixel.
 *         -9999 = no clear observation (R treats it as missing).
 * Next  : copy the CSVs to <project>/data/gee (CFG$s2_dir), run P4_01.
 ************************************************************/

// ==========================================================
// 1. CONFIGURATION -- edit here only
// ==========================================================
var ASSET_ROOT     = 'projects/YOUR-GEE-PROJECT/assets';   // your Earth Engine asset root
var OUT_DIR        = ASSET_ROOT + '/chapter_3/output/P4v2';
var TRAIN_TABLE    = OUT_DIR + '/P4_gee_training';               // from P4_06
var DRIVE_FOLDER   = 'P4v3_GEE';
print('P4 script set 2026-11-14 (v3c)');

var S2_LEVEL        = 'L2A';         // 'L2A' surface reflectance; 'L1C' for a sensitivity run only
var MAX_SCENE_CLOUD = 100;           // 100 = every scene; e.g. 60 = drop scenes with
                                     // CLOUDY_PIXEL_PERCENTAGE >= 60 (thesis Chapter 4)
var USE_CLOUDSCORE  = true;          // per-pixel Cloud Score+ mask
var CS_THRESHOLD    = 0.60;          // keep pixels with cs_cdf >= this
var ONE_PER_PASS    = true;          // merge overlapping tiles of the same pass
var BANDS           = ['B2', 'B3', 'B4', 'B8', 'B11'];
var PCTS            = [10, 25, 50, 75, 90];
var TRAIN_YEARS     = [2017, 2018, 2019, 2020, 2021, 2022, 2023];
var SCALE           = 10;
var NA_VALUE        = -9999;
var AGC_ONLY        = true;          // only samples used by the AGC models
var CHECK_CLEAR     = false;         // true = also count samples with a clear value in the
                                     // console (slow; P4_01 reports the same from the CSVs)


// ==========================================================
// 2. HELPERS
// ==========================================================
var COLL_ID = {L1C: 'COPERNICUS/S2_HARMONIZED', L2A: 'COPERNICUS/S2_SR_HARMONIZED'};
var S2_COLL = COLL_ID[S2_LEVEL];
var CS_COLL = 'GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED';

// Output names in Earth Engine order (band-major: B2_p10 ... B2_p90, B3_p10 ...)
var RAW_NAMES = [], S2_NAMES = [];
BANDS.forEach(function (b) {
  PCTS.forEach(function (p) {
    RAW_NAMES.push(b + '_p' + p);
    S2_NAMES.push('S2_' + b + '_p' + p);
  });
});

function withPoint(f) {
  return f.setGeometry(ee.Geometry.Point([ee.Number(f.get('xcoord')), ee.Number(f.get('ycoord'))]));
}

// Scenes of one calendar year over the sample points
function s2Scenes(year, region, collId) {
  var start = ee.Date.fromYMD(year, 1, 1);
  var col = ee.ImageCollection(collId || S2_COLL)
    .filterDate(start, start.advance(1, 'year'))
    .filterBounds(region);
  if (MAX_SCENE_CLOUD < 100) col = col.filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', MAX_SCENE_CLOUD));
  return col;
}

// One image per satellite pass: tiles of the same pass (spacecraft, date, orbit)
// overlap at their edges and would count one observation two to four times
function onePerPass(col) {
  col = col.map(function (img) {
    var key = ee.String(img.get('SPACECRAFT_NAME')).cat('_')
      .cat(ee.Date(img.get('system:time_start')).format('YYYYMMdd')).cat('_')
      .cat(ee.Number(img.get('SENSING_ORBIT_NUMBER')).format('%d'));
    return img.set('pass_id', key);
  });
  var ids = ee.List(col.aggregate_array('pass_id')).distinct();
  return ee.ImageCollection.fromImages(ids.map(function (id) {
    return ee.Image(col.filter(ee.Filter.eq('pass_id', id)).mosaic()).set('pass_id', id);
  }));
}

function s2Collection(year, region) {
  var col = s2Scenes(year, region);
  if (USE_CLOUDSCORE) {
    col = col.linkCollection(ee.ImageCollection(CS_COLL), ['cs_cdf'])
      .map(function (img) {
        return img.updateMask(img.select('cs_cdf').gte(CS_THRESHOLD));
      });
  }
  col = col.select(BANDS);
  if (ONE_PER_PASS) col = onePerPass(col);
  // A fully masked dummy keeps the band names when a year has no image,
  // so the composite exists and every sample is exported as -9999.
  var dummy = ee.Image.constant(ee.List.repeat(0, BANDS.length)).rename(BANDS)
    .toUint16().updateMask(ee.Image.constant(0));
  return col.merge(ee.ImageCollection([dummy]));
}

function s2Composite(col) {
  var pct  = col.reduce(ee.Reducer.percentile(PCTS)).select(RAW_NAMES).rename(S2_NAMES);
  var nobs = col.select('B2').count().rename('S2_nobs');
  return pct.addBands(nobs).unmask(NA_VALUE, false);
}


// ==========================================================
// 3. EXTRACTION, one export per survey year
// ==========================================================
var table = ee.FeatureCollection(TRAIN_TABLE);
if (AGC_ONLY) table = table.filter(ee.Filter.eq('is_agc', 1));
table = table.map(withPoint);
print('Samples (R: 6616 AGC samples)', table.size());

TRAIN_YEARS.forEach(function (y) {
  var pts = table.filter(ee.Filter.eq('year', y));
  var region = pts.geometry();                      // the sample points, not their bounding box
  var col = s2Collection(y, region);
  var img = s2Composite(col);

  // Console check: L2A and L1C scenes over the sample points in the survey year
  // (L2A clearly below L1C = the L2A gap; those samples leave the comparison)
  print(ee.String('Survey year ' + y + ': samples ').cat(ee.Number(pts.size()).format('%d'))
          .cat(' | L2A scenes ').cat(ee.Number(s2Scenes(y, region, COLL_ID.L2A).size()).format('%d'))
          .cat(' | L1C scenes ').cat(ee.Number(s2Scenes(y, region, COLL_ID.L1C).size()).format('%d')));

  var res = img.reduceRegions({
    collection: pts.select(['rid', 'gee_id', 'year']),
    reducer: ee.Reducer.first(),
    scale: SCALE,
    tileScale: 4
  });

  if (CHECK_CLEAR) {
    var withValue = res.filter(ee.Filter.gt('S2_nobs', 0)).size();
    print(ee.String('  ' + S2_LEVEL + ' ' + y + ': samples with a clear value ')
            .cat(ee.Number(withValue).format('%d')));
  }

  Export.table.toDrive({
    collection:  res,
    description: 'P4_gee_s2_extract_' + y,
    folder:      DRIVE_FOLDER,
    fileFormat:  'CSV',
    selectors:   ['rid', 'gee_id', 'year'].concat(S2_NAMES).concat(['S2_nobs'])
  });
});

print('Exports queued: P4_gee_s2_extract_YYYY -> <project>/data/gee (CFG$s2_dir), then P4_01 (v3)');
print('Predictors (' + S2_NAMES.length + '):', S2_NAMES);
print('If an export runs out of memory or time: set ONE_PER_PASS = false (some passes then count',
      'two to four times at tile edges) or raise tileScale to 8.');
