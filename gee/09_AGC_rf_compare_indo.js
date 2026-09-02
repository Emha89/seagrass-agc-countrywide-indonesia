/************************************************************
 * 09_AGC_rf_compare_indo.js
 * Total AGC per ROI: Method A (pixel-based RF) vs Method B (Tier-1
 * area-based). Revised version.
 *
 * Revisions from the prior version:
 *   - Method A: read directly from the 07_AGC asset (no re-running
 *     RF, no Monte Carlo)
 *   - Method B: mean of field samples x seagrass area
 *   - SD, CI, and their derivatives removed -- the national-scale
 *     uncertainty formula isn't representative at individual ROI/local
 *     scale
 *   - Output: total AGC for A and B, their difference, and percent
 *     difference -- formatted for Table S3
 *
 * Requires export_ROI_collection.js to have been run first, since
 * this script filters that asset by loc_id via TARGET_LOC_ID below.
 ************************************************************/


// ==========================================================
// CONFIG -- edit here only
// ==========================================================
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var TARGET_LOC_ID  = 'ROI-095';
var TARGET_YEAR    = 2017;
var PROB_THRESHOLD = 0.6;
var AGC_VIS_MIN    = 10;
var AGC_VIS_MAX    = 60;
var S2_CLOUD_MAX   = 20;


// ==========================================================
// CONSTANTS
// ==========================================================
var SCALE       = 10;
var AGC_PALETTE = ['#9e0142','#f46d43','#ffffbf','#66c2a5','#3288bd'];


// ==========================================================
// ASSETS
// ==========================================================
var FIELD_PREFIX     = ASSET_ROOT + '/chapter_3/TrainingPoint/indo_trainingAGC_complete_';
var PA_PREFIX        = ASSET_ROOT + '/chapter_3/output/02_PA/indo_RF_probability2_';
var AGC_ASSET_PREFIX = ASSET_ROOT + '/chapter_3/output/07_AGC/indo_RF_AGCsimulated_reduceVar_limitCI2_';
var ROI_ASSET        = ASSET_ROOT + '/chapter_2/InputArea/trainingAGC_COMPLETE_for_GEE_buffer_split3';


// ==========================================================
// MASKS + AOI
// ==========================================================
var probImg = ee.Image(PA_PREFIX + TARGET_YEAR).rename('prob');
var sgMask  = probImg.gte(PROB_THRESHOLD).selfMask();

var AOI = ee.FeatureCollection(ROI_ASSET)
            .filter(ee.Filter.eq('loc_id', TARGET_LOC_ID))
            .first().geometry();

Map.centerObject(AOI, 18);


// ==========================================================
// SENTINEL-2 -- best single scene (lowest cloud cover)
// ==========================================================
var s2Best = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterDate(
    ee.Date.fromYMD(TARGET_YEAR, 1, 1),
    ee.Date.fromYMD(TARGET_YEAR + 1, 1, 1)
  )
  .filterBounds(AOI)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', S2_CLOUD_MAX))
  .sort('CLOUDY_PIXEL_PERCENTAGE')
  .first()
  .select(['B4', 'B3', 'B2'])
  .clip(AOI);

Map.addLayer(
  s2Best,
  {bands: ['B4', 'B3', 'B2'], min: 0, max: 2500},
  'S2 RGB best scene - ' + TARGET_YEAR
);

s2Best.date().evaluate(function(d) {
  print('S2 best scene date: ' + d);
});


// ==========================================================
// AGC ASSET -- Method A map display (pre-computed 07_AGC)
// ==========================================================
var agcAsset = ee.Image(AGC_ASSET_PREFIX + TARGET_YEAR)
  .select(0).rename('AGC_pred')
  .updateMask(sgMask)
  .clip(AOI);

Map.addLayer(
  agcAsset,
  {min: AGC_VIS_MIN, max: AGC_VIS_MAX, palette: AGC_PALETTE},
  'AGC asset - ' + TARGET_YEAR
);


// ==========================================================
// PRINT CONFIG SUMMARY
// ==========================================================
print('=== CONFIG ===');
print('  LOC_ID    : ' + TARGET_LOC_ID);
print('  Year      : ' + TARGET_YEAR);
print('  Threshold : ' + PROB_THRESHOLD);
print('==============');


// ==========================================================
// FIELD POINTS -- METHOD B
// Parse AGC_pred (stored as String in asset) -> numeric AGC_val
// ==========================================================
var fieldRaw = ee.FeatureCollection(FIELD_PREFIX + TARGET_YEAR)
  .filterBounds(AOI)
  .map(function(f) {
    var raw  = f.get('AGC_pred');
    var safe = ee.Algorithms.If(
      ee.Algorithms.IsEqual(raw, 'NA'), null,
      ee.Number.parse(ee.String(raw))
    );
    return f.set('AGC_val', safe);
  })
  .filter(ee.Filter.notNull(['AGC_val']));

// Keep only points within seagrass mask (prob >= threshold)
var fieldFiltered = probImg.sampleRegions({
  collection: fieldRaw, scale: SCALE, geometries: true
}).filter(ee.Filter.gte('prob', PROB_THRESHOLD));


// ==========================================================
// MAP LAYERS -- probability, mask, field points
// ==========================================================
Map.addLayer(
  probImg.clip(AOI),
  {min: 0, max: 1, palette: ['#ffffff', '#00aa55']},
  'Probability - ' + TARGET_YEAR, false
);

Map.addLayer(
  sgMask.clip(AOI),
  {palette: ['#00ff88'], opacity: 0.4},
  'Seagrass mask (prob >= ' + PROB_THRESHOLD + ')'
);

Map.addLayer(
  fieldRaw,
  {color: '#aaaaaa', pointSize: 4},
  'Field pts - all - ' + TARGET_YEAR
);

Map.addLayer(
  fieldFiltered,
  {color: '#ff6f00', pointSize: 5},
  'Field pts - prob >= ' + PROB_THRESHOLD
);


// ==========================================================
// GUARD -- stop early if no field points pass the threshold
// ==========================================================
fieldFiltered.size().evaluate(function(n) {

  print('=== FIELD POINTS CHECK ===');
  print('  Field pts (>= threshold) : ' + n);
  fieldRaw.size().evaluate(function(nAll) {
    print('  Field pts (all in AOI)   : ' + nAll);
  });

  if (n < 1) {
    print('');
    print('STOPPED -- no field points found after threshold filter.');
    print('  Options:');
    print('    1. Lower PROB_THRESHOLD in CONFIG');
    print('    2. Choose a different TARGET_YEAR');
    print('    3. Choose a different TARGET_LOC_ID');
    return;
  }

  // --------------------------------------------------------
  // SEAGRASS AREA (m2)
  // --------------------------------------------------------
  var pixelArea    = ee.Image.pixelArea();
  var seagrassArea = ee.Number(
    sgMask.multiply(pixelArea)
      .reduceRegion({
        reducer:    ee.Reducer.sum(),
        geometry:   AOI,
        scale:      SCALE,
        maxPixels:  1e13,
        bestEffort: true
      }).get('prob')
  );

  // --------------------------------------------------------
  // METHOD A -- pixel-based, from the national 07_AGC asset
  //   Formula: sum(AGC_pixel x pixel_area)  [gC/m2 x m2 = gC] -> ton C
  // --------------------------------------------------------
  var totalA_ton = ee.Number(
    agcAsset.multiply(pixelArea)
      .reduceRegion({
        reducer:    ee.Reducer.sum(),
        geometry:   AOI,
        scale:      SCALE,
        maxPixels:  1e13,
        bestEffort: true
      }).get('AGC_pred')
  ).divide(1e6);   // gC -> ton C (1 ton = 1e6 gC)

  // --------------------------------------------------------
  // METHOD A -- pixel statistics (min, max, mean) from the asset
  // --------------------------------------------------------
  var agcStats = agcAsset.reduceRegion({
    reducer:    ee.Reducer.min().combine(ee.Reducer.max(),  '', true)
                               .combine(ee.Reducer.mean(), '', true),
    geometry:   AOI,
    scale:      SCALE,
    maxPixels:  1e13,
    bestEffort: true
  });

  var pixelMin  = ee.Number(agcStats.get('AGC_pred_min'));
  var pixelMax  = ee.Number(agcStats.get('AGC_pred_max'));
  var pixelMean = ee.Number(agcStats.get('AGC_pred_mean'));

  // --------------------------------------------------------
  // METHOD B -- Tier-1 area-based, from field samples
  //   Formula: mean(AGC_field)  [gC/m2] x seagrass_area [m2] -> ton C
  // --------------------------------------------------------
  var rMean  = fieldFiltered.reduceColumns(ee.Reducer.mean(),  ['AGC_val']);
  var rN     = fieldFiltered.reduceColumns(ee.Reducer.count(), ['AGC_val']);

  var meanB      = ee.Number(rMean.get('mean'));
  var nB         = ee.Number(rN.get('count'));
  var totalB_ton = meanB.multiply(seagrassArea).divide(1e6);

  // --------------------------------------------------------
  // PRINT RESULTS -- waterfall evaluate (avoid null returns)
  // --------------------------------------------------------
  seagrassArea.evaluate(function(sm) {

    totalA_ton.evaluate(function(tA) {

      ee.Dictionary({ min: pixelMin, max: pixelMax, mean: pixelMean })
        .evaluate(function(dStats) {

          ee.Dictionary({ meanB: meanB, nB: nB, tB: totalB_ton })
            .evaluate(function(dB) {

          if (!dB) {
            print('WARNING: Method B reduceColumns failed');
            return;
          }

          var mB   = dB.meanB;
          var n    = dB.nB;
          var tB   = dB.tB;

          var diff = tA - tB;
          var pct  = tB > 0 ? (diff / tB * 100) : null;
          var pStr = pct !== null
                     ? (pct >= 0 ? '+' : '') + pct.toFixed(2) + '%'
                     : 'n/a';
          var flag = pct === null        ? ''
                   : Math.abs(pct) < 5  ? 'consistent (<5%)'
                   : diff > 0           ? 'A > B (pixel > Tier-1)'
                   :                      'A < B (pixel < Tier-1)';

          print('==================================================');
          print('  RESULTS -- ' + TARGET_LOC_ID + ' | Year: ' + TARGET_YEAR);
          print('  Seagrass area : ' + (sm / 10000).toFixed(4) + ' ha');
          print('');
          print('  Method A -- Pixel-based (asset 07_AGC)');
          print('    Pixel min   : ' + (dStats ? dStats.min.toFixed(4)  : 'n/a') + ' gC/m2');
          print('    Pixel max   : ' + (dStats ? dStats.max.toFixed(4)  : 'n/a') + ' gC/m2');
          print('    Pixel mean  : ' + (dStats ? dStats.mean.toFixed(4) : 'n/a') + ' gC/m2');
          print('    Total AGC   : ' + tA.toFixed(4) + ' ton C');
          print('');
          print('  Method B -- Tier-1 area-based (n = ' + n + ')');
          print('    Mean AGC    : ' + mB.toFixed(4) + ' gC/m2');
          print('    Total AGC   : ' + tB.toFixed(4) + ' ton C');
          print('');
          print('  A - B          : ' + diff.toFixed(4) + ' ton C');
          print('  % diff (A-B)/B : ' + pStr + '  ' + flag);
          print('==================================================');

          // Table S3 copy-paste helper
          // Columns: loc_id | year | seagrass_extent_ha | n_sample
          //          | tAGC_pixel_tonC | Mean_AGC_gCm2
          //          | tAGC_area_tonC | difference_tonC | pct_diff
          print('');
          print('-- Table S3 row --------------------------------');
          print(
            TARGET_LOC_ID +
            '\t' + TARGET_YEAR +
            '\t' + (sm / 10000).toFixed(4) +
            '\t' + n +
            '\t' + tA.toFixed(4) +
            '\t' + mB.toFixed(4) +
            '\t' + tB.toFixed(4) +
            '\t' + diff.toFixed(4) +
            '\t' + (pct !== null ? pct.toFixed(2) : 'NA')
          );
          print('Columns: loc_id | year | seagrass_ha | n_sample |' +
                ' tAGC_pixel | Mean_gCm2 | tAGC_area | diff | pct_diff');
          print('--------------------------------------------------');

            }); // end dB.evaluate
        });     // end dStats.evaluate
    });         // end totalA_ton.evaluate
  });           // end seagrassArea.evaluate

}); // end guard
