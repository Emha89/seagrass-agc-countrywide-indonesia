/************************************************************
 * 02b_seagrassExtent_indo.js
 * Seagrass extent / persistence mask -- full Indonesia.
 *
 * Purpose:
 *   Build a seagrass persistence mask from yearly PA probability
 *   outputs (02_PROB_rf_modelDev_indo.js). This mask is used as
 *   PERSISTENCE_MASK in all downstream scripts (03 MORPH3, 04 SPC,
 *   05 AGB, 06 cIndex, 07 AGC) for display and spatial filtering.
 *
 * Run order:
 *   02_PROB_rf_modelDev_indo.js  -> [wait for all exports to complete]
 *   02b_seagrassExtent_indo.js   -> [wait for exports to complete]
 *   Scripts 03-07
 *
 * Adopted from the proven chapter_2 (Study 1-era) script
 * (02_extent_2017-2025_v2).
 * Changes:
 *   [1] AOI: R2_case_study -> ecoregion_indo_diss (full Indonesia)
 *   [2] PROB_ASSET_PREFIX: chapter_2 -> chapter_3/output/02_PA/
 *   [3] Output: chapter_3/output/02_PA/
 *   [4] Area stats removed (no sub-location features for full Indonesia)
 *
 * TRAINING-YEAR NOTE: the QA check in Section C below samples known
 * seagrass points from 2018-2023 only (6 years) -- matching
 * 02_PROB_rf_modelDev_indo.js's own `available_years`. 01_trainingData_prep_indo.js
 * exports training points for 2017-2023 (7 years), but 2017's export
 * isn't actually loaded when training the PA model, and 2024 was never
 * exported as training points at all. So only 2018-2023 (6 years) of
 * ground truth feeds model training; 2017 and 2024 predictions rely on
 * the trained model applied to that year's GSE embeddings, with no
 * matching training labels of their own.
 *
 * Outputs:
 *   indo_seagrass_persistence_mask_byYearCount  <- used by scripts 03-07
 *   indo_seagrass_presence_count                <- optional, for visualisation
 *   indo_Seagrass_MaxProb_YYYY_mask             <- used by scripts 03-07 (display)
 *   indo_Seagrass_MaxProb_YYYY_maxprob_image    <- optional
 *
 * NAMING TO VERIFY: check that 03_MORPHO_rf_modelDev_indo.js (and 04-07)
 * reference the exact export name produced below
 * (indo_Seagrass_MaxProb_0_7_mask, given MAX_PROB_THRESHOLD=0.7) --
 * the equivalent scripts in the Study 2 repository had a naming
 * mismatch between this extent script's export and what the
 * downstream scripts actually looked for.
 ************************************************************/

// ==========================================================
// CONFIGURATION
// ==========================================================
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var START_YEAR = 2017;
var END_YEAR   = 2024;

// [2] Input: PA probability outputs from 02_PROB_rf_modelDev_indo.js
var PROB_ASSET_PREFIX =
  ASSET_ROOT + '/chapter_3/output/02_PA/indo_RF_probability2_';

// [3] Output folder
var EXPORT_BASE =
  ASSET_ROOT + '/chapter_3/output/02_PA';

// Persistence threshold: minimum years with seagrass to be considered persistent
var PROB_PRESENCE_THRESHOLD = 0.8;  // per-year probability threshold: 0.8
var MIN_YEARS_PRESENT       = 8;    // min years present out of 8 (2017-2024)

// Max probability mask threshold (for display + downstream masking)
var MAX_PROB_THRESHOLD = 0.7;       // original: 0.8
var THRESHOLD_STR      = String(MAX_PROB_THRESHOLD).replace('.', '_');

// [1] Full Indonesia AOI
var AOI_FC   = ee.FeatureCollection(
  ASSET_ROOT + '/packard_fieldData/ecoregion_indo_diss'
);
var AOI_GEOM = AOI_FC.geometry();

// Year list
var yearList = [];
for (var y = START_YEAR; y <= END_YEAR; y++) { yearList.push(y); }

print('Year range:', START_YEAR, '-', END_YEAR);
print('Prob threshold per year:', PROB_PRESENCE_THRESHOLD);
print('Min years for persistence:', MIN_YEARS_PRESENT);
print('Max prob mask threshold:', MAX_PROB_THRESHOLD);

Map.addLayer(AOI_GEOM, {color: 'yellow'}, 'AOI (Full Indonesia)', false);
Map.centerObject(AOI_GEOM, 5);


// ==========================================================
// DISPLAY: YEARLY PA PROBABILITY MAPS (preview, off by default)
// ==========================================================
yearList.forEach(function(y) {
  Map.addLayer(
    ee.Image(PROB_ASSET_PREFIX + y),
    {min: 0, max: 1, palette: ['#d73027', '#fee08b', '#1a9850']},
    'Prob ' + y, false
  );
});


// ==========================================================
// SECTION A: YEAR-COUNT PERSISTENCE MASK
// Pixel is "persistent" if prob >= PROB_PRESENCE_THRESHOLD
// in at least MIN_YEARS_PRESENT years
// ==========================================================
var presenceImages = yearList.map(function(y) {
  return ee.Image(PROB_ASSET_PREFIX + y).gte(PROB_PRESENCE_THRESHOLD);
});

var presenceCount = ee.ImageCollection(presenceImages).sum()
  .rename('presence_count');

Map.addLayer(presenceCount,
  {min: 0, max: yearList.length, palette: ['white', 'yellow', 'green']},
  'Seagrass year count (prob >= ' + PROB_PRESENCE_THRESHOLD + ')', false);

var persistenceMask = presenceCount.gte(MIN_YEARS_PRESENT)
  .rename('seagrass_persistence')
  .selfMask();

Map.addLayer(persistenceMask,
  {palette: ['#1a9850']},
  'Persistence mask (>= ' + MIN_YEARS_PRESENT + ' years)', true);

// Export: persistence mask (primary output for scripts 03-07)
Export.image.toAsset({
  image:       persistenceMask.clip(AOI_GEOM),
  description: 'export_indo_seagrass_persistence_mask',
  assetId:     EXPORT_BASE + '/indo_seagrass_persistence_mask_byYearCount',
  region:      AOI_GEOM,
  scale:       10,
  crs:         'EPSG:4326',
  maxPixels:   1e13
});
print('Export queued: persistence mask (year count)');

// Export: presence count image (optional, for visualisation)
Export.image.toAsset({
  image:       presenceCount.clip(AOI_GEOM),
  description: 'export_indo_seagrass_presence_count',
  assetId:     EXPORT_BASE + '/indo_seagrass_presence_count',
  region:      AOI_GEOM,
  scale:       10,
  crs:         'EPSG:4326',
  maxPixels:   1e13
});
print('Export queued: presence count image');


// ==========================================================
// SECTION B: MAXIMUM PROBABILITY MASK
// Maps areas ever predicted as seagrass across all years.
// This is the recommended PERSISTENCE_MASK for scripts 03-07
// as it is less restrictive than the year-count mask.
// ==========================================================
var probImages = yearList.map(function(y) {
  return ee.Image(PROB_ASSET_PREFIX + y).rename('prob');
});

var maxProbImage = ee.ImageCollection(probImages)
  .reduce(ee.Reducer.max())
  .rename('max_probability');

Map.addLayer(maxProbImage,
  {min: 0, max: 1, palette: ['#d73027', '#fee08b', '#1a9850']},
  'Max probability (2017-' + END_YEAR + ')', false);

var seagrassMask = maxProbImage
  .gte(MAX_PROB_THRESHOLD)
  .rename('seagrass_presence')
  .selfMask();

Map.addLayer(seagrassMask,
  {palette: ['#ffff00'], opacity: 0.5},
  'Seagrass habitat (MaxProb >= ' + MAX_PROB_THRESHOLD + ')', false);

// Export: max prob image
Export.image.toAsset({
  image:       maxProbImage.clip(AOI_GEOM),
  description: 'export_indo_MaxProb_' + THRESHOLD_STR + '_maxprob_image',
  assetId:     EXPORT_BASE + '/indo_Seagrass_MaxProb_' + THRESHOLD_STR + '_maxprob_image',
  region:      AOI_GEOM,
  scale:       10,
  crs:         'EPSG:4326',
  maxPixels:   1e13
});
print('Export queued: max prob image');

// Export: max prob mask (recommended for scripts 03-07 PERSISTENCE_MASK)
Export.image.toAsset({
  image:       seagrassMask.clip(AOI_GEOM),
  description: 'export_indo_MaxProb_' + THRESHOLD_STR + '_mask',
  assetId:     EXPORT_BASE + '/indo_Seagrass_MaxProb_' + THRESHOLD_STR + '_mask',
  region:      AOI_GEOM,
  scale:       10,
  crs:         'EPSG:4326',
  maxPixels:   1e13
});
print('Export queued: max prob mask (use as PERSISTENCE_MASK in scripts 03-07)');


// ==========================================================
// SECTION C: MINIMUM PROBABILITY (QA -- check field sites coverage)
// Verifies that known seagrass field sites always have prob > 0
// ==========================================================
var minProbImage = ee.ImageCollection(probImages)
  .reduce(ee.Reducer.min())
  .rename('min_probability');

Map.addLayer(minProbImage,
  {min: 0, max: 1, palette: ['#d7191c', '#fdae61', '#1a9641']},
  'Min probability (2017-' + END_YEAR + ')', false);

// Sample min prob at training points that have AGB data (known seagrass sites)
var trainingPoints = ee.FeatureCollection([]);
[2018, 2019, 2020, 2021, 2022, 2023].forEach(function(y) {
  trainingPoints = trainingPoints.merge(
    ee.FeatureCollection(
      ASSET_ROOT + '/chapter_3/TrainingPoint/indo_trainingAGC_complete_' + y
    )
  );
});

var seagrassPoints = trainingPoints.filter(ee.Filter.notNull(['AGB_pred']));
print('Known seagrass points (AGB_pred not null):', seagrassPoints.size());

var sampledMinProb = minProbImage.sampleRegions({
  collection: seagrassPoints,
  scale:      10,
  geometries: true
});

var minProbValue   = sampledMinProb.aggregate_min('min_probability');
var belowThreshold = sampledMinProb.filter(
  ee.Filter.lt('min_probability', PROB_PRESENCE_THRESHOLD)
);

print('Min probability across known seagrass sites:', minProbValue);
print('Known seagrass points below threshold (' + PROB_PRESENCE_THRESHOLD + '):', belowThreshold.size());
Map.addLayer(belowThreshold, {color: 'red'}, 'Below-threshold field points', false);


// ==========================================================
// SUMMARY PRINT
// ==========================================================
print('==============================================');
print('SCRIPT 02b COMPLETE -- update PERSISTENCE_MASK in scripts 03-07:');
print('Option A (year-count, strict):');
print('  ' + EXPORT_BASE + '/indo_seagrass_persistence_mask_byYearCount');
print('Option B (max prob, recommended):');
print('  ' + EXPORT_BASE + '/indo_Seagrass_MaxProb_' + THRESHOLD_STR + '_mask');
print('==============================================');
