/****************************************************
 * 01_trainingData_prep_indo.js
 * Training data prep -- full Indonesia.
 *
 * CONCEPT:
 *   - Input: training_COMPLETE_for_GEE (output of 22_compile_training_GEE.R)
 *     -> already contains ALL label columns: PA, tSPC, morph3,
 *        AGB_pred, AGC_pred, carbon_index, P_morph3, etc.
 *     -> does NOT contain GSE/env predictors (extracted here)
 *
 *   - Extracts GSE + depth at each training point using
 *     reduceRegion() per feature (mirrors the approach used in
 *     the original R-extraction script that produced this
 *     project's earlier, chapter_2-labelled training data
 *     successfully -- see NOTE ON ASSET PATHS below)
 *
 *   - Label columns including NAs are carried through;
 *     each model script (02-07) applies its own notNull filter
 *
 * OUTPUT:
 *   indo_training_complete_YYYY (FeatureCollection asset per year)
 *
 * EXTRACTION METHOD:
 *   Uses point.map + reduceRegion (NOT sampleRegions) to avoid
 *   projection-snapping issues. Each layer is sampled at its
 *   native resolution -- GEE resolves the CRS per-layer as in
 *   the original extraction workflow.
 *   GSE bands are prefixed with 'GSE_' to match the R-side naming.
 *
 * NOTE ON ASSET PATHS: this project's GEE asset folders are named
 * "chapter_2" and "chapter_3" -- these are internal GEE project
 * folder labels from an earlier stage of the work and do not
 * correspond to the thesis chapter numbers. This script's own data
 * is under "chapter_3" despite this being the thesis's Chapter 5
 * pipeline. Replace ASSET_ROOT with your own project path below.
 *
 * NOTE ON ASSET NAME: TRAINING_COMPLETE points to an asset named
 * "trainingAGC_COMPLETE_for_GEE" (the commented-out line above it
 * shows "training_COMPLETE_for_GEE", matching 22_compile_training_GEE.R's
 * literal output filename) -- likely just renamed slightly on
 * upload to GEE; kept as in the original for clarity.
 *
 * OPEN QUESTION: TRAIN_YEARS below covers 2017-2023, not 2017-2024
 * like every R-side script in this repository. Worth confirming
 * whether 2024 is processed separately or was an oversight here.
 ****************************************************/

// ==========================================================
// CONFIGURATION
// ==========================================================
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var SCALE       = 10;
var TRAIN_YEARS = [2017, 2018, 2019, 2020, 2021, 2022, 2023];

// Input: complete training asset from 22_compile_training_GEE.R (uploaded to GEE)
// var TRAINING_COMPLETE = ee.FeatureCollection(
//   ASSET_ROOT + '/chapter_3/TrainingPoint/training_COMPLETE_for_GEE');
var TRAINING_COMPLETE = ee.FeatureCollection(
  ASSET_ROOT + '/chapter_3/TrainingPoint/trainingAGC_COMPLETE_for_GEE'
);

// Full Indonesia AOI
var AOI_FC = ee.FeatureCollection(
  ASSET_ROOT + '/packard_fieldData/ecoregion_indo_diss'
);
var AOI = AOI_FC.geometry();

var EMB_COLL = 'GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL';

// Output folder -- separate from local (5-site) results to avoid overwriting
// var EXPORT_BASE =
//   ASSET_ROOT + '/chapter_3/TrainingPoint/indo_training_complete_';
var EXPORT_BASE =
  ASSET_ROOT + '/chapter_3/TrainingPoint/indo_trainingAGC_complete_';


// ==========================================================
// DIAGNOSTICS
// ==========================================================
print('================================================');
print('TRAINING DATA PREP -- Full Indonesia (unified)');
print('Total points (all years):', TRAINING_COMPLETE.size());
print('PA histogram (all years):', TRAINING_COMPLETE.aggregate_histogram('PA'));
print('Year distribution:', TRAINING_COMPLETE.aggregate_histogram('year'));

// Inspect first feature to confirm year column type (numeric vs string)
var firstFeature = ee.Feature(TRAINING_COMPLETE.first());
print('First feature -- property names:', firstFeature.propertyNames());
print('First feature -- year value:', firstFeature.get('year'));
print('================================================');


// ==========================================================
// LABEL COLUMNS TO CARRY THROUGH extraction
// All label columns from 22_compile_training_GEE.R -- NAs are
// retained here; each model script applies its own notNull filter
// downstream
// ==========================================================
var LABEL_COLS = [
  // Identity
  'gee_id', 'year', 'xcoord', 'ycoord', 'location', 'compositio',
  // Presence-absence
  'PA', 'PA_prob', 'PA_pred',
  // Seagrass percent cover
  'tSPC', 'tSPC_pred',
  // Leaf morphology
  'morph3', 'sg_morpho',
  'P_mixed_long', 'P_mixed_short_plus_mono_short', 'P_mono_Ea',
  // Above-ground biomass
  'AGB_pred', 'AGB_low', 'AGB_up', 'AGB_CIwidt', 'AGB_simulated', 'tAGB_pred',
  // Above-ground carbon
  'AGC_pred', 'AGC_low', 'AGC_up', 'AGC_CIwidt', 'AGC_simulated', 'tAGC_pred',
  // Carbon index
  'carbon_index', 'carbon_index_pred'
];


// ==========================================================
// STATIC PREDICTOR LAYERS
// No .reproject() applied -- each layer is sampled at its
// native CRS/resolution via reduceRegion (mirrors the approach
// used in the original R-extraction script)
// ==========================================================
var depth = ee.Image(
  ASSET_ROOT + '/chapter_2/InputArea/raster_INA_ACA_Bathymetry'
).multiply(-1).rename('depth');
// Note: .multiply(-1) matches the original script convention
// (bathymetry stored as negative, converted to positive depth)

// distToLand REMOVED -- incomplete coverage causes sample loss.
// This is why every R-side model script in this repository trains
// on depth alone (env_vars <- "depth"), not depth + distToLand.


// ==========================================================
// GSE BUILDER -- matches the R-side naming (GSE_ prefix)
// unmask(-9999, false) preserves no-data as -9999 so the feature
// is still emitted (point not dropped)
// ==========================================================
function buildGSE(year) {
  var img = ee.ImageCollection(EMB_COLL)
    .filterDate(
      ee.Date.fromYMD(year, 1, 1),
      ee.Date.fromYMD(ee.Number(year).add(1), 1, 1)
    )
    .filterBounds(AOI)
    .mosaic()
    .unmask(-9999, false);

  // Rename A00..A63 -> GSE_A00..GSE_A63 to match the R-side convention
  return img.rename(
    img.bandNames().map(function(b) {
      return ee.String('GSE_').cat(ee.String(b));
    })
  );
}


// ==========================================================
// EXPORT PER YEAR
// ==========================================================
TRAIN_YEARS.forEach(function(year) {

  var yearStr = String(year);

  // Filter points for this year
  // ee.Filter.or handles both numeric and string year column types
  var pointsYear = TRAINING_COMPLETE.filter(
    ee.Filter.or(
      ee.Filter.eq('year', year),
      ee.Filter.eq('year', yearStr)
    )
  );

  print('----------------------------------');
  print('YEAR:', yearStr);
  print('Total points this year:', pointsYear.size());
  print('PA histogram:', pointsYear.aggregate_histogram('PA'));
  print('morph3 histogram:', pointsYear.aggregate_histogram('morph3'));

  // Build GSE image for this year
  var gseImg = buildGSE(year);

  // Full predictor stack: GSE + depth only
  var stack = gseImg
    .addBands(depth);

  // ----------------------------------------------------------
  // EXTRACTION: reduceRegion per point (mirrors original R-extraction script)
  // Avoids projection-snapping issues caused by sampleRegions + reproject combo.
  // All existing label properties are preserved via point.set(values).
  // ----------------------------------------------------------
  var sampled = pointsYear.map(function(point) {
    var values = stack.reduceRegion({
      reducer:   ee.Reducer.first(),
      geometry:  point.geometry(),
      scale:     SCALE,
      maxPixels: 1e13
    });
    return point.set(values);
  });

  // Diagnostics: count points with valid GSE and depth values
  var validGSE   = sampled.filter(ee.Filter.notNull(['GSE_A00']));
  var validDepth = sampled.filter(ee.Filter.notNull(['depth']));

  print('Sampled (all points):', sampled.size());
  print('  -> with valid GSE_A00:', validGSE.size());
  print('  -> with valid depth:', validDepth.size());

  // Map preview
  Map.addLayer(sampled.filter(ee.Filter.eq('PA', 1)),
    {color: '0066ff'}, 'PA=1 ' + yearStr, false);
  Map.addLayer(sampled.filter(ee.Filter.eq('PA', 0)),
    {color: 'ff0000'}, 'PA=0 ' + yearStr, false);

  // Export to asset
  Export.table.toAsset({
    collection:  sampled,
    description: 'export_indo_training_complete_' + yearStr,
    assetId:     EXPORT_BASE + yearStr
  });

  print('Export queued:', EXPORT_BASE + yearStr);
});

Map.centerObject(AOI, 5);
