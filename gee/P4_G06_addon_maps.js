/************************************************************
 * P4_G06_addon_maps.js
 * ADD-ON TO P4_G01 (added 2 Nov 2026). Read this before using the P4v2 maps.
 *
 * WHY THIS SCRIPT EXISTS
 *   Extent principle of Paper 4: shallow water (ACA depth 0-500) with
 *   maximum presence probability >= 0.7, minus ACA geomorphic class 16
 *   (removed to limit false seagrass in turbid water, e.g. north of Java).
 *
 *   P4_G01 copied the geomorphic mask of the thesis script 09b:
 *       ee.Image(ACA).select('geomorphic').neq(16)  ->  updateMask(...)
 *   In GEE, neq(16) on a pixel WITHOUT ACA geomorphic data returns no data,
 *   and updateMask then drops that pixel. So P4_G01 also removed every
 *   shallow pixel outside the ACA geomorphic layer (all of Derawan, most of
 *   Malacca Strait and Papua), which the principle did not intend.
 *   P4_G05 measured this at about 34% of the intended extent.
 *
 * WHAT THIS SCRIPT DOES
 *   Maps AGC and cover ONLY for the missing pixels ("add-on" area):
 *       MaxProb >= 0.7 AND depth 0-500 AND geomorphic != 16 (no data kept)
 *       AND NOT P4_extent_mask (the P4_G01 extent, so nothing is mapped twice)
 *   with the same model specification as P4_G01 (same training table, same
 *   RF settings, same seed) on the same grid. P4_G01 assets are not changed.
 *
 * FINAL PRODUCT OF PAPER 4 = CORE (P4_G01) + ADD-ON (this script)
 *   extent : P4_extent_mask OR P4_extent_addon_mask
 *   maps   : ee.ImageCollection([ee.Image(P4_AGC_COV_YYYY),
 *                                ee.Image(P4_AGC_COV_ADDON_YYYY)]).mosaic()
 *   tables : P4_G07_addon_tables.js computes the same sums as P4_G02/G03/G04
 *            for the add-on area; P4_07_gee_results.R adds them to the core
 *            sums (all are sums, so core + add-on = total).
 *
 * RUN ORDER
 *   P4_G01 (done) -> this script -> P4_G07_addon_tables.js -> R: P4_07
 *
 * EXPORTS (assets in OUT_DIR)
 *   P4_extent_addon_mask
 *   P4_AGC_COV_ADDON_2017 ... P4_AGC_COV_ADDON_2024 (bands AGC_pred, tSPC_pred)
 ************************************************************/

// ==========================================================
// 1. CONFIGURATION -- edit here only
// ==========================================================
var ASSET_ROOT   = 'projects/YOUR-GEE-PROJECT/assets';   // your Earth Engine asset root
print('P4 script set 2026-11-02 (add-on)');

var OUT_DIR            = ASSET_ROOT + '/chapter_3/output/P4v2';
var TRAIN_TABLE        = OUT_DIR + '/P4_gee_training';
var CORE_EXTENT_ASSET  = OUT_DIR + '/P4_extent_mask';            // from P4_G01
var ADDON_EXTENT_ASSET = OUT_DIR + '/P4_extent_addon_mask';      // written here
var ADDON_MAP_PREFIX   = OUT_DIR + '/P4_AGC_COV_ADDON_';          // written here, + year

var MAXPROB_MASK = ASSET_ROOT + '/chapter_3/output/02_PA/indo_Seagrass_MaxProb_0_7_mask';
var DEPTH_ASSET  = ASSET_ROOT + '/chapter_2/InputArea/raster_INA_ACA_Bathymetry';
var ACA_ASSET    = ASSET_ROOT + '/chapter_2/InputArea/acav2_indoarch_stack';
var AOI_ASSET    = ASSET_ROOT + '/packard_fieldData/ecoregion_indo_diss';
var EMB_COLL     = 'GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL';

var APPLY_YEARS  = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];
var SEED         = 42;
var NA_VALUE     = -9999;

// MUST be identical to P4_G01 (checked against the asset properties of
// P4_AGC_COV_2017 on 2 Nov 2026: rf_agc 20/5/0.5, rf_cover 12/5/0.9)
var RF_AGC = {numberOfTrees: 500, variablesPerSplit: 20, minLeafPopulation: 5,
              bagFraction: 0.5, seed: SEED};
var RF_COV = {numberOfTrees: 500, variablesPerSplit: 12, minLeafPopulation: 5,
              bagFraction: 0.9, seed: SEED};

var RUN = {
  extent: true,   // 1 task : P4_extent_addon_mask
  maps:   true    // 8 tasks: P4_AGC_COV_ADDON_YYYY
};


// ==========================================================
// 2. HELPERS, AOI, GRID (same as P4_G01)
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

var AOI  = ee.FeatureCollection(AOI_ASSET).geometry();
var GRID = ee.Image(MAXPROB_MASK).projection().getInfo();
var GRID_CRS = GRID.crs || GRID.wkt;


// ==========================================================
// 3. ADD-ON EXTENT
// unmask(x, false): pixels outside an asset's footprint get x as well
// ==========================================================
var probOK  = ee.Image(MAXPROB_MASK).select('seagrass_presence').unmask(0, false).gt(0);
var depth   = ee.Image(DEPTH_ASSET).select(0).unmask(0, false);
var depthOK = depth.gt(0).and(depth.lte(500));
var not16   = ee.Image(ACA_ASSET).select('geomorphic').unmask(0, false).neq(16);  // no data kept
var core    = ee.Image(CORE_EXTENT_ASSET).select(0).unmask(0, false).gt(0);

var addon = probOK.and(depthOK).and(not16).and(core.not())
  .selfMask()
  .rename('extent_addon')
  .toByte();

Map.centerObject(AOI, 5);
Map.addLayer(core.selfMask(), {palette: ['00aa00']}, 'Core extent (P4_G01)', false);
Map.addLayer(addon, {palette: ['ff8800']}, 'Add-on extent (this script)', false);


// ==========================================================
// 4. MODELS (same specification as P4_G01)
// ==========================================================
function withPoint(f) {
  return f.setGeometry(ee.Geometry.Point([ee.Number(f.get('xcoord')), ee.Number(f.get('ycoord'))]));
}
var table = ee.FeatureCollection(TRAIN_TABLE).map(withPoint);

var agcTrain = table
  .filter(ee.Filter.eq('is_agc', 1))
  .filter(ee.Filter.gt('AGC_pred', 0))
  .filter(ee.Filter.neq('GSE_A00', NA_VALUE));
var covTrain = table
  .filter(ee.Filter.eq('is_cov', 1))
  .filter(ee.Filter.gte('tSPC', 0))
  .filter(ee.Filter.neq('GSE_A00', NA_VALUE));

print('AGC samples (P4_G01: 6616)', agcTrain.size());
print('Cover samples (P4_G01: 7014)', covTrain.size());

var rfAgc = ee.Classifier.smileRandomForest(RF_AGC)
  .setOutputMode('REGRESSION')
  .train({features: agcTrain, classProperty: 'AGC_pred', inputProperties: GSE_BANDS});
var rfCov = ee.Classifier.smileRandomForest(RF_COV)
  .setOutputMode('REGRESSION')
  .train({features: covTrain, classProperty: 'tSPC', inputProperties: GSE_BANDS});


// ==========================================================
// 5. EXPORTS
// ==========================================================
if (RUN.extent) {
  Export.image.toAsset({
    image:            addon.set({part: 'add-on to P4_extent_mask', script: 'P4_G06_addon_maps.js'}),
    description:      'P4_extent_addon_mask',
    assetId:          ADDON_EXTENT_ASSET,
    region:           AOI,
    crs:              GRID_CRS,
    crsTransform:     GRID.transform,
    maxPixels:        1e13,
    pyramidingPolicy: {'.default': 'mode'}
  });
  print('Export queued: P4_extent_addon_mask');
}

if (RUN.maps) {
  APPLY_YEARS.forEach(function (y) {
    var emb = gseImage(y).updateMask(addon);
    var img = emb.classify(rfAgc).rename('AGC_pred')
      .addBands(emb.classify(rfCov).rename('tSPC_pred'))
      .toFloat()
      .set({
        year: y,
        part: 'add-on to P4_AGC_COV_' + y,
        script: 'P4_G06_addon_maps.js',
        predictors: 'GSE_A00-GSE_A63',
        rf_agc: JSON.stringify(RF_AGC),
        rf_cover: JSON.stringify(RF_COV),
        extent: 'MaxProb_0_7 AND depth 0-500 AND geomorphic != 16 (no data kept) AND NOT P4_extent_mask'
      });

    Export.image.toAsset({
      image:            img,
      description:      'P4_AGC_COV_ADDON_' + y,
      assetId:          ADDON_MAP_PREFIX + y,
      region:           AOI,
      crs:              GRID_CRS,
      crsTransform:     GRID.transform,
      maxPixels:        1e13,
      pyramidingPolicy: {'.default': 'mean'}
    });
  });
  print('Export queued: ' + APPLY_YEARS.length + ' add-on map assets');
}

print('NEXT: when these assets finish, run P4_G07_addon_tables.js');
