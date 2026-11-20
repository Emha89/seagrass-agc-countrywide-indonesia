/************************************************************
 * P4_G03_extent_coverage.js
 * Paper 4 (revised design): how the field samples and the mapped extent
 * relate to each other.
 *   A. Field -> map : which samples fall inside the national extent, and
 *      which mask drops the others (MaxProb, ACA depth, ACA geomorphic).
 *   B. Extent check : presence probability at the PA training points per
 *      year, for the best-F1 threshold in R. These points trained the PA
 *      model, so F1 here is apparent (in-sample); it checks the threshold,
 *      not the accuracy of the extent.
 *   C. Map -> field : distance from every extent pixel to the nearest AGC
 *      sample (all years pooled), area and pixel-based AGC per distance
 *      class and ecoregion. The class limits include the median
 *      validation distances of the five designs in P4_04 (6 m, 37 m,
 *      460 m, 8.6 km, 396 km), so each part of the extent can be matched
 *      to the validation design that represents it.
 *   The distance classes extend an earlier class set with two classes
 *   beyond 5 km.
 *
 * EXPORTS (Google Drive, folder DRIVE_FOLDER)
 *   P4_gee_samples_extent.csv  rid, gee_id, year, is_agc, is_cov, AGC_pred, eco_id,
 *                              maxprob, depth_raw, geomorphic, in_extent
 *   P4_gee_pa_points.csv       year_pa, PA, PA01, prob
 *   P4_gee_distance_area.csv   eco_id, dist_class, dist_label, area_m2,
 *                              agc_g_2017 ... agc_g_2024
 *   P4_gee_distance_pct.csv    eco_id (-1 = all), p10, p50, p90, mean, count
 *
 * NEEDS: P4_extent_mask and the 8 map assets (P4_G01), P4_meow_indo
 ************************************************************/

// ==========================================================
// 1. CONFIGURATION -- edit here only
// ==========================================================
var ASSET_ROOT    = 'projects/YOUR-GEE-PROJECT/assets';   // your Earth Engine asset root
print('P4 script set 2026-10-31');

// All new assets (inputs and outputs) live in OUT_DIR with a P4_ prefix; the
// existing chapter_2 / chapter_3 assets are only read, never written.
var OUT_DIR       = ASSET_ROOT + '/chapter_3/output/P4v2';
var TRAIN_TABLE   = OUT_DIR + '/P4_gee_training';               // from P4_06
var MAP_PREFIX    = OUT_DIR + '/P4_AGC_COV_';
var EXTENT_ASSET  = OUT_DIR + '/P4_extent_mask';
var MEOW_ASSET    = OUT_DIR + '/P4_meow_indo';
var DEPTH_ASSET   = ASSET_ROOT + '/chapter_2/InputArea/raster_INA_ACA_Bathymetry';
var ACA_ASSET     = ASSET_ROOT + '/chapter_2/InputArea/acav2_indoarch_stack';
var AOI_ASSET     = ASSET_ROOT + '/packard_fieldData/ecoregion_indo_diss';
var DRIVE_FOLDER  = 'P4v2_GEE';

// PA training points and annual probability images (scripts 01 and 02)
var PA_TRAIN_PREFIX = ASSET_ROOT + '/chapter_3/TrainingPoint/indo_trainingAGC_complete_';
var PA_PROB_PREFIX  = ASSET_ROOT + '/chapter_3/output/02_PA/indo_RF_probability2_';
var PA_YEARS        = [2017, 2018, 2019, 2020, 2021, 2022, 2023];

var MAP_YEARS     = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];
var NA_VALUE      = -9999;

// Distances are computed up to SEARCH_RADIUS; pixels farther away fall in
// the last class. Class k covers (LIMITS[k-1], LIMITS[k]].
var SEARCH_RADIUS  = 50000;
var DIST_MAX_ERROR = 5;
var DIST_LIMITS = [10, 50, 100, 250, 500, 1000, 2000, 5000, 20000, SEARCH_RADIUS];
var DIST_LABELS = ['<=10 m', '10-50 m', '50-100 m', '100-250 m', '250-500 m',
                   '500 m-1 km', '1-2 km', '2-5 km', '5-20 km', '20-50 km', '>50 km'];

// null = native 10 m grid of the extent asset (recommended); e.g. 30 if a
// distance task fails for memory
var REDUCE_SCALE = null;
var TILE_SCALE   = 8;

var RUN = {
  samplesInExtent: true,   // A
  paPoints:        true,   // B
  distanceArea:    true,   // C (area and AGC per class)
  distancePct:     false   // C (exact percentiles; optional and very slow, >1 day:
                           //    P4_07 takes the percentile classes from distance_area)
};


// ==========================================================
// 2. INPUTS AND GRID
// ==========================================================
var AOI    = ee.FeatureCollection(AOI_ASSET).geometry();
var extent = ee.Image(EXTENT_ASSET).select(0);
var GRID   = extent.projection().getInfo();
var GRID_CRS = GRID.crs || GRID.wkt;

function gridArgs(o) {
  if (REDUCE_SCALE) { o.crs = 'EPSG:4326'; o.scale = REDUCE_SCALE; }
  else { o.crs = GRID_CRS; o.crsTransform = GRID.transform; }
  return o;
}

var ecoImg = ee.FeatureCollection(MEOW_ASSET)
  .reduceToImage({properties: ['eco_id'], reducer: ee.Reducer.first()})
  .unmask(0, false).toInt().rename('eco_id');        // 0 = outside the MEOW polygons

// Point geometry rebuilt from xcoord/ycoord (EPSG:4326), so the point
// locations do not depend on the X/Y columns chosen when the CSV was uploaded
function withPoint(f) {
  return f.setGeometry(ee.Geometry.Point([ee.Number(f.get('xcoord')), ee.Number(f.get('ycoord'))]));
}

var table  = ee.FeatureCollection(TRAIN_TABLE)
  .map(withPoint)
  .select(['rid', 'gee_id', 'year', 'is_agc', 'is_cov', 'AGC_pred', 'eco_id']);
var agcPts = table.filter(ee.Filter.eq('is_agc', 1));
print('AGC samples (R: 6616)', agcPts.size());


// ==========================================================
// A. FIELD SAMPLES INSIDE THE EXTENT
// ==========================================================
if (RUN.samplesInExtent) {
  // Maximum presence probability 2017-2024, computed here as in 02b
  // (the 02b max-probability image was not exported as an asset; the
  // MaxProb_0_7 mask is this image >= 0.7)
  var maxProb = ee.ImageCollection(MAP_YEARS.map(function (y) {
    return ee.Image(PA_PROB_PREFIX + y).select(0).rename('prob');
  })).max().rename('maxprob');

  // unmask(x, false): points outside an asset's footprint get the fill value too
  var diag = maxProb.unmask(NA_VALUE, false)
    .addBands(ee.Image(DEPTH_ASSET).select(0).rename('depth_raw').unmask(NA_VALUE, false))
    .addBands(ee.Image(ACA_ASSET).select('geomorphic').rename('geomorphic').unmask(NA_VALUE, false))
    .addBands(extent.unmask(0, false).rename('in_extent'));

  var inExt = diag.reduceRegions(gridArgs({
    collection: table,
    reducer:    ee.Reducer.first(),
    tileScale:  4
  }));

  Export.table.toDrive({
    collection:  inExt,
    description: 'P4_gee_samples_extent',
    folder:      DRIVE_FOLDER,
    fileFormat:  'CSV',
    selectors:   ['rid', 'gee_id', 'year', 'is_agc', 'is_cov', 'AGC_pred', 'eco_id',
                  'maxprob', 'depth_raw', 'geomorphic', 'in_extent']
  });
  print('Export queued: P4_gee_samples_extent');
}


// ==========================================================
// B. PRESENCE PROBABILITY AT THE PA TRAINING POINTS
// PA is cast to 0/1 as in 02_PROB_rf_modelDev_indo.js; points without PA
// are dropped first.
// ==========================================================
if (RUN.paPoints) {
  var paRows = ee.FeatureCollection(PA_YEARS.map(function (y) {
    var pts = ee.FeatureCollection(PA_TRAIN_PREFIX + y)
      .filter(ee.Filter.notNull(['PA']))
      .map(function (f) {
        return f.set('PA01', ee.Algorithms.If(f.get('PA'), 1, 0), 'year_pa', y);
      });
    return ee.Image(PA_PROB_PREFIX + y).select(0).rename('prob').unmask(NA_VALUE, false)
      .reduceRegions({collection: pts, reducer: ee.Reducer.first().setOutputs(['prob']),
                       scale: 10, tileScale: 4});   // one band: name the output, else it is 'first'
  })).flatten();

  Export.table.toDrive({
    collection:  paRows,
    description: 'P4_gee_pa_points',
    folder:      DRIVE_FOLDER,
    fileFormat:  'CSV',
    selectors:   ['year_pa', 'PA', 'PA01', 'prob']
  });
  print('Export queued: P4_gee_pa_points');
}


// ==========================================================
// C. DISTANCE FROM EACH EXTENT PIXEL TO THE NEAREST AGC SAMPLE
// ==========================================================
var dist = agcPts.distance(SEARCH_RADIUS, DIST_MAX_ERROR)
  .unmask(SEARCH_RADIUS * 2, false)   // beyond the search radius -> last class (kept in the area table)
  .rename('dist_m')
  .updateMask(extent);

var dcls = ee.Image(DIST_LIMITS.length);
for (var k = DIST_LIMITS.length - 1; k >= 0; k--) {
  dcls = dcls.where(dist.lte(DIST_LIMITS[k]), k);
}
dcls = dcls.updateMask(dist.mask()).rename('dist_class').toInt();

var area = ee.Image.pixelArea();
var SUM_NAMES = ['area_m2'];
var sumImg = area.rename('area_m2');
MAP_YEARS.forEach(function (y) {
  SUM_NAMES.push('agc_g_' + y);
  sumImg = sumImg.addBands(
    ee.Image(MAP_PREFIX + y).select('AGC_pred').multiply(area).rename('agc_g_' + y));
});

if (RUN.distanceArea) {
  var key = ecoImg.multiply(100).add(dcls).rename('key');
  var res = sumImg.addBands(key).updateMask(extent).reduceRegion(gridArgs({
    reducer:   ee.Reducer.sum().repeat(SUM_NAMES.length)
                 .group({groupField: SUM_NAMES.length, groupName: 'key'}),
    geometry:  AOI,
    maxPixels: 1e13,
    tileScale: TILE_SCALE
  }));

  var distRows = ee.FeatureCollection(ee.List(res.get('groups')).map(function (g) {
    g = ee.Dictionary(g);
    var kk = ee.Number(g.get('key')).round();
    var c  = kk.mod(100);
    return ee.Feature(null, ee.Dictionary.fromLists(SUM_NAMES, ee.List(g.get('sum'))))
      .set({eco_id: kk.divide(100).floor(), dist_class: c,
            dist_label: ee.List(DIST_LABELS).get(c)});
  }));

  Export.table.toDrive({
    collection:  distRows,
    description: 'P4_gee_distance_area',
    folder:      DRIVE_FOLDER,
    fileFormat:  'CSV',
    selectors:   ['eco_id', 'dist_class', 'dist_label'].concat(SUM_NAMES)
  });
  print('Export queued: P4_gee_distance_area');
}

if (RUN.distancePct) {
  var pctReducer = ee.Reducer.percentile([10, 50, 90])
    .combine({reducer2: ee.Reducer.mean(), sharedInputs: true})
    .combine({reducer2: ee.Reducer.count(), sharedInputs: true});

  var pctByEco = function (groupImg) {
    var r = dist.addBands(groupImg).reduceRegion(gridArgs({
      reducer:   pctReducer.group({groupField: 1, groupName: 'eco_id'}),
      geometry:  AOI,
      maxPixels: 1e13,
      tileScale: TILE_SCALE
    }));
    return ee.FeatureCollection(ee.List(r.get('groups')).map(function (g) {
      return ee.Feature(null, ee.Dictionary(g));
    }));
  };
  var pctRows = pctByEco(ecoImg)
    .merge(pctByEco(ee.Image(-1).toInt().rename('eco_id')));   // -1 = all ecoregions

  Export.table.toDrive({
    collection:  pctRows,
    description: 'P4_gee_distance_pct',
    folder:      DRIVE_FOLDER,
    fileFormat:  'CSV',
    selectors:   ['eco_id', 'p10', 'p50', 'p90', 'mean', 'count']
  });
  print('Export queued: P4_gee_distance_pct');
}


// ==========================================================
// MAP CHECK
// ==========================================================
Map.centerObject(AOI, 5);
Map.addLayer(extent, {palette: ['00aa00']}, 'Extent', false);
Map.addLayer(dist, {min: 0, max: 5000, palette: ['0571b0', 'f7f7f7', 'ca0020']},
             'Distance to nearest AGC sample (m)', false);
Map.addLayer(agcPts, {color: 'black'}, 'AGC samples', false);
