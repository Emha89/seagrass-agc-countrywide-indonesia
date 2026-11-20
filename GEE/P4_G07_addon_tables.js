/************************************************************
 * P4_G07_addon_tables.js
 * ADD-ON TO P4_G02, P4_G03 (distance area) and P4_G04 (CV), added 2 Nov 2026.
 * Read the header of P4_G06_addon_maps.js first.
 *
 * WHAT THIS SCRIPT DOES
 *   Computes, for the ADD-ON area only (P4_extent_addon_mask, mapped by
 *   P4_G06), the same sums as the core scripts. Every table here has the same
 *   columns as its core table, so P4_07_gee_results.R simply adds the rows:
 *       core (P4_G02/G03/G04) + add-on (this script) = total for the
 *       intended extent (shallow, MaxProb >= 0.7, only class 16 removed).
 *   Not repeated for the add-on area (not needed): samples inside the extent
 *   (computed in R from P4_gee_samples_extent), PA points, GSE similarity
 *   (sample drawn from the core extent; reported as such).
 *
 * EXPORTS
 *   Drive (folder DRIVE_FOLDER):
 *     P4_gee_totals_addon_YYYY.csv     as P4_gee_totals_YYYY (P4_G02)
 *     P4_gee_roi_totals_addon.csv      as P4_gee_roi_totals (P4_G02)
 *     P4_gee_distance_area_addon.csv   as P4_gee_distance_area (P4_G03)
 *     P4_gee_cv_summary_addon.csv      as P4_gee_cv_summary (P4_G04)
 *   Asset: P4_AGC_CV_ADDON_2017_2024     as P4_AGC_CV_2017_2024 (P4_G04)
 *
 * NEEDS: the 8 P4_AGC_COV_ADDON_YYYY assets and P4_extent_addon_mask (P4_G06)
 *
 * NOTE (5 Nov 2026): the CV summary computed on the fly from the 8 add-on maps
 * timed out after 12 h. It now reads the stored CV asset of this script
 * (P4_AGC_CV_ADDON_2017_2024; same mean, sample SD and CV), so run the
 * cvAsset task first and the cv task after it has finished.
 ************************************************************/

// ==========================================================
// 1. CONFIGURATION -- edit here only (paths and classes as in P4_G02-G04)
// ==========================================================
var ASSET_ROOT   = 'projects/YOUR-GEE-PROJECT/assets';   // your Earth Engine asset root
print('P4 script set 2026-11-02 (add-on)');

var OUT_DIR       = ASSET_ROOT + '/chapter_3/output/P4v2';
var TRAIN_TABLE   = OUT_DIR + '/P4_gee_training';
var ADDON_PREFIX  = OUT_DIR + '/P4_AGC_COV_ADDON_';
var ADDON_EXTENT  = OUT_DIR + '/P4_extent_addon_mask';
var ADDON_CV      = OUT_DIR + '/P4_AGC_CV_ADDON_2017_2024';
var MEOW_ASSET    = OUT_DIR + '/P4_meow_indo';
var ROI_ASSET     = ASSET_ROOT + '/chapter_2/InputArea/trainingAGC_COMPLETE_for_GEE_buffer_split3';
var ROI_ID_FIELD  = 'loc_id';
var AOI_ASSET     = ASSET_ROOT + '/packard_fieldData/ecoregion_indo_diss';
var DRIVE_FOLDER  = 'P4v2_GEE';

var YEARS         = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];
var COVER_BREAKS  = [30, 60];                                     // as P4_G02
var SEARCH_RADIUS = 50000;                                        // as P4_G03
var DIST_MAX_ERROR = 5;
var DIST_LIMITS = [10, 50, 100, 250, 500, 1000, 2000, 5000, 20000, SEARCH_RADIUS];
var DIST_LABELS = ['<=10 m', '10-50 m', '50-100 m', '100-250 m', '250-500 m',
                   '500 m-1 km', '1-2 km', '2-5 km', '5-20 km', '20-50 km', '>50 km'];
var CV_BREAKS   = [10, 20, 30, 50, 100];                          // as P4_G04
var TILE_SCALE  = 8;

var RUN = {
  totals:   true,   // 8 tasks
  roi:      true,   // 1 task
  distance: true,   // 1 task (the slowest one)
  cv:       true,   // 1 task (reads the CV asset: run cvAsset first)
  cvAsset:  true    // 1 task (for the CV map figure and the cv task)
};

// true = the CV summary reads P4_AGC_CV_ADDON_2017_2024 (fast);
// false = computes CV from the 8 maps on the fly (timed out on 5 Nov 2026)
var CV_FROM_ASSET = true;


// ==========================================================
// 2. INPUTS AND GRID (the add-on maps share the grid of the core maps)
// ==========================================================
var AOI    = ee.FeatureCollection(AOI_ASSET).geometry();
var extent = ee.Image(ADDON_EXTENT).select(0);
var GRID   = ee.Image(ADDON_PREFIX + YEARS[0]).projection().getInfo();
var GRID_CRS = GRID.crs || GRID.wkt;
function gridArgs(o) { o.crs = GRID_CRS; o.crsTransform = GRID.transform; return o; }

var ecoImg = ee.FeatureCollection(MEOW_ASSET)
  .reduceToImage({properties: ['eco_id'], reducer: ee.Reducer.first()})
  .unmask(0, false).toInt().rename('eco_id');
var area = ee.Image.pixelArea();

function classImage(cov) {
  return ee.Image(1).where(cov.gte(COVER_BREAKS[0]), 2).where(cov.gte(COVER_BREAKS[1]), 3)
    .updateMask(cov.mask()).rename('cover_class').toInt();
}
function sumStack(y) {
  var agc = ee.Image(ADDON_PREFIX + y).select('AGC_pred');
  return agc.multiply(area).rename('agc_g')
    .addBands(area.rename('area_m2'))
    .addBands(ee.Image(1).rename('n_px'))
    .updateMask(agc.mask());
}


// ==========================================================
// 3. TOTALS: year x ecoregion x cover class (as P4_G02)
// ==========================================================
if (RUN.totals) {
  YEARS.forEach(function (y) {
    var cls = classImage(ee.Image(ADDON_PREFIX + y).select('tSPC_pred'));
    var res = sumStack(y).addBands(ecoImg.multiply(10).add(cls).rename('key')).reduceRegion(gridArgs({
      reducer: ee.Reducer.sum().repeat(3).group({groupField: 3, groupName: 'key'}),
      geometry: AOI, maxPixels: 1e13, tileScale: TILE_SCALE
    }));
    var rows = ee.FeatureCollection(ee.List(res.get('groups')).map(function (g) {
      g = ee.Dictionary(g);
      var k = ee.Number(g.get('key')).round();
      var s = ee.List(g.get('sum'));
      return ee.Feature(null, {year: y, eco_id: k.divide(10).floor(), cover_class: k.mod(10),
                               agc_g: s.get(0), area_m2: s.get(1), n_px: s.get(2)});
    }));
    Export.table.toDrive({
      collection: rows, description: 'P4_gee_totals_addon_' + y, folder: DRIVE_FOLDER,
      fileFormat: 'CSV', selectors: ['year', 'eco_id', 'cover_class', 'area_m2', 'agc_g', 'n_px']
    });
  });
  print('Export queued: 8 x P4_gee_totals_addon_YYYY');
}


// ==========================================================
// 4. ROI TOTALS (as P4_G02)
// ==========================================================
if (RUN.roi) {
  var roi = ee.FeatureCollection(ROI_ASSET).select([ROI_ID_FIELD]);
  var roiRows = ee.FeatureCollection(YEARS.map(function (y) {
    return sumStack(y).reduceRegions(gridArgs({
      collection: roi, reducer: ee.Reducer.sum(), tileScale: TILE_SCALE
    })).map(function (f) { return f.set('year', y); });
  })).flatten();
  Export.table.toDrive({
    collection: roiRows, description: 'P4_gee_roi_totals_addon', folder: DRIVE_FOLDER,
    fileFormat: 'CSV', selectors: [ROI_ID_FIELD, 'year', 'area_m2', 'agc_g', 'n_px']
  });
  print('Export queued: P4_gee_roi_totals_addon');
}


// ==========================================================
// 5. DISTANCE TO THE NEAREST AGC SAMPLE (as P4_G03 part C)
// ==========================================================
if (RUN.distance) {
  var agcPts = ee.FeatureCollection(TRAIN_TABLE)
    .map(function (f) {
      return f.setGeometry(ee.Geometry.Point([ee.Number(f.get('xcoord')), ee.Number(f.get('ycoord'))]));
    })
    .filter(ee.Filter.eq('is_agc', 1));

  var dist = agcPts.distance(SEARCH_RADIUS, DIST_MAX_ERROR)
    .unmask(SEARCH_RADIUS * 2, false)
    .rename('dist_m')
    .updateMask(extent);
  var dcls = ee.Image(DIST_LIMITS.length);
  for (var k = DIST_LIMITS.length - 1; k >= 0; k--) dcls = dcls.where(dist.lte(DIST_LIMITS[k]), k);
  dcls = dcls.updateMask(dist.mask()).rename('dist_class').toInt();

  var SUM_NAMES = ['area_m2'];
  var sumImg = area.rename('area_m2');
  YEARS.forEach(function (y) {
    SUM_NAMES.push('agc_g_' + y);
    sumImg = sumImg.addBands(ee.Image(ADDON_PREFIX + y).select('AGC_pred').multiply(area).rename('agc_g_' + y));
  });

  var dres = sumImg.addBands(ecoImg.multiply(100).add(dcls).rename('key')).updateMask(extent)
    .reduceRegion(gridArgs({
      reducer: ee.Reducer.sum().repeat(SUM_NAMES.length).group({groupField: SUM_NAMES.length, groupName: 'key'}),
      geometry: AOI, maxPixels: 1e13, tileScale: TILE_SCALE
    }));
  var distRows = ee.FeatureCollection(ee.List(dres.get('groups')).map(function (g) {
    g = ee.Dictionary(g);
    var kk = ee.Number(g.get('key')).round();
    var c  = kk.mod(100);
    return ee.Feature(null, ee.Dictionary.fromLists(SUM_NAMES, ee.List(g.get('sum'))))
      .set({eco_id: kk.divide(100).floor(), dist_class: c, dist_label: ee.List(DIST_LABELS).get(c)});
  }));
  Export.table.toDrive({
    collection: distRows, description: 'P4_gee_distance_area_addon', folder: DRIVE_FOLDER,
    fileFormat: 'CSV', selectors: ['eco_id', 'dist_class', 'dist_label'].concat(SUM_NAMES)
  });
  print('Export queued: P4_gee_distance_area_addon');
}


// ==========================================================
// 6. AGC CV 2017-2024 (as P4_G04 part A)
// ==========================================================
var agcCol  = ee.ImageCollection(YEARS.map(function (y) { return ee.Image(ADDON_PREFIX + y).select('AGC_pred'); }));
var agcMean = agcCol.mean().rename('AGC_mean');
var agcSd   = agcCol.reduce(ee.Reducer.sampleStdDev()).rename('AGC_sd');
var agcCv   = agcSd.divide(agcMean).multiply(100).rename('AGC_cv');

var cvAssetExists = true;
try { ee.data.getAsset(ADDON_CV); } catch (err) { cvAssetExists = false; }

if (RUN.cv && CV_FROM_ASSET && !cvAssetExists) {
  print('CV summary skipped: ' + ADDON_CV + ' does not exist yet. Run the cvAsset task first.');
} else if (RUN.cv) {
  if (CV_FROM_ASSET) agcCv = ee.Image(ADDON_CV).select('AGC_cv');
  var cvCls = ee.Image(CV_BREAKS.length);
  for (var j = CV_BREAKS.length - 1; j >= 0; j--) cvCls = cvCls.where(agcCv.lt(CV_BREAKS[j]), j);
  cvCls = cvCls.updateMask(agcCv.mask()).toInt();

  var cres = area.rename('area_m2')
    .addBands(agcCv.multiply(area).rename('cv_x_area'))
    .addBands(ecoImg.multiply(10).add(cvCls).rename('key'))
    .updateMask(extent)
    .reduceRegion(gridArgs({
      reducer: ee.Reducer.sum().repeat(2).group({groupField: 2, groupName: 'key'}),
      geometry: AOI, maxPixels: 1e13, tileScale: TILE_SCALE
    }));
  var cvRows = ee.FeatureCollection(ee.List(cres.get('groups')).map(function (g) {
    g = ee.Dictionary(g);
    var kk = ee.Number(g.get('key')).round();
    var s  = ee.List(g.get('sum'));
    return ee.Feature(null, {eco_id: kk.divide(10).floor(), cv_class: kk.mod(10),
                             area_m2: s.get(0), cv_x_area: s.get(1)});
  }));
  Export.table.toDrive({
    collection: cvRows, description: 'P4_gee_cv_summary_addon', folder: DRIVE_FOLDER,
    fileFormat: 'CSV', selectors: ['eco_id', 'cv_class', 'area_m2', 'cv_x_area']
  });
  print('Export queued: P4_gee_cv_summary_addon');
}

if (RUN.cvAsset) {
  Export.image.toAsset({
    image:            agcMean.addBands(agcSd).addBands(agcCv).toFloat().updateMask(extent)
                        .set({years: '2017-2024', part: 'add-on to P4_AGC_CV_2017_2024'}),
    description:      'P4_AGC_CV_ADDON_2017_2024',
    assetId:          ADDON_CV,
    region:           AOI,
    crs:              GRID_CRS,
    crsTransform:     GRID.transform,
    maxPixels:        1e13,
    pyramidingPolicy: {'.default': 'mean'}
  });
  print('Export queued: P4_AGC_CV_ADDON_2017_2024 (asset)');
}
