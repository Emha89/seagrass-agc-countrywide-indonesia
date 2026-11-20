/************************************************************
 * P4_G08_figure_grid.js
 * Paper 4 (revised design), added 6 Nov 2026: aggregated national grid for
 * the map figures. A national 10 m map cannot be read at the scale of
 * Indonesia (about 5,000 km wide), so the national panels of the figures
 * show sums on a ~1 km grid; 10 m detail is shown in zoom windows (part B,
 * a separate script, after this grid has been used to choose the Banda Sea
 * window).
 *
 * WHAT THIS SCRIPT DOES
 *   Sums 10 m pixels of the FINAL extent (core P4_G01 + add-on P4_G06, i.e.
 *   shallow, MaxProb >= 0.7, only ACA class 16 removed) into cells of
 *   CELL_PX x CELL_PX native pixels (50 x 50 = ~500 m), on the grid of the
 *   map assets, and exports one row per non-empty cell. R aggregates the
 *   cells further for the national panels (e.g. 5-10 km). Same inputs, masks and distance
 *   method as P4_G02, P4_G03, P4_G04 and P4_G07, so the cell sums add up to
 *   the national tables (checked in R: area sum = 4,707.68 km2).
 *   No asset is written; only Drive tables with new names.
 *
 * EXPORTS (Google Drive, folder DRIVE_FOLDER; one row per ~500 m cell;
 * the file names carry the cell size, e.g. P4_grid500m_agcmean)
 *   P4_grid500m_agcmean.csv   lon, lat (cell centre), area_m2, area_core_m2,
 *                             agc_mean_g
 *                             area_m2      = extent area in the cell (core + add-on)
 *                             area_core_m2 = part of it in the core extent
 *                             agc_mean_g   = sum(mean AGC 2017-2024 x pixel area), g C,
 *                                            from the AGC_mean band of the CV assets
 *   P4_grid500m_agc_YYYY.csv  lon, lat, area_m2, agc_g   (one file per year,
 *                             optional, for the annual maps in the Supplement)
 *                             agc_g = sum(AGC_pred x pixel area), g C
 *   P4_grid500m_dist.csv lon, lat, area_d00 ... area_d10
 *                        extent area per distance class to the nearest AGC
 *                        sample (classes as P4_G03: <=10 m ... >50 km)
 *   P4_grid500m_cv.csv   lon, lat, area_cv_m2, cv_x_area
 *                        area-weighted AGC CV 2017-2024 = cv_x_area / area_cv_m2
 *
 * NEEDS: P4_extent_mask, P4_AGC_COV_YYYY, P4_AGC_CV_2017_2024 (P4_G01, G04)
 *        P4_extent_addon_mask, P4_AGC_COV_ADDON_YYYY, P4_AGC_CV_ADDON_2017_2024
 *        (P4_G06, G07), P4_gee_training (P4_06)
 *
 * RUN: agcMean, dist and cv = 3 tasks; agcYears = 8 more tasks (one per year).
 *
 * NOTE (6 Nov 2026): the first version used ~1 km cells (100 x 100 = 10,000
 * native pixels per cell) and all three tasks failed within seconds with
 * "User memory limit exceeded": every tile had to hold the 10 m pixels of all
 * its cells for 10-12 bands at once. Cells of 50 x 50 pixels (4x fewer pixels
 * per cell) and TILE_SCALE 16 (the maximum; 16x smaller tiles) reduce the
 * pixels held per tile by about 64x. If a task still fails for memory, set
 * CELL_PX = 25 (~250 m); the export names change with it.
 *
 * NOTE (6 Nov 2026, version 2026-11-06c): with 500 m cells the dist (7 min)
 * and cv (2 min) tasks finished, but the agc task (8 years in one task, 16
 * map assets read per tile) failed with "Computation timed out" after 21 min.
 * It is now split: agcMean reads the two CV assets only (same load as the cv
 * task), and agcYears gives one task per year (two map assets each).
 * dist and cv need not be rerun: set them to false.
 ************************************************************/

// ==========================================================
// 1. CONFIGURATION -- edit here only (paths and classes as in P4_G02-G07)
// ==========================================================
var ASSET_ROOT   = 'projects/YOUR-GEE-PROJECT/assets';   // your Earth Engine asset root
print('P4 script set 2026-11-06c (figure grid, 500 m cells)');

var OUT_DIR       = ASSET_ROOT + '/chapter_3/output/P4v2';
var TRAIN_TABLE   = OUT_DIR + '/P4_gee_training';
var CORE_PREFIX   = OUT_DIR + '/P4_AGC_COV_';
var ADDON_PREFIX  = OUT_DIR + '/P4_AGC_COV_ADDON_';
var CORE_EXTENT   = OUT_DIR + '/P4_extent_mask';
var ADDON_EXTENT  = OUT_DIR + '/P4_extent_addon_mask';
var CORE_CV       = OUT_DIR + '/P4_AGC_CV_2017_2024';
var ADDON_CV      = OUT_DIR + '/P4_AGC_CV_ADDON_2017_2024';
var AOI_ASSET     = ASSET_ROOT + '/packard_fieldData/ecoregion_indo_diss';
var DRIVE_FOLDER  = 'P4v2_GEE';

var YEARS         = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];

// Cell size in native pixels (50 x 10 m = ~500 m). Memory per tile grows with
// CELL_PX^2: 100 (~1 km) failed for memory on 6 Nov 2026; use 25 if 50 fails.
var CELL_PX       = 50;
var CELL_LABEL    = 'P4_grid' + (CELL_PX * 10) + 'm_';   // export names, e.g. P4_grid500m_agcmean

var SEARCH_RADIUS  = 50000;                                       // as P4_G03
var DIST_MAX_ERROR = 5;
var DIST_LIMITS = [10, 50, 100, 250, 500, 1000, 2000, 5000, 20000, SEARCH_RADIUS];

var TILE_SCALE    = 16;   // maximum; smaller tiles, less memory per tile

var RUN = {
  agcMean:  true,    // 1 task : P4_grid500m_agcmean (main map figure)
  agcYears: true,    // 8 tasks: P4_grid500m_agc_YYYY (annual maps, Supplement)
  dist:     false,   // 1 task : P4_grid500m_dist (finished 6 Nov 2026)
  cv:       false    // 1 task : P4_grid500m_cv   (finished 6 Nov 2026)
};


// ==========================================================
// 2. INPUTS, GRID AND FINAL EXTENT
// ==========================================================
var AOI  = ee.FeatureCollection(AOI_ASSET).geometry();
var GRID = ee.Image(CORE_PREFIX + YEARS[0]).projection().getInfo();
var GRID_CRS = GRID.crs || GRID.wkt;
var NATIVE   = ee.Projection(GRID_CRS, GRID.transform);

// Output cell grid: same origin as the native grid, pixel size x CELL_PX
var T = GRID.transform;   // [xScale, xShear, xOrigin, yShear, yScale, yOrigin]
var CELL_TRANSFORM = [T[0] * CELL_PX, T[1], T[2], T[3], T[4] * CELL_PX, T[5]];
var CELL = ee.Projection(GRID_CRS, CELL_TRANSFORM);
print('Native grid', GRID_CRS, GRID.transform);
print('Cell grid (' + CELL_PX + ' native pixels per side)', CELL_TRANSFORM);

// Final extent = core or add-on (the two masks do not overlap; P4_G06)
var coreExt  = ee.Image(CORE_EXTENT).select(0).unmask(0, false).gt(0);
var addonExt = ee.Image(ADDON_EXTENT).select(0).unmask(0, false).gt(0);
var extentAll = coreExt.or(addonExt).selfMask();

var area = ee.Image.pixelArea();

// Core and add-on images of one asset pair as one image. The mosaic has no
// native projection of its own, so it is set back to the map grid: this keeps
// every reduction below at the 10 m pixels.
function combined(coreId, addonId, band) {
  return ee.ImageCollection([ee.Image(coreId).select(band), ee.Image(addonId).select(band)])
    .mosaic()
    .setDefaultProjection(NATIVE);
}

// Sum 10 m pixels into cells and export one row per non-empty cell
function exportCells(stack, names, description) {
  var cells = stack
    .setDefaultProjection(NATIVE)
    .reduceResolution({reducer: ee.Reducer.sum().unweighted(), maxPixels: 65536})
    .reproject(CELL);
  var rows = cells.addBands(ee.Image.pixelLonLat().rename(['lon', 'lat']))
    .sample({region: AOI, projection: CELL, dropNulls: true, geometries: false,
             tileScale: TILE_SCALE});
  Export.table.toDrive({
    collection:  rows,
    description: description,
    folder:      DRIVE_FOLDER,
    fileFormat:  'CSV',
    selectors:   ['lon', 'lat'].concat(names)
  });
  print('Export queued: ' + description);
}


// ==========================================================
// 3. AREA AND AGC
// Mean 2017-2024: AGC_mean band of the stored CV assets (P4_G04, P4_G07), which
// is the mean of the 8 annual AGC_pred maps, so sum(agc_mean_g) = mean of the
// annual national totals. Per year: one task per year (two map assets each).
// ==========================================================
if (RUN.agcMean) {
  var agcMean = combined(CORE_CV, ADDON_CV, 'AGC_mean');
  var mStack = area.rename('area_m2')
    .addBands(area.multiply(coreExt).rename('area_core_m2'))
    .addBands(agcMean.multiply(area).rename('agc_mean_g'));
  exportCells(mStack.updateMask(extentAll), ['area_m2', 'area_core_m2', 'agc_mean_g'],
              CELL_LABEL + 'agcmean');
}

if (RUN.agcYears) {
  YEARS.forEach(function (y) {
    var yStack = area.rename('area_m2')
      .addBands(combined(CORE_PREFIX + y, ADDON_PREFIX + y, 'AGC_pred').multiply(area).rename('agc_g'));
    exportCells(yStack.updateMask(extentAll), ['area_m2', 'agc_g'], CELL_LABEL + 'agc_' + y);
  });
}


// ==========================================================
// 4. AREA PER DISTANCE CLASS TO THE NEAREST AGC SAMPLE (as P4_G03 part C)
// ==========================================================
if (RUN.dist) {
  var agcPts = ee.FeatureCollection(TRAIN_TABLE)
    .map(function (f) {
      return f.setGeometry(ee.Geometry.Point([ee.Number(f.get('xcoord')), ee.Number(f.get('ycoord'))]));
    })
    .filter(ee.Filter.eq('is_agc', 1));
  print('AGC samples (R: 6616)', agcPts.size());

  var dist = agcPts.distance(SEARCH_RADIUS, DIST_MAX_ERROR)
    .unmask(SEARCH_RADIUS * 2, false)          // beyond the search radius -> last class
    .setDefaultProjection(NATIVE);
  var dcls = ee.Image(DIST_LIMITS.length);
  for (var k = DIST_LIMITS.length - 1; k >= 0; k--) dcls = dcls.where(dist.lte(DIST_LIMITS[k]), k);

  var dNames = [];
  var dStack = null;
  for (var c = 0; c <= DIST_LIMITS.length; c++) {
    var nm = 'area_d' + (c < 10 ? '0' : '') + c;
    dNames.push(nm);
    var b = area.multiply(dcls.eq(c)).rename(nm);
    dStack = dStack === null ? b : dStack.addBands(b);
  }
  exportCells(dStack.updateMask(extentAll), dNames, CELL_LABEL + 'dist');
}


// ==========================================================
// 5. AGC CV 2017-2024 (from the stored CV assets of P4_G04 and P4_G07)
// ==========================================================
if (RUN.cv) {
  var cv = combined(CORE_CV, ADDON_CV, 'AGC_cv');
  var cStack = area.updateMask(cv.mask()).rename('area_cv_m2')
    .addBands(cv.multiply(area).rename('cv_x_area'));
  exportCells(cStack.updateMask(extentAll), ['area_cv_m2', 'cv_x_area'], CELL_LABEL + 'cv');
}


// ==========================================================
// MAP CHECK (display only; nothing here is exported)
// ==========================================================
Map.centerObject(AOI, 5);
Map.addLayer(extentAll, {palette: ['00aa00']}, 'Final extent (core + add-on)', false);
Map.addLayer(addonExt.selfMask(), {palette: ['ff8800']}, 'Add-on part', false);
