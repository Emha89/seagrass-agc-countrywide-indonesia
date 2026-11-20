/************************************************************
 * P4_G05_extent_sensitivity.js
 * Paper 4 (revised design): how much the geomorphic mask changes the
 * national extent, and a check of the old 09b area method. Optional, light
 * (one area-only reduction, no RF and no GSE).
 *
 * Geomorphic mask variants (all with MaxProb_0_7 AND ACA depth 0-500):
 *   current : geomorphic != 16, and pixels WITHOUT ACA geomorphic data also
 *             removed. This is what 09b and P4_G01 do: neq(16) on a masked
 *             pixel stays masked, and updateMask then drops that pixel.
 *   only16  : only class 16 removed; pixels without ACA geomorphic kept
 *             (geomorphic.unmask(0, false).neq(16)).
 *   nogeo   : no geomorphic mask.
 *
 * EXPORT (Drive, folder DRIVE_FOLDER)
 *   P4_gee_extent_sensitivity[_30m].csv  eco_id, a_current_m2, a_only16_m2, a_nogeo_m2
 *   At 10 m, a_current summed over eco_id equals the extent area of P4_G02;
 *   at 30 m it is close, and P4_07 uses the ratios between the masks.
 *
 * PRINT (Code Editor, no task): the old 09b area method on the 2023 map
 *   (count of 100 m cells x 1 ha) next to the pixel area at 100 m.
 ************************************************************/

// ==========================================================
// 1. CONFIGURATION -- edit here only
// ==========================================================
var ASSET_ROOT   = 'projects/YOUR-GEE-PROJECT/assets';   // your Earth Engine asset root
print('P4 script set 2026-11-02');

var OUT_DIR      = ASSET_ROOT + '/chapter_3/output/P4v2';
var MAP_PREFIX   = OUT_DIR + '/P4_AGC_COV_';
var MEOW_ASSET   = OUT_DIR + '/P4_meow_indo';
var MAXPROB_MASK = ASSET_ROOT + '/chapter_3/output/02_PA/indo_Seagrass_MaxProb_0_7_mask';
var DEPTH_ASSET  = ASSET_ROOT + '/chapter_2/InputArea/raster_INA_ACA_Bathymetry';
var ACA_ASSET    = ASSET_ROOT + '/chapter_2/InputArea/acav2_indoarch_stack';
var AOI_ASSET    = ASSET_ROOT + '/packard_fieldData/ecoregion_indo_diss';
var DRIVE_FOLDER = 'P4v2_GEE';
var TILE_SCALE   = 8;

// The masks are computed on the fly over the whole AOI, so a 10 m run is slow
// (hours). At 30 m the task is ~9x faster; P4_07 uses the ratios between the
// masks and scales them to the 10 m extent area of P4_G02. null = 10 m grid.
var REDUCE_SCALE = 30;

var RUN = {
  oldMethodPrint: false,  // done on 2 Nov: 8,459 km2 (old method) vs 3,452 km2 (pixel area)
  sensitivity:    true    // 1 task
};


// ==========================================================
// 2. INPUTS
// ==========================================================
var AOI  = ee.FeatureCollection(AOI_ASSET).geometry();
var GRID = ee.Image(MAP_PREFIX + 2023).projection().getInfo();
var GRID_CRS = GRID.crs || GRID.wkt;

var ecoImg = ee.FeatureCollection(MEOW_ASSET)
  .reduceToImage({properties: ['eco_id'], reducer: ee.Reducer.first()})
  .unmask(0, false).toInt().rename('eco_id');

// 0/1 layers, filled everywhere (unmask(x, false)), so no band is ever masked
var probOK  = ee.Image(MAXPROB_MASK).select('seagrass_presence').unmask(0, false).gt(0);
var depth   = ee.Image(DEPTH_ASSET).select(0).unmask(0, false);
var depthOK = depth.gt(0).and(depth.lte(500));
var geo     = ee.Image(ACA_ASSET).select('geomorphic');

var base     = probOK.and(depthOK);                           // nogeo
var vCurrent = base.and(geo.unmask(16, false).neq(16));       // as 09b / P4_G01
var vOnly16  = base.and(geo.unmask(0, false).neq(16));        // only class 16 removed


// ==========================================================
// 3. OLD 09b AREA METHOD vs PIXEL AREA (print)
// ==========================================================
if (RUN.oldMethodPrint) {
  var m = ee.Image(MAP_PREFIX + 2023).select('AGC_pred');
  var args = {geometry: AOI, scale: 100, maxPixels: 1e13, bestEffort: true, tileScale: TILE_SCALE};
  var cnt  = m.reduceRegion(Object.assign({reducer: ee.Reducer.count()}, args));
  var parea = ee.Image.pixelArea().updateMask(m.mask())
    .reduceRegion(Object.assign({reducer: ee.Reducer.sum()}, args));
  print('Old 09b method (100 m cells x 1 ha), km2:', ee.Number(cnt.get('AGC_pred')).multiply(0.01));
  print('Pixel area at 100 m (weighted), km2:', ee.Number(parea.get('area')).divide(1e6));
  print('P4_G02 (10 m pixel area), km2: 3447.65');
}


// ==========================================================
// 4. AREA PER ECOREGION UNDER THE THREE MASKS (1 task)
// ==========================================================
if (RUN.sensitivity) {
  var area  = ee.Image.pixelArea();
  var stack = area.multiply(vCurrent).rename('a_current')
    .addBands(area.multiply(vOnly16).rename('a_only16'))
    .addBands(area.multiply(base).rename('a_nogeo'))
    .addBands(ecoImg)
    .updateMask(base);                                         // largest variant only

  var redArgs = {
    reducer:   ee.Reducer.sum().repeat(3).group({groupField: 3, groupName: 'eco_id'}),
    geometry:  AOI,
    maxPixels: 1e13,
    tileScale: TILE_SCALE
  };
  if (REDUCE_SCALE) { redArgs.crs = 'EPSG:4326'; redArgs.scale = REDUCE_SCALE; }
  else { redArgs.crs = GRID_CRS; redArgs.crsTransform = GRID.transform; }
  var res = stack.reduceRegion(redArgs);

  var rows = ee.FeatureCollection(ee.List(res.get('groups')).map(function (g) {
    g = ee.Dictionary(g);
    var s = ee.List(g.get('sum'));
    return ee.Feature(null, {eco_id: g.get('eco_id'), a_current_m2: s.get(0),
                             a_only16_m2: s.get(1), a_nogeo_m2: s.get(2)});
  }));

  Export.table.toDrive({
    collection:  rows,
    description: 'P4_gee_extent_sensitivity' + (REDUCE_SCALE ? '_' + REDUCE_SCALE + 'm' : ''),
    folder:      DRIVE_FOLDER,
    fileFormat:  'CSV',
    selectors:   ['eco_id', 'a_current_m2', 'a_only16_m2', 'a_nogeo_m2']
  });
  print('Export queued: P4_gee_extent_sensitivity');
}

Map.centerObject(AOI, 5);
Map.addLayer(vCurrent.selfMask(), {palette: ['00aa00']}, 'Extent: current (09b)', false);
Map.addLayer(vOnly16.and(vCurrent.not()).selfMask(), {palette: ['ff8800']}, 'Added if no-ACA pixels kept', false);
Map.addLayer(base.and(vOnly16.not()).selfMask(), {palette: ['ff0000']}, 'Class 16 removed', false);
