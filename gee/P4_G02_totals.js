/************************************************************
 * P4_G02_totals.js
 * Paper 4 (revised design): areas and pixel-based AGC sums from the masked
 * annual maps of P4_G01, split by year x MEOW ecoregion x cover class, and
 * per ROI. National, ecoregion and Tier 2 class-area figures all come from
 * the same rows, so every method uses one extent and one grid.
 *
 * Replaces 09b / 08_AGC_rf_modelUncer_indo.js and 09_AGC_rf_compare_indo.js:
 *   - no AGC 1-60 gC m-2 filter
 *   - area from ee.Image.pixelArea() (EPSG:4326 pixels shrink with latitude)
 *   - one grouped reduction per year at the native 10 m grid
 *   - CIs are computed in R (P4_07_gee_results.R): the pixel-based CI uses
 *     the cluster-robust SE of the mean residual from the cluster design
 *   - ROI totals use the national extent, not per-ROI probability thresholds
 *
 * EXPORTS (Google Drive, folder DRIVE_FOLDER)
 *   P4_gee_totals_YYYY.csv   year, eco_id, cover_class, area_m2, agc_g, n_px
 *       agc_g = sum(AGC_pred x pixel area) in g C (/ 1e6 = t C)
 *       eco_id 0 = pixel outside the MEOW polygons; cover_class 1/2/3 =
 *       <30 / 30-59.9 / >=60 % from the cover map (Tier 2 strata)
 *   P4_gee_roi_totals.csv    <ROI_ID_FIELD>, year, area_m2, agc_g, n_px
 *   P4_roi_polygons.geojson  ROI polygons for R (CFG$roi_path)
 *
 * NEEDS: the 8 map assets and P4_meow_indo (uploaded from P4_06_gee_prepare.R)
 ************************************************************/

// ==========================================================
// 1. CONFIGURATION -- edit here only
// ==========================================================
var ASSET_ROOT   = 'projects/YOUR-GEE-PROJECT/assets';   // your Earth Engine asset root
print('P4 script set 2026-10-31');

// All new assets (inputs and outputs) live in OUT_DIR with a P4_ prefix; the
// existing chapter_2 / chapter_3 assets are only read, never written.
var OUT_DIR      = ASSET_ROOT + '/chapter_3/output/P4v2';
var MAP_PREFIX   = OUT_DIR + '/P4_AGC_COV_';
var MEOW_ASSET   = OUT_DIR + '/P4_meow_indo';     // from P4_06 (eco_id, ECOREGION)
var ROI_ASSET    = ASSET_ROOT + '/chapter_2/InputArea/trainingAGC_COMPLETE_for_GEE_buffer_split3';
var ROI_ID_FIELD = 'loc_id';
var ROI_IDS      = [];          // empty = every ROI in ROI_ASSET; or e.g. ['ROI-095', 'ROI-101']
var AOI_ASSET    = ASSET_ROOT + '/packard_fieldData/ecoregion_indo_diss';
var DRIVE_FOLDER = 'P4v2_GEE';

var YEARS        = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];
var COVER_BREAKS = [30, 60];    // same classes as P4_G01 and R (cut, right = FALSE)

// null = native 10 m grid of the map assets (recommended).
// Set e.g. 30 only if a task fails for memory; the sums then use 30 m pixels.
var REDUCE_SCALE = null;
var TILE_SCALE   = 4;

var RUN = {
  national:    true,   // 8 tasks: P4_gee_totals_YYYY
  roi:         true,   // 1 task : P4_gee_roi_totals
  roiPolygons: true    // 1 task : P4_roi_polygons (GeoJSON)
};


// ==========================================================
// 2. INPUTS AND GRID
// ==========================================================
var AOI  = ee.FeatureCollection(AOI_ASSET).geometry();
var GRID = ee.Image(MAP_PREFIX + YEARS[0]).projection().getInfo();
var GRID_CRS = GRID.crs || GRID.wkt;

function gridArgs(o) {
  if (REDUCE_SCALE) { o.crs = 'EPSG:4326'; o.scale = REDUCE_SCALE; }
  else { o.crs = GRID_CRS; o.crsTransform = GRID.transform; }
  return o;
}

var meow = ee.FeatureCollection(MEOW_ASSET);
var ecoImg = meow
  .reduceToImage({properties: ['eco_id'], reducer: ee.Reducer.first()})
  .unmask(0, false)          // 0 = outside the MEOW polygons; such pixels stay in the totals
  .toInt()
  .rename('eco_id');

var roi = ee.FeatureCollection(ROI_ASSET);
if (ROI_IDS.length > 0) roi = roi.filter(ee.Filter.inList(ROI_ID_FIELD, ROI_IDS));
roi = roi.select([ROI_ID_FIELD]);

print('MEOW ecoregions (eco_id, ECOREGION):',
      meow.sort('eco_id').reduceColumns(ee.Reducer.toList(2), ['eco_id', 'ECOREGION']).get('list'));
print('ROIs:', roi.size());


// ==========================================================
// 3. PER-PIXEL BANDS
// ==========================================================
// Cover class 1/2/3 from the cover map
function classImage(cov) {
  return ee.Image(1)
    .where(cov.gte(COVER_BREAKS[0]), 2)
    .where(cov.gte(COVER_BREAKS[1]), 3)
    .updateMask(cov.mask())
    .rename('cover_class')
    .toInt();
}

// agc_g (g C per pixel), area_m2, n_px for one year, masked by the extent
function sumStack(y) {
  var agc  = ee.Image(MAP_PREFIX + y).select('AGC_pred');
  var area = ee.Image.pixelArea();
  return agc.multiply(area).rename('agc_g')
    .addBands(area.rename('area_m2'))
    .addBands(ee.Image(1).rename('n_px'))
    .updateMask(agc.mask());
}


// ==========================================================
// 4. YEAR x ECOREGION x COVER CLASS (one grouped reduction per year)
// key = eco_id x 10 + cover_class
// ==========================================================
function nationalRows(y) {
  var cls = classImage(ee.Image(MAP_PREFIX + y).select('tSPC_pred'));
  var key = ecoImg.multiply(10).add(cls).rename('key');

  var res = sumStack(y).addBands(key).reduceRegion(gridArgs({
    reducer:   ee.Reducer.sum().repeat(3).group({groupField: 3, groupName: 'key'}),
    geometry:  AOI,
    maxPixels: 1e13,
    tileScale: TILE_SCALE
  }));

  return ee.FeatureCollection(ee.List(res.get('groups')).map(function (g) {
    g = ee.Dictionary(g);
    var k = ee.Number(g.get('key')).round();
    var s = ee.List(g.get('sum'));
    return ee.Feature(null, {
      year:        y,
      eco_id:      k.divide(10).floor(),
      cover_class: k.mod(10),
      agc_g:       s.get(0),
      area_m2:     s.get(1),
      n_px:        s.get(2)
    });
  }));
}

if (RUN.national) {
  YEARS.forEach(function (y) {
    Export.table.toDrive({
      collection:  nationalRows(y),
      description: 'P4_gee_totals_' + y,
      folder:      DRIVE_FOLDER,
      fileFormat:  'CSV',
      selectors:   ['year', 'eco_id', 'cover_class', 'area_m2', 'agc_g', 'n_px']
    });
  });
  print('Export queued: ' + YEARS.length + ' x P4_gee_totals_YYYY');
}


// ==========================================================
// 5. ROI TOTALS (national extent inside each ROI polygon)
// ==========================================================
if (RUN.roi) {
  var roiRows = ee.FeatureCollection(YEARS.map(function (y) {
    return sumStack(y).reduceRegions(gridArgs({
      collection: roi,
      reducer:    ee.Reducer.sum(),
      tileScale:  TILE_SCALE
    })).map(function (f) { return f.set('year', y); });
  })).flatten();

  Export.table.toDrive({
    collection:  roiRows,
    description: 'P4_gee_roi_totals',
    folder:      DRIVE_FOLDER,
    fileFormat:  'CSV',
    selectors:   [ROI_ID_FIELD, 'year', 'area_m2', 'agc_g', 'n_px']
  });
  print('Export queued: P4_gee_roi_totals');
}

if (RUN.roiPolygons) {
  Export.table.toDrive({
    collection:  roi,
    description: 'P4_roi_polygons',
    folder:      DRIVE_FOLDER,
    fileFormat:  'GeoJSON'
  });
  print('Export queued: P4_roi_polygons (GeoJSON) -> set CFG$roi_path in R');
}


// ==========================================================
// 6. MAP CHECK
// ==========================================================
var m23 = ee.Image(MAP_PREFIX + 2023);
Map.centerObject(AOI, 5);
Map.addLayer(ecoImg.randomVisualizer(), {}, 'MEOW eco_id', false);
Map.addLayer(m23.select('AGC_pred'),
  {min: 0, max: 80, palette: ['ffffcc', 'a1dab4', '41b6c4', '2c7fb8', '253494']}, 'AGC 2023', false);
Map.addLayer(classImage(m23.select('tSPC_pred')),
  {min: 1, max: 3, palette: ['fee391', '78c679', '006837']}, 'Cover class 2023', false);
Map.addLayer(roi.style({color: 'ff0000', fillColor: '00000000'}), {}, 'ROIs', false);
