/************************************************************
 * P4_G08b_zoom_windows.js
 * Paper 4 (revised design), added 6 Nov 2026: 10 m zoom windows for the map
 * figures (part B of P4_G08). The national panels use the ~500 m grid of
 * P4_G08_figure_grid.js; these windows show the 10 m pixels at four sites.
 *
 * WINDOWS (0.06 x 0.06 deg, ~6.7 km; centres chosen in R from the 500 m grid
 * and the AGC samples)
 *   Rote2    centred on the main Rote sample cluster (1,534 AGC samples,
 *            2021/2022; Lesser Sunda, clear water)              [7 Nov 2026]
 *   Derawan2 centred on the Derawan sample cluster (412 samples, 2017/2021;
 *            small-leaf meadows, low field AGC)                  [7 Nov 2026]
 *   Seribu   largest extent among windows with >= 20 AGC samples (749
 *            samples around Pulau Pari, 2021/2022)               [6 Nov 2026]
 *   Selayar  largest extent on the Selayar coast; every pixel > 50 km from
 *            the nearest AGC sample (Banda Sea, no field data)   [6 Nov 2026]
 *   The first Rote and Derawan windows (6 Nov) had the samples at the window
 *   edge; Rote2 and Derawan2 replace them (new names, nothing overwritten).
 *
 * WHAT THIS SCRIPT DOES
 *   For each window, on the native 10 m grid of the map assets:
 *   - layers (masked to the final extent = core + add-on):
 *       AGC_mean     mean AGC 2017-2024 (gC m-2), from the CV assets (P4_G04, G07)
 *       tSPC_mean    mean predicted cover 2017-2024 (%), from the 8 annual maps
 *       cover_class  1/2/3 = <30 / 30-59.9 / >=60 % of tSPC_mean (Tier 2 strata)
 *       dist_class   0-10, distance to the nearest AGC sample (classes as P4_G03)
 *       extent_part  1 = core extent (P4_G01), 2 = add-on (P4_G06)
 *   - Sentinel-2 surface reflectance true colour (B4, B3, B2), median of
 *     S2_YEARS after Cloud Score+ masking, not masked (background).
 *   Nothing is written to assets; only Drive GeoTIFFs with new names.
 *
 * EXPORTS (Google Drive, folder DRIVE_FOLDER)
 *   P4_zoom_<site>_layers.tif      5 bands above, Float32, nodata outside the extent
 *   P4_zoom_<site>_s2rgb.tif       B4, B3, B2, UInt16 (reflectance x 10,000)
 *   P4_zoom_<site>_agc_annual.tif  AGC_2017 ... AGC_2024 (gC m-2), Float32,
 *                                  nodata outside the extent  [added 7 Nov 2026]
 *   P4_zoom_<site>_maskdiag.tif    why a pixel is in or out of the extent, not
 *                                  masked [added 7 Nov 2026]:
 *                                  mask_code = 1 x (MaxProb_0_7 mask absent)
 *                                            + 2 x (ACA depth not in 0-500, incl. no data)
 *                                            + 4 x (ACA geomorphic class 16)
 *                                  (0 = in the final extent), extent_final (0/1),
 *                                  depth_raw (-9999 = no data), geomorphic (-9999 = no data)
 *   layers and s2rgb are exported only for windows with redo: true (the
 *   Seribu and Selayar files of 6 Nov are kept); agc_annual for all windows.
 *
 * WHY THE MASK DIAGNOSTIC (7 Nov 2026): around Rote, 677 of 1,538 field
 * samples (mean cover 76 %) are outside the extent only because of ACA
 * geomorphic class 16 (P4_gee_samples_extent). The diagnostic shows which
 * mask removes which pixels. The first Rote window (6 Nov) is included
 * because dark, seagrass-like areas were seen in its northern part.
 *
 * WHY THE ANNUAL BANDS (7 Nov 2026): at the sample pixels, the map of the
 * survey year reproduces the field values (in-sample: Derawan 3.5 vs field
 * 2.8, Seribu 44.8 vs 44.9 gC m-2), while the 2017-2024 mean is pulled toward
 * the national mean (Derawan 14.2, Seribu 36.6). The annual bands show this.
 *
 * NEEDS: P4_extent_mask, P4_AGC_COV_YYYY, P4_AGC_CV_2017_2024 (P4_G01, G04),
 *        P4_extent_addon_mask, P4_AGC_COV_ADDON_YYYY, P4_AGC_CV_ADDON_2017_2024
 *        (P4_G06, G07), P4_gee_training (P4_06)
 ************************************************************/

// ==========================================================
// 1. CONFIGURATION -- edit here only (paths and classes as in P4_G02-G08)
// ==========================================================
var ASSET_ROOT   = 'projects/YOUR-GEE-PROJECT/assets';   // your Earth Engine asset root
print('P4 script set 2026-11-07b (zoom windows v2 + mask diagnostic)');

var OUT_DIR       = ASSET_ROOT + '/chapter_3/output/P4v2';
var TRAIN_TABLE   = OUT_DIR + '/P4_gee_training';
var CORE_PREFIX   = OUT_DIR + '/P4_AGC_COV_';
var ADDON_PREFIX  = OUT_DIR + '/P4_AGC_COV_ADDON_';
var CORE_EXTENT   = OUT_DIR + '/P4_extent_mask';
var ADDON_EXTENT  = OUT_DIR + '/P4_extent_addon_mask';
var CORE_CV       = OUT_DIR + '/P4_AGC_CV_2017_2024';
var ADDON_CV      = OUT_DIR + '/P4_AGC_CV_ADDON_2017_2024';
// Mask inputs, as P4_G01 / P4_G05 / P4_G06 (for the mask diagnostic only)
var MAXPROB_MASK  = ASSET_ROOT + '/chapter_3/output/02_PA/indo_Seagrass_MaxProb_0_7_mask';
var DEPTH_ASSET   = ASSET_ROOT + '/chapter_2/InputArea/raster_INA_ACA_Bathymetry';
var ACA_ASSET     = ASSET_ROOT + '/chapter_2/InputArea/acav2_indoarch_stack';
var DRIVE_FOLDER  = 'P4v2_GEE';

var YEARS         = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];
var COVER_BREAKS  = [30, 60];                                     // as P4_G02
var SEARCH_RADIUS  = 50000;                                       // as P4_G03
var DIST_MAX_ERROR = 5;
var DIST_LIMITS = [10, 50, 100, 250, 500, 1000, 2000, 5000, 20000, SEARCH_RADIUS];

// Window centres (lon, lat) and size. expectKm2 = extent area in the window
// from the 500 m grid (cells whose centre falls inside), for the console check.
var WINDOW_DEG = 0.06;
// redo: export layers and s2rgb again (false = the 6 Nov files are kept)
// annual: export the 8 annual AGC bands; every window gets the mask diagnostic
var WINDOWS = [
  {site: 'Rote2',    lon: 122.820, lat: -10.899, expectKm2: 2.53, redo: true,  annual: true},
  {site: 'Derawan2', lon: 118.246, lat:   2.283, expectKm2: 0.57, redo: true,  annual: true},
  {site: 'Seribu',   lon: 106.610, lat:  -5.875, expectKm2: 5.67, redo: false, annual: true},
  {site: 'Selayar',  lon: 120.455, lat:  -6.300, expectKm2: 9.20, redo: false, annual: true},
  {site: 'Rote',     lon: 122.800, lat: -10.860, expectKm2: 4.11, redo: false, annual: false}  // 6 Nov window
];

// Sentinel-2 background: median of these years after Cloud Score+ masking
var S2_START  = '2020-01-01';
var S2_END    = '2025-01-01';
var CS_THRESH = 0.60;          // keep pixels with cs_cdf >= 0.60

var RUN = {
  layers: true,   // P4_zoom_<site>_layers      (windows with redo: true)
  s2rgb:  true,   // P4_zoom_<site>_s2rgb       (windows with redo: true)
  annual: true,   // P4_zoom_<site>_agc_annual  (windows with annual: true)
  maskdiag: true, // P4_zoom_<site>_maskdiag    (all windows)
  check:  true    // console: extent area per window (small reductions, seconds)
};


// ==========================================================
// 2. INPUTS, GRID AND FINAL EXTENT (as P4_G08_figure_grid.js)
// ==========================================================
var GRID = ee.Image(CORE_PREFIX + YEARS[0]).projection().getInfo();
var GRID_CRS = GRID.crs || GRID.wkt;
var NATIVE   = ee.Projection(GRID_CRS, GRID.transform);

var coreExt  = ee.Image(CORE_EXTENT).select(0).unmask(0, false).gt(0);
var addonExt = ee.Image(ADDON_EXTENT).select(0).unmask(0, false).gt(0);
var extentAll = coreExt.or(addonExt).selfMask();

function combined(coreId, addonId, band) {
  return ee.ImageCollection([ee.Image(coreId).select(band), ee.Image(addonId).select(band)])
    .mosaic()
    .setDefaultProjection(NATIVE);
}

function windowGeom(w) {
  var h = WINDOW_DEG / 2;
  return ee.Geometry.Rectangle([w.lon - h, w.lat - h, w.lon + h, w.lat + h], null, false);
}


// ==========================================================
// 3. LAYERS
// ==========================================================
var agcMean = combined(CORE_CV, ADDON_CV, 'AGC_mean').rename('AGC_mean');

var tspcMean = ee.ImageCollection(YEARS.map(function (y) {
  return combined(CORE_PREFIX + y, ADDON_PREFIX + y, 'tSPC_pred');
})).mean().setDefaultProjection(NATIVE).rename('tSPC_mean');

var coverClass = ee.Image(1)
  .where(tspcMean.gte(COVER_BREAKS[0]), 2)
  .where(tspcMean.gte(COVER_BREAKS[1]), 3)
  .rename('cover_class');

var agcPts = ee.FeatureCollection(TRAIN_TABLE)
  .map(function (f) {
    return f.setGeometry(ee.Geometry.Point([ee.Number(f.get('xcoord')), ee.Number(f.get('ycoord'))]));
  })
  .filter(ee.Filter.eq('is_agc', 1));
print('AGC samples (R: 6616)', agcPts.size());

var dist = agcPts.distance(SEARCH_RADIUS, DIST_MAX_ERROR)
  .unmask(SEARCH_RADIUS * 2, false)          // beyond the search radius -> last class
  .setDefaultProjection(NATIVE);
var distClass = ee.Image(DIST_LIMITS.length);
for (var k = DIST_LIMITS.length - 1; k >= 0; k--) distClass = distClass.where(dist.lte(DIST_LIMITS[k]), k);
distClass = distClass.rename('dist_class');

var extentPart = ee.Image(0).where(coreExt, 1).where(addonExt, 2).rename('extent_part');

var agcAnnual = ee.Image(YEARS.map(function (y) {      // ee.Image(list) = one band per element
  return combined(CORE_PREFIX + y, ADDON_PREFIX + y, 'AGC_pred').rename('AGC_' + y);
})).toFloat().updateMask(extentAll);

// Mask diagnostic: the three masks of the extent, filled everywhere
var probFail  = ee.Image(MAXPROB_MASK).select('seagrass_presence').unmask(0, false).gt(0).not();
var depthRaw  = ee.Image(DEPTH_ASSET).select(0);
var depth0    = depthRaw.unmask(0, false);                       // no data -> 0 -> fails, as P4_G01
var depthFail = depth0.gt(0).and(depth0.lte(500)).not();
var geoRaw    = ee.Image(ACA_ASSET).select('geomorphic');
var class16   = geoRaw.unmask(0, false).eq(16);
var maskDiag = probFail.multiply(1).add(depthFail.multiply(2)).add(class16.multiply(4)).rename('mask_code')
  .addBands(extentAll.unmask(0, false).rename('extent_final'))
  .addBands(depthRaw.unmask(-9999, false).rename('depth_raw'))
  .addBands(geoRaw.unmask(-9999, false).rename('geomorphic'))
  .toFloat()
  .setDefaultProjection(NATIVE);

var layers = agcMean
  .addBands(tspcMean)
  .addBands(coverClass)
  .addBands(distClass)
  .addBands(extentPart)
  .toFloat()
  .updateMask(extentAll);


// ==========================================================
// 4. SENTINEL-2 BACKGROUND (Cloud Score+ masked median)
// ==========================================================
function s2Rgb(geom) {
  var cs = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
  return ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(geom)
    .filterDate(S2_START, S2_END)
    .linkCollection(cs, ['cs_cdf'])
    .map(function (img) { return img.updateMask(img.select('cs_cdf').gte(CS_THRESH)); })
    .select(['B4', 'B3', 'B2'])
    .median()
    .toUint16();
}


// ==========================================================
// 5. EXPORTS AND CHECKS
// ==========================================================
WINDOWS.forEach(function (w) {
  var geom = windowGeom(w);

  if (RUN.layers && w.redo) {
    Export.image.toDrive({
      image:        layers,
      description:  'P4_zoom_' + w.site + '_layers',
      folder:       DRIVE_FOLDER,
      region:       geom,
      crs:          GRID_CRS,
      crsTransform: GRID.transform,
      maxPixels:    1e9,
      fileFormat:   'GeoTIFF'
    });
  }
  if (RUN.s2rgb && w.redo) {
    Export.image.toDrive({
      image:        s2Rgb(geom),
      description:  'P4_zoom_' + w.site + '_s2rgb',
      folder:       DRIVE_FOLDER,
      region:       geom,
      crs:          GRID_CRS,
      crsTransform: GRID.transform,
      maxPixels:    1e9,
      fileFormat:   'GeoTIFF'
    });
  }
  if (RUN.maskdiag) {
    Export.image.toDrive({
      image:        maskDiag,
      description:  'P4_zoom_' + w.site + '_maskdiag',
      folder:       DRIVE_FOLDER,
      region:       geom,
      crs:          GRID_CRS,
      crsTransform: GRID.transform,
      maxPixels:    1e9,
      fileFormat:   'GeoTIFF'
    });
  }
  if (RUN.annual && w.annual) {
    Export.image.toDrive({
      image:        agcAnnual,
      description:  'P4_zoom_' + w.site + '_agc_annual',
      folder:       DRIVE_FOLDER,
      region:       geom,
      crs:          GRID_CRS,
      crsTransform: GRID.transform,
      maxPixels:    1e9,
      fileFormat:   'GeoTIFF'
    });
  }
  if (RUN.check) {
    var a = ee.Image.pixelArea().updateMask(extentAll).reduceRegion({
      reducer: ee.Reducer.sum(), geometry: geom, crs: GRID_CRS, crsTransform: GRID.transform,
      maxPixels: 1e9
    }).values().get(0);
    print(w.site + ': extent area in window (km2), expected about ' + w.expectKm2,
          ee.Number(a).divide(1e6));
  }
});
var nRedo = WINDOWS.filter(function (w) { return w.redo; }).length;
var nAnnual = WINDOWS.filter(function (w) { return w.annual; }).length;
print('Export tasks queued: ' + ((RUN.layers ? nRedo : 0) + (RUN.s2rgb ? nRedo : 0) +
      (RUN.annual ? nAnnual : 0) + (RUN.maskdiag ? WINDOWS.length : 0)));


// ==========================================================
// MAP CHECK (display only)
// ==========================================================
var first = windowGeom(WINDOWS[0]);
Map.centerObject(first, 13);
Map.addLayer(s2Rgb(first), {bands: ['B4', 'B3', 'B2'], min: 0, max: 1500}, 'Sentinel-2 (first window)');
Map.addLayer(layers.select('AGC_mean'), {min: 0, max: 60, palette: ['ffffcc', '78c679', '006837']}, 'AGC mean');
Map.addLayer(layers.select('dist_class'), {min: 0, max: 10, palette: ['2166ac', 'f7f7f7', 'b2182b']}, 'Distance class', false);
WINDOWS.forEach(function (w) { Map.addLayer(windowGeom(w), {color: 'red'}, 'Window ' + w.site, false); });
