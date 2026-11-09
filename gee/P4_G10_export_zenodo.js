/************************************************************
 * P4_G10_export_zenodo.js
 * Exports the final seagrass products as GeoTIFFs to Google Drive, for the
 * open data deposit (Zenodo). Added 9 Nov 2026. Reads the P4_G01 and P4_G06
 * assets only; no asset is created or overwritten.
 *
 * FILES (Drive folder DRIVE_FOLDER, Cloud Optimized GeoTIFF)
 *   P4_seagrass_extent.tif
 *       1 band 'extent', UInt8: 1 = core extent (P4_G01),
 *       2 = add-on extent (P4_G06, shallow pixels without ACA geomorphic data).
 *       NoData = 0.
 *   P4_seagrass_AGC_cover_YYYY.tif   (one file per year, 2017-2024)
 *       2 bands, UInt16, scaled by 10 (divide by 10 to read the value):
 *         AGC_x10    above-ground carbon (gC m-2)
 *         cover_x10  seagrass percent cover (%)
 *       NoData = 65535. Pixels outside the final extent are NoData.
 *   Grid: the native grid of the maps (EPSG:4326, about 10 m), so every
 *   exported pixel is one mapped pixel. Areas must be computed per pixel
 *   (cell area changes with latitude), e.g. ee.Image.pixelArea() or
 *   cellSize() in R terra.
 *   GEE splits a large file into tiles (suffix -0000000000-0000000000);
 *   tiles without seagrass are skipped. Put all tiles of one product in the
 *   deposit; a GIS reads them as one layer (e.g. a VRT).
 *
 * NEEDS: P4_extent_mask, P4_AGC_COV_YYYY (P4_G01), P4_extent_addon_mask,
 *   P4_AGC_COV_ADDON_YYYY (P4_G06), the AOI asset.
 * RUN: press Run, then start the tasks in the Tasks tab (9 tasks).
 ************************************************************/

// ==========================================================
// 1. CONFIGURATION -- edit here only (paths as in P4_G08)
// ==========================================================
var ASSET_ROOT   = 'projects/YOUR-GEE-PROJECT/assets';   // your Earth Engine asset root
print('P4 script set 2026-11-09 (Zenodo GeoTIFF export)');

var OUT_DIR       = ASSET_ROOT + '/chapter_3/output/P4v2';
var CORE_PREFIX   = OUT_DIR + '/P4_AGC_COV_';
var ADDON_PREFIX  = OUT_DIR + '/P4_AGC_COV_ADDON_';
var CORE_EXTENT   = OUT_DIR + '/P4_extent_mask';
var ADDON_EXTENT  = OUT_DIR + '/P4_extent_addon_mask';
var AOI_ASSET     = ASSET_ROOT + '/packard_fieldData/ecoregion_indo_diss';
var DRIVE_FOLDER  = 'P4v2_ZENODO';          // new folder, separate from P4v2_GEE

var YEARS         = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];
var SCALE_FACTOR  = 10;                     // stored value = value x 10 (0.1 precision)
var NODATA_MAP    = 65535;
var FILE_DIM      = 65536;                  // pixels per file side (multiple of 256)

var RUN = {
  extent: true,     // 1 task : P4_seagrass_extent
  maps:   true      // 8 tasks: P4_seagrass_AGC_cover_YYYY
};


// ==========================================================
// 2. GRID, REGION AND FINAL EXTENT
// ==========================================================
var AOI    = ee.FeatureCollection(AOI_ASSET).geometry();
var REGION = AOI.bounds();
var GRID   = ee.Image(CORE_PREFIX + YEARS[0]).projection().getInfo();
var GRID_CRS = GRID.crs || GRID.wkt;
var NATIVE   = ee.Projection(GRID_CRS, GRID.transform);
print('Export grid (native map grid)', GRID_CRS, GRID.transform);

// Final extent = core or add-on (the two masks do not overlap; P4_G06)
var coreExt  = ee.Image(CORE_EXTENT).select(0).unmask(0, false).gt(0);
var addonExt = ee.Image(ADDON_EXTENT).select(0).unmask(0, false).gt(0);
var extentCode = coreExt.multiply(1).add(addonExt.multiply(2))
  .selfMask()
  .rename('extent')
  .toUint8()
  .setDefaultProjection(NATIVE);

// Core and add-on maps of one year as one image, on the native grid
function mapsOfYear(y) {
  return ee.ImageCollection([
      ee.Image(CORE_PREFIX + y).select(['AGC_pred', 'tSPC_pred']),
      ee.Image(ADDON_PREFIX + y).select(['AGC_pred', 'tSPC_pred'])])
    .mosaic()
    .setDefaultProjection(NATIVE)
    .updateMask(extentCode.mask())
    .multiply(SCALE_FACTOR).round()
    .rename(['AGC_x10', 'cover_x10'])
    .toUint16();
}

Map.centerObject(REGION, 5);
Map.addLayer(extentCode, {min: 1, max: 2, palette: ['1a9850', 'fdae61']}, 'Final extent (1 core, 2 add-on)');
Map.addLayer(mapsOfYear(2023).select('AGC_x10'), {min: 150, max: 400, palette: ['440154', '21908C', 'FDE725']},
             'AGC 2023 x 10', false);


// ==========================================================
// 3. EXPORTS
// ==========================================================
if (RUN.extent) {
  Export.image.toDrive({
    image:          extentCode,
    description:    'P4_seagrass_extent',
    fileNamePrefix: 'P4_seagrass_extent',
    folder:         DRIVE_FOLDER,
    region:         REGION,
    crs:            GRID_CRS,
    crsTransform:   GRID.transform,
    maxPixels:      1e13,
    fileDimensions: FILE_DIM,
    skipEmptyTiles: true,
    fileFormat:     'GeoTIFF',
    formatOptions:  {cloudOptimized: true, noData: 0}
  });
  print('Export queued: P4_seagrass_extent');
}

if (RUN.maps) {
  YEARS.forEach(function (y) {
    Export.image.toDrive({
      image:          mapsOfYear(y),
      description:    'P4_seagrass_AGC_cover_' + y,
      fileNamePrefix: 'P4_seagrass_AGC_cover_' + y,
      folder:         DRIVE_FOLDER,
      region:         REGION,
      crs:            GRID_CRS,
      crsTransform:   GRID.transform,
      maxPixels:      1e13,
      fileDimensions: FILE_DIM,
      skipEmptyTiles: true,
      fileFormat:     'GeoTIFF',
      formatOptions:  {cloudOptimized: true, noData: NODATA_MAP}
    });
  });
  print('Export queued: ' + YEARS.length + ' annual AGC and cover files');
}
