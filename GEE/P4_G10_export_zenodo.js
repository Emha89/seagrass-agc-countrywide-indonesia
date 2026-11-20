/************************************************************
 * P4_G10_export_zenodo.js
 * Exports the final seagrass products as GeoTIFFs to Google Drive, for the
 * open data deposit (Zenodo). Added 9 Nov 2026. Reads the P4_G01 and P4_G06
 * assets only; no asset is created or overwritten.
 * v3 (20 Nov 2026): AGC only. The cover band is no longer exported, because
 * the cover model is not part of design v3. New: the mean AGC map of
 * 2017-2024 (Fig. 3 of the paper). New Drive folder P4v3_ZENODO.
 * Only map products are exported: no field sample, sample location or
 * distance-to-sample layer leaves Earth Engine.
 *
 * FILES (Drive folder DRIVE_FOLDER, Cloud Optimized GeoTIFF)
 *   P4_seagrass_extent.tif
 *       1 band 'extent', UInt8: 1 = core extent (P4_G01),
 *       2 = add-on extent (P4_G06, shallow pixels without ACA geomorphic data).
 *       NoData = 0.
 *   P4_seagrass_AGC_YYYY.tif   (one file per year, 2017-2024)
 *       1 band 'AGC_x10', UInt16: above-ground carbon (gC m-2) x 10
 *       (divide by 10 to read the value).
 *       NoData = 65535. Pixels outside the final extent are NoData.
 *   P4_seagrass_AGC_mean_2017_2024.tif
 *       1 band 'AGC_x10', UInt16: mean of the eight annual AGC maps x 10.
 *       NoData = 65535.
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
 * RUN: press Run, then start the tasks in the Tasks tab (10 tasks).
 * THEN: download the Drive folder into data/zenodo_maps of the R project and
 *   run R/P4_11_zenodo_package.R, which zips the tiles of each product.
 ************************************************************/

// ==========================================================
// 1. CONFIGURATION -- edit here only (paths as in P4_G08)
// ==========================================================
var ASSET_ROOT   = 'projects/YOUR-GEE-PROJECT/assets';   // your Earth Engine asset root
print('P4 script set 2026-11-20 (v3, Zenodo GeoTIFF export, AGC only)');

var OUT_DIR       = ASSET_ROOT + '/chapter_3/output/P4v2';
var CORE_PREFIX   = OUT_DIR + '/P4_AGC_COV_';
var ADDON_PREFIX  = OUT_DIR + '/P4_AGC_COV_ADDON_';
var CORE_EXTENT   = OUT_DIR + '/P4_extent_mask';
var ADDON_EXTENT  = OUT_DIR + '/P4_extent_addon_mask';
var AOI_ASSET     = ASSET_ROOT + '/packard_fieldData/ecoregion_indo_diss';
var DRIVE_FOLDER  = 'P4v3_ZENODO';          // v3 folder (P4v2_ZENODO held the cover band)

var YEARS         = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];
var SCALE_FACTOR  = 10;                     // stored value = value x 10 (0.1 precision)
var NODATA_MAP    = 65535;
var FILE_DIM      = 65536;                  // pixels per file side (multiple of 256)

var RUN = {
  extent: true,     // 1 task : P4_seagrass_extent
  maps:   true,     // 8 tasks: P4_seagrass_AGC_YYYY
  mean:   true      // 1 task : P4_seagrass_AGC_mean_2017_2024
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

// Core and add-on AGC of one year as one image (gC m-2), on the native grid
function agcOfYear(y) {
  return ee.ImageCollection([
      ee.Image(CORE_PREFIX + y).select('AGC_pred'),
      ee.Image(ADDON_PREFIX + y).select('AGC_pred')])
    .mosaic()
    .setDefaultProjection(NATIVE)
    .updateMask(extentCode.mask());
}

// Stored values: AGC x 10 as UInt16 (0.1 gC m-2 precision)
function toStored(img) {
  return img.multiply(SCALE_FACTOR).round().rename('AGC_x10').toUint16();
}
function mapsOfYear(y) { return toStored(agcOfYear(y)); }

// Mean of the eight annual maps (Fig. 3)
var agcMean = toStored(ee.ImageCollection(YEARS.map(agcOfYear)).mean()
  .setDefaultProjection(NATIVE));

Map.centerObject(REGION, 5);
Map.addLayer(extentCode, {min: 1, max: 2, palette: ['1a9850', 'fdae61']}, 'Final extent (1 core, 2 add-on)');
Map.addLayer(mapsOfYear(2023), {min: 150, max: 400, palette: ['440154', '21908C', 'FDE725']},
             'AGC 2023 x 10', false);
Map.addLayer(agcMean, {min: 150, max: 400, palette: ['440154', '21908C', 'FDE725']},
             'AGC mean 2017-2024 x 10', false);


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
      description:    'P4_seagrass_AGC_' + y,
      fileNamePrefix: 'P4_seagrass_AGC_' + y,
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
  print('Export queued: ' + YEARS.length + ' annual AGC files');
}

if (RUN.mean) {
  Export.image.toDrive({
    image:          agcMean,
    description:    'P4_seagrass_AGC_mean_2017_2024',
    fileNamePrefix: 'P4_seagrass_AGC_mean_2017_2024',
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
  print('Export queued: P4_seagrass_AGC_mean_2017_2024');
}
