/************************************************************
 * P4_G04_temporal.js
 * Paper 4 (revised design): what changes between years.
 *   A. Per-pixel variation of mapped AGC, 2017-2024 (mean, SD, CV %):
 *      asset for the CV map, and area per CV class x ecoregion.
 *   B. GSE inter-annual similarity. GSE vectors are unit length (checked in
 *      P4_01: norm 0.993-1.007), so the dot product of two years is their
 *      cosine similarity (1 = same embedding). Consecutive-year pairs
 *      2017-2018 ... 2023-2024 inside the extent:
 *        B1 percentiles of the mean and minimum similarity per ecoregion
 *        B2 stratified random pixels (N_PER_ECO per ecoregion; drawn from the
 *           AGC maps, similarity computed at those points only) with every
 *           pair, AGC per year and AGC CV, for R
 *   C. Similarity at the AGC sample points between the survey-year GSE
 *      (the training values) and the GSE of every map year 2017-2024.
 *      This is the predictor change a Tier 1 factor ignores when a survey
 *      from one year stands for another year.
 *
 * EXPORTS
 *   Asset : P4_AGC_CV_2017_2024 (AGC_mean, AGC_sd, AGC_cv)
 *   Drive : P4_gee_cv_summary.csv        eco_id, cv_class, area_m2, cv_x_area
 *           P4_gee_gse_similarity_eco.csv eco_id (-1 = all), band, p5..p95, mean, count
 *           P4_gee_gse_sim_sample.csv    one row per sampled pixel
 *           P4_gee_gse_sim_points.csv    rid, gee_id, year, eco_id, AGC_pred, sim_table,
 *                                        sim_2017..sim_2024
 *
 * NEEDS: P4_extent_mask and the 8 map assets (P4_G01), P4_meow_indo,
 *        P4_gee_training
 ************************************************************/

// ==========================================================
// 1. CONFIGURATION -- edit here only
// ==========================================================
var ASSET_ROOT   = 'projects/YOUR-GEE-PROJECT/assets';   // your Earth Engine asset root
print('P4 script set 2026-10-31');

// All new assets (inputs and outputs) live in OUT_DIR with a P4_ prefix; the
// existing chapter_2 / chapter_3 assets are only read, never written.
var OUT_DIR      = ASSET_ROOT + '/chapter_3/output/P4v2';
var TRAIN_TABLE  = OUT_DIR + '/P4_gee_training';               // from P4_06
var MAP_PREFIX   = OUT_DIR + '/P4_AGC_COV_';
var EXTENT_ASSET = OUT_DIR + '/P4_extent_mask';
var CV_ASSET     = OUT_DIR + '/P4_AGC_CV_2017_2024';
var MEOW_ASSET   = OUT_DIR + '/P4_meow_indo';
var AOI_ASSET    = ASSET_ROOT + '/packard_fieldData/ecoregion_indo_diss';
var EMB_COLL     = 'GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL';
var DRIVE_FOLDER = 'P4v2_GEE';

var MAP_YEARS    = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];
var SCALE        = 10;
var SEED         = 42;
var NA_VALUE     = -9999;

// CV classes (%): 0 <10, 1 10-20, 2 20-30, 3 30-50, 4 50-100, 5 >=100
var CV_BREAKS    = [10, 20, 30, 50, 100];
var N_PER_ECO    = 3000;       // B2: points per ecoregion (eco_id 0 included), ~45k in total
var SAMPLE_SCALE = 100;        // B2: the sample locations are drawn on a 100 m grid (the scan
                               //     of the whole AOI at 10 m took >12 h); values are then
                               //     read at the 10 m pixel under each point

var REDUCE_SCALE = null;       // null = native 10 m grid of the extent asset
var TILE_SCALE   = 8;

var RUN = {
  cvAsset:   true,   // A  asset
  cvSummary: true,   // A  table
  simEco:    false,  // B1 exact percentiles over every extent pixel: optional and
                     //    very slow (like P4_gee_distance_pct); P4_07 takes the
                     //    percentiles from the B2 sample instead
  simSample: true,   // B2
  simPoints: true    // C
};


// ==========================================================
// 2. INPUTS, GRID, HELPERS
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

var AOI    = ee.FeatureCollection(AOI_ASSET).geometry();

// Values of a reduceRegion(first) result in key order; null -> NA_VALUE
// (a point outside every GSE tile returns null even after unmask)
function safeValues(dict, keys) {
  return ee.List(ee.Dictionary(dict).values(keys)).map(function (x) {
    return ee.Algorithms.If(ee.Algorithms.IsEqual(x, null), NA_VALUE, x);
  });
}
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

var area = ee.Image.pixelArea();

// Point geometry rebuilt from xcoord/ycoord (EPSG:4326), so the point
// locations do not depend on the X/Y columns chosen when the CSV was uploaded
function withPoint(f) {
  return f.setGeometry(ee.Geometry.Point([ee.Number(f.get('xcoord')), ee.Number(f.get('ycoord'))]));
}


// ==========================================================
// A. PER-PIXEL AGC VARIATION 2017-2024
// ==========================================================
var agcCol = ee.ImageCollection(MAP_YEARS.map(function (y) {
  return ee.Image(MAP_PREFIX + y).select('AGC_pred');
}));
var agcYears = ee.Image.cat(MAP_YEARS.map(function (y) {
  return ee.Image(MAP_PREFIX + y).select('AGC_pred').rename('AGC_' + y);
}));

var agcMean = agcCol.mean().rename('AGC_mean');
var agcSd   = agcCol.reduce(ee.Reducer.sampleStdDev()).rename('AGC_sd');
var agcCv   = agcSd.divide(agcMean).multiply(100).rename('AGC_cv');
var cvImg   = agcMean.addBands(agcSd).addBands(agcCv).toFloat().updateMask(extent);

if (RUN.cvAsset) {
  Export.image.toAsset({
    image:            cvImg.set({years: '2017-2024', sd: 'sample SD (n - 1)'}),
    description:      'P4_AGC_CV_2017_2024',
    assetId:          CV_ASSET,
    region:           AOI,
    crs:              GRID_CRS,
    crsTransform:     GRID.transform,
    maxPixels:        1e13,
    pyramidingPolicy: {'.default': 'mean'}
  });
  print('Export queued: P4_AGC_CV_2017_2024 (asset)');
}

if (RUN.cvSummary) {
  var cvCls = ee.Image(CV_BREAKS.length);
  for (var k = CV_BREAKS.length - 1; k >= 0; k--) {
    cvCls = cvCls.where(agcCv.lt(CV_BREAKS[k]), k);
  }
  cvCls = cvCls.updateMask(agcCv.mask()).toInt();
  var cvKey = ecoImg.multiply(10).add(cvCls).rename('key');

  var cvRes = area.rename('area_m2')
    .addBands(agcCv.multiply(area).rename('cv_x_area'))
    .addBands(cvKey)
    .updateMask(extent)
    .reduceRegion(gridArgs({
      reducer:   ee.Reducer.sum().repeat(2).group({groupField: 2, groupName: 'key'}),
      geometry:  AOI,
      maxPixels: 1e13,
      tileScale: TILE_SCALE
    }));

  var cvRows = ee.FeatureCollection(ee.List(cvRes.get('groups')).map(function (g) {
    g = ee.Dictionary(g);
    var kk = ee.Number(g.get('key')).round();
    var s  = ee.List(g.get('sum'));
    return ee.Feature(null, {eco_id: kk.divide(10).floor(), cv_class: kk.mod(10),
                             area_m2: s.get(0), cv_x_area: s.get(1)});
  }));

  Export.table.toDrive({
    collection:  cvRows,
    description: 'P4_gee_cv_summary',
    folder:      DRIVE_FOLDER,
    fileFormat:  'CSV',
    selectors:   ['eco_id', 'cv_class', 'area_m2', 'cv_x_area']
  });
  print('Export queued: P4_gee_cv_summary');
}


// ==========================================================
// B. GSE SIMILARITY BETWEEN CONSECUTIVE YEARS (inside the extent)
// ==========================================================
var simBands = [];
var SIM_NAMES = [];
for (var p = 0; p < MAP_YEARS.length - 1; p++) {
  var y0 = MAP_YEARS[p], y1 = MAP_YEARS[p + 1];
  var nm = 'sim_' + y0 + '_' + y1;
  SIM_NAMES.push(nm);
  simBands.push(gseImage(y0).multiply(gseImage(y1)).reduce(ee.Reducer.sum()).rename(nm));
}
var simImg  = ee.Image.cat(simBands);
var simMean = simImg.reduce(ee.Reducer.mean()).rename('sim_mean');
var simMin  = simImg.reduce(ee.Reducer.min()).rename('sim_min');
var simAll  = simImg.addBands(simMean).addBands(simMin).updateMask(extent);

if (RUN.simEco) {
  var pctReducer = ee.Reducer.percentile([5, 25, 50, 75, 95])
    .combine({reducer2: ee.Reducer.mean(), sharedInputs: true})
    .combine({reducer2: ee.Reducer.count(), sharedInputs: true});

  var simPct = function (band, groupImg) {
    var r = simAll.select(band).addBands(groupImg).reduceRegion(gridArgs({
      reducer:   pctReducer.group({groupField: 1, groupName: 'eco_id'}),
      geometry:  AOI,
      maxPixels: 1e13,
      tileScale: TILE_SCALE
    }));
    return ee.FeatureCollection(ee.List(r.get('groups')).map(function (g) {
      return ee.Feature(null, ee.Dictionary(g)).set('band', band);
    }));
  };
  var allImg = ee.Image(-1).toInt().rename('eco_id');
  var simEcoRows = simPct('sim_mean', ecoImg)
    .merge(simPct('sim_mean', allImg))
    .merge(simPct('sim_min', ecoImg))
    .merge(simPct('sim_min', allImg));

  Export.table.toDrive({
    collection:  simEcoRows,
    description: 'P4_gee_gse_similarity_eco',
    folder:      DRIVE_FOLDER,
    fileFormat:  'CSV',
    selectors:   ['eco_id', 'band', 'p5', 'p25', 'p50', 'p75', 'p95', 'mean', 'count']
  });
  print('Export queued: P4_gee_gse_similarity_eco');
}

if (RUN.simSample) {
  // Step 1: draw the locations on a coarse grid from a one-band image
  // (ecoregion inside the extent), so the scan of the AOI stays cheap.
  var locs = ecoImg.updateMask(extent).stratifiedSample({
    numPoints:  N_PER_ECO,
    classBand:  'eco_id',
    region:     AOI,
    projection: ee.Projection('EPSG:4326').atScale(SAMPLE_SCALE),
    seed:       SEED,
    dropNulls:  true,
    tileScale:  16,
    geometries: true
  });

  // Step 2: read every value at the 10 m pixel under each point (extent grid).
  // A 100 m cell can be partly outside the extent: its centre pixel then has
  // no values and P4_07 drops that row.
  var valueImg = simAll
    .addBands(agcMean).addBands(agcSd).addBands(agcCv)
    .addBands(agcYears);

  var sampSim = locs.map(function (f) {
    var v = valueImg.reduceRegion({
      reducer: ee.Reducer.first(), geometry: f.geometry(),
      crs: GRID_CRS, crsTransform: GRID.transform, maxPixels: 1e13
    });
    var ll = f.geometry().coordinates();
    return f.set(v).set({longitude: ll.get(0), latitude: ll.get(1)});
  });

  var sampCols = ['eco_id', 'longitude', 'latitude'].concat(SIM_NAMES)
    .concat(['sim_mean', 'sim_min', 'AGC_mean', 'AGC_sd', 'AGC_cv'])
    .concat(MAP_YEARS.map(function (y) { return 'AGC_' + y; }));

  Export.table.toDrive({
    collection:  sampSim,
    description: 'P4_gee_gse_sim_sample',
    folder:      DRIVE_FOLDER,
    fileFormat:  'CSV',
    selectors:   sampCols
  });
  print('Export queued: P4_gee_gse_sim_sample');
}


// ==========================================================
// C. SIMILARITY AT THE AGC SAMPLE POINTS: survey year vs every map year
// Both vectors come from the GSE images with the same extraction
// (reduceRegion per point, scale 10, as 01_trainingData_prep_indo.js), so
// sim at lag 0 is 1 by construction. About a fifth of the table vectors
// sit on a neighbouring pixel (P4_G01 check); sim_table records that match
// (table vector x image vector, survey year). No GSE -> NA_VALUE.
// ==========================================================
if (RUN.simPoints) {
  var KEEP = ['rid', 'gee_id', 'year', 'eco_id', 'AGC_pred'];
  var simCols = MAP_YEARS.map(function (y) { return 'sim_' + y; });
  var imgs = {};
  MAP_YEARS.forEach(function (y) { imgs[y] = gseImage(y).unmask(NA_VALUE, false); });

  var vecAt = function (img, geom) {
    return safeValues(img.reduceRegion({
      reducer: ee.Reducer.first(), geometry: geom, scale: SCALE, maxPixels: 1e13
    }), GSE_BANDS);
  };
  var dotOf = function (a, b) {
    return ee.Array(a).multiply(ee.Array(b)).reduce(ee.Reducer.sum(), [0]).get([0]);
  };

  var pts = ee.FeatureCollection(TRAIN_TABLE).map(withPoint)
    .filter(ee.Filter.eq('is_agc', 1))
    .map(function (f) {
      var g   = f.geometry();
      var a   = vecAt(gseImage(f.get('year')).unmask(NA_VALUE, false), g);   // survey year, image
      var aNA = ee.Number(a.get(0)).eq(NA_VALUE);
      var props = {
        sim_table: ee.Algorithms.If(aNA, NA_VALUE,
                                    dotOf(f.toDictionary(GSE_BANDS).values(GSE_BANDS), a))
      };
      MAP_YEARS.forEach(function (y) {
        var b = vecAt(imgs[y], g);
        props['sim_' + y] = ee.Algorithms.If(aNA.or(ee.Number(b.get(0)).eq(NA_VALUE)),
                                             NA_VALUE, dotOf(a, b));
      });
      return ee.Feature(null, f.toDictionary(KEEP)).set(props);
    });

  Export.table.toDrive({
    collection:  pts,
    description: 'P4_gee_gse_sim_points',
    folder:      DRIVE_FOLDER,
    fileFormat:  'CSV',
    selectors:   KEEP.concat(['sim_table']).concat(simCols)
  });
  print('Export queued: P4_gee_gse_sim_points');
}


// ==========================================================
// MAP CHECK
// ==========================================================
Map.centerObject(AOI, 5);
Map.addLayer(cvImg.select('AGC_cv'), {min: 0, max: 50, palette: ['2166ac', 'f7f7f7', 'b2182b']},
             'AGC CV 2017-2024 (%)', false);
Map.addLayer(simMean.updateMask(extent), {min: 0.6, max: 1, palette: ['d7191c', 'ffffbf', '1a9641']},
             'GSE similarity, mean of consecutive years', false);
