/************************************************************
 * export_ROI_collection.js
 * Export named ROI geometries to a single asset (1 FeatureCollection).
 *
 * UPDATE WORKFLOW: add a new ROI to ROI_LIST, delete the old asset,
 * run the script, Tasks -> Run -> done.
 *
 * REQUIRES GEE IMPORTS: the 19 geometry variables referenced in
 * ROI_LIST below (rote01, pari, ayau01, padaido, komodo, kaledupa01,
 * karimun, baranglompo, kemujan, rote02, rote03, parang, nias,
 * teluk_bakau, berakit, tana_merah, derawan, lembongan, tapil) are not
 * defined in this file -- they come from the GEE Code Editor's
 * "Imports" panel (hand-drawn or imported boundary geometries specific
 * to each site). You'll need to recreate or import equivalent
 * geometries under those same variable names before running this
 * script.
 *
 * OUTPUT ASSET:
 *   .../AOI_AGCcompare/ROI_collection
 *
 * This produces the ROI_collection asset that
 * 09_AGC_rf_compare_indo.js reads via ROI_ASSET, filtering to one
 * ROI at a time by loc_id.
 ************************************************************/

var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var AOI_ASSET = ASSET_ROOT + '/chapter_3/output/AOI_AGCcompare/ROI_collection';

// ==========================================================
// ROI LIST -- geometries referenced directly from Imports
// ==========================================================
var ROI_LIST = [
  { id: 'rote',         label: 'Rote - Flores Alor',        geom: rote01 },
  { id: 'pari',         label: 'Pari - Jakarta Bay',         geom: pari },
  { id: 'ayau',         label: 'Ayau - Raja Ampat',          geom: ayau01 },
  { id: 'padaido',      label: 'Padaido - Biak',             geom: padaido },
  { id: 'komodo',       label: 'Komodo - NTT',               geom: komodo },
  { id: 'kaledupa',     label: 'Kaledupa - Wakatobi',        geom: kaledupa01 },
  { id: 'karimunjawa',  label: 'Karimunjawa - Jawa',         geom: karimun },
  { id: 'baranglompo',  label: 'Barang Lompo - Makassar',    geom: baranglompo },
  { id: 'kemujan',      label: 'Kemujan - Karimunjawa',      geom: kemujan },
  { id: 'rote02',       label: 'Rote 02',                    geom: rote02 },
  { id: 'rote03',       label: 'Rote 03',                    geom: rote03 },
  { id: 'parang',       label: 'Parang',                     geom: parang },
  { id: 'nias',         label: 'Nias',                       geom: nias },
  { id: 'teluk_bakau',  label: 'Teluk Bakau',                geom: teluk_bakau },
  { id: 'berakit',      label: 'Berakit',                    geom: berakit },
  { id: 'tana_merah',   label: 'Tana Merah',                 geom: tana_merah },
  { id: 'derawan',      label: 'Derawan',                    geom: derawan },
  { id: 'lembongan',    label: 'Lembongan',                  geom: lembongan },
  { id: 'tapil',        label: 'Tapil',                      geom: tapil },
];

// ==========================================================
// BUILD FEATURE COLLECTION
// ==========================================================
var roiFC = ee.FeatureCollection(
  ROI_LIST.map(function(r) {
    return ee.Feature(r.geom, { aoi_id: r.id, label: r.label });
  })
);

// ==========================================================
// DISPLAY
// ==========================================================
Map.setCenter(118, -4, 5);
Map.addLayer(roiFC,
  {color:'FF6D00', strokeWidth:2, fillColor:'FF6D0020'},
  'All ROIs', true);

// Print area for each ROI
print('=== ROI LIST (' + ROI_LIST.length + ' areas) ===');
ROI_LIST.forEach(function(r) {
  r.geom.area({maxError:100}).divide(1e6).evaluate(function(a) {
    var flag = a > 1000 ? ' WARNING: >1000 km2' : ' OK';
    print(r.id + ' | ' + r.label + ' | ' + a.toFixed(1) + ' km2' + flag);
  });
});

// ==========================================================
// EXPORT -- 1 file, all ROIs
// ==========================================================
Export.table.toAsset({
  collection:  roiFC,
  description: 'export_ROI_collection',
  assetId:     AOI_ASSET
});

print('');
print('Go to Tasks -> Run -> export_ROI_collection');
print('Asset: ' + AOI_ASSET);

// ==========================================================
// TRAINING POINTS -- show points where AGC_pred is not NA
// ==========================================================
var trainingPts = ee.FeatureCollection(
  ASSET_ROOT + '/chapter_3/TrainingPoint/trainingAGC_COMPLETE_for_GEE'
);

// Filter: AGC_pred not null AND not the literal string 'NA'
var ptsFiltered = trainingPts
  .filter(ee.Filter.notNull(['AGC_pred']))
  .filter(ee.Filter.neq('AGC_pred', 'NA'));

Map.addLayer(
  ptsFiltered,
  { color: '1E88E5', pointSize: 4, pointShape: 'circle' },
  'Training Points (AGC_pred valid)',
  true
);

print('Training points (AGC_pred valid):', ptsFiltered.size());
