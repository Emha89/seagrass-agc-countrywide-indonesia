/************************************************************ 
 * Script 09 - AGC TIME SERIES + FORECAST (GEE APPS)
 *            [v9: DEFF/ICC uncertainty method]
 *
 * CHANGE FROM v8r3:
 *   Uncertainty propagation switched from the spatial-autocorrelation
 *   (L=30, A_corr) method to the cluster design effect (DEFF/ICC)
 *   method, matching the adopted thesis methodology:
 *     Eq 3: SDtotal = sqrt(RMSE_yr^2 + SDdata^2)            [unchanged]
 *     Eq 7: SDC_adj = SDtotal x Apx_total / sqrt(N_EFFECTIVE) / 1e6
 *     Eq 8: CI half = 1.96 x SDC_adj                        [unchanged]
 *   N_EFFECTIVE = 201.3 is the pooled DEFF/ICC constant from
 *   recompute_N_chapter5_DEFF.R (33 DBSCAN clusters, ICC=0.188,
 *   DEFF=21.1) -- the same value used in 08_AGC_rf_modelUncer_indo.js
 *   (the national total AGC script). L and A_corr are removed.
 *
 * RETAINED FROM v8r3:
 *   - All UI layout, layer helpers, pixel inspector
 *   - loadAGC() with AGC_MIN/AGC_MAX range filter
 *   - Forecast functions (linear, recent3, recent4)
 *   - Debounce 600ms on sliders
 *   - Map layers reset ONLY on new polygon / re-analyse
 ************************************************************/

// ==========================================================
// 1. CONFIGURATION
// ==========================================================
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var SCALE         = 10;
var APPLY_YEARS   = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];
var FUTURE_START  = 2025;
var FUTURE_END    = 2027;
var CLAMP_MIN0    = true;
var MAX_AOI_KM2   = 1000;
var AGC_PALETTE   = ['#9e0142','#f46d43','#ffffbf','#66c2a5','#3288bd'];
var MASK_PALETTE  = ['#ffff00'];

// ----------------------------------------------------------
// AGC PREDICTION RANGE FILTER (gC/m2)
// Pixels outside [AGC_MIN, AGC_MAX] are masked before display
// and excluded from total AGC calculation.
// Rationale: model error is large below 10 and above 60 gC/m2
// (see F3_rf_evalModel_AGC.R BLOCK 5 / Figure S21).
// ----------------------------------------------------------
var AGC_MIN = 10;   // minimum valid AGC pred (gC/m2)
var AGC_MAX = 60;   // maximum valid AGC pred (gC/m2)

// ----------------------------------------------------------
// UNCERTAINTY PROPAGATION PARAMETERS (DEFF/ICC method)
// SDdata      : median SD from field AGC CI bounds
//               Source: calc_SDdata_AGC.R -> SDdata_AGC_for_GEE.csv
// N_EFFECTIVE : cluster design effect-adjusted effective sample size
//               Source: recompute_N_chapter5_DEFF.R (pooled, 33
//               DBSCAN clusters, ICC=0.188, DEFF=21.1, N=4,250)
// ----------------------------------------------------------
var SDdata      = 5.5046;   // gC/m2 -- field measurement uncertainty
var N_EFFECTIVE = 201.3;    // DEFF/ICC pooled effective sample size

var PROB_PREFIX   = ASSET_ROOT + '/chapter_3/output/02_PA/indo_RF_probability2_';

var YEARLY_STATS  = {
  2017: { n: 980,  RMSE: 8.724,  MAE: 6.279,  R2:  0.705 },
  2018: { n: 493,  RMSE: 14.693, MAE: 10.537, R2:  0.370 },
  2019: { n: 1396, RMSE: 11.451, MAE: 8.538,  R2:  0.511 },
  2020: { n: 10,   RMSE: 9.275,  MAE: 8.031,  R2:  0.541 },
  2021: { n: 2058, RMSE: 9.363,  MAE: 6.552,  R2:  0.647 },
  2022: { n: 972,  RMSE: 12.876, MAE: 9.768,  R2:  0.267 },
  2023: { n: 330,  RMSE: 10.785, MAE: 8.789,  R2:  0.599 }
};

var MC_SUMMARY = {
  RMSE_mean: 8.53, RMSE_sd: 0.03, RMSE_lower: 8.47, RMSE_upper: 8.57,
  MAE_mean:  6.26, MAE_sd:  0.02, MAE_lower:  6.22, MAE_upper:  6.29,
  R2_mean:   0.79, R2_sd:   0.00, R2_lower:   0.79, R2_upper:   0.79
};

// ==========================================================
// 2. STATIC DATA ASSETS
// ==========================================================
var AGC_PREFIX = ASSET_ROOT + '/chapter_3/output/07_AGC/indo_RF_AGCsimulated_reduceVar_limitCI2_';

var MPA_FC     = ee.FeatureCollection(
  ASSET_ROOT + '/chapter_2/InputArea/ZonasiKawasanKonservasi_AR_50K_2021_Rev2022_diss'
);
var depthRaw   = ee.Image(
  ASSET_ROOT + '/chapter_2/InputArea/raster_INA_ACA_Bathymetry'
);
var depthMask  = depthRaw.unmask(0).gt(0).and(depthRaw.unmask(0).lte(500));
var pixelArea  = ee.Image.pixelArea();

var probImageCollection = ee.ImageCollection(
  APPLY_YEARS.map(function(y) {
    return ee.Image(PROB_PREFIX + y).rename('prob');
  })
);
var maxProbImage = probImageCollection.reduce(ee.Reducer.max()).rename('max_probability');

var trainingAllFC = ee.FeatureCollection(
  ASSET_ROOT + '/chapter_3/TrainingPoint/training_COMPLETE_for_GEE'
);

// ==========================================================
// 3. STATE VARIABLES
// ==========================================================
var currentAOI  = null;    // ee.Geometry of latest valid polygon
var inspectMode = false;   // false = Draw active, true = Inspect active

// Pixel inspector references - updated by renderResults()
var _pixelInfoLabel  = null;
var _pixelChartPanel = null;
var _dynamicMask     = null;

// ==========================================================
// 4. DYNAMIC MASK BUILDER
// ==========================================================
function buildDynamicMask(threshold) {
  return maxProbImage
    .gte(ee.Number(threshold))
    .rename('seagrass_presence')
    .selfMask();
}

// ==========================================================
// 5. ROOT LAYOUT
// ==========================================================
ui.root.clear();

var leftPanel = ui.Panel({
  style: { width: '20%', padding: '10px', backgroundColor: '#ffffff' }
});
var mapPanel = ui.Map();
mapPanel.setOptions('SATELLITE');
mapPanel.setCenter(118, -2, 5);
mapPanel.style().set({ width: '50%' });
mapPanel.setControlVisibility({ layerList: false });

var rightPanel = ui.Panel({
  style: { width: '30%', padding: '10px', backgroundColor: '#fafafa' }
});

ui.root.add(ui.Panel({
  widgets: [leftPanel, mapPanel, rightPanel],
  layout: ui.Panel.Layout.flow('horizontal'),
  style: { stretch: 'both' }
}));

// ==========================================================
// 6. BUTTON STYLE HELPERS
// ==========================================================
var STYLE_BTN_ACTIVE = {
  width: '100%', margin: '0 0 4px 0',
  backgroundColor: '#212121', color: '#111111'
};
var STYLE_BTN_INACTIVE = {
  width: '100%', margin: '0 0 4px 0',
  backgroundColor: '#9e9e9e', color: '#111111'
};
var STYLE_BTN_DANGER = {
  width: '100%', margin: '0 0 4px 0',
  backgroundColor: '#e65100', color: '#111111'
};

var togglePanel = ui.Panel({
  style: { margin: '0 0 4px 0', padding: '0',
           backgroundColor: 'rgba(0,0,0,0)' }
});

function secLabel(txt) {
  return ui.Label(txt, {
    fontWeight:'bold', fontSize:'11px', color:'#333', margin:'4px 0 2px 0'
  });
}

// ==========================================================
// 7. LEFT PANEL - static widgets
// ==========================================================
leftPanel.add(ui.Label('AGC Viewer', {
  fontWeight: 'bold', fontSize: '16px', color: '#111111', margin: '4px 0 1px 0'
}));
leftPanel.add(ui.Label('Seagrass Carbon - Indonesia', {
  fontSize: '10px', color: '#777', margin: '0 0 1px 0'
}));
leftPanel.add(ui.Label('Draw polygon AOI - max ' + MAX_AOI_KM2 + ' km2', {
  fontSize: '9px', color: '#bbb', fontStyle: 'italic', margin: '0 0 7px 0'
}));
leftPanel.add(ui.Label('------------------------------', { color:'#ddd', margin:'0 0 6px 0' }));

leftPanel.add(secLabel('Cursor Mode:'));
leftPanel.add(ui.Label(
  '1. Zoom to target area\n2. Activate Draw -> click vertices\n3. Click first vertex again to close',
  { fontSize:'9px', color:'#555', margin:'0 0 4px 2px', whiteSpace:'pre' }
));
leftPanel.add(ui.Label(
  'Tip: at zoom ~9, one screen ~ 500-1000 km2',
  { fontSize:'9px', color:'#1565c0', fontStyle:'italic', margin:'0 0 6px 2px' }
));

leftPanel.add(togglePanel);

var btnClear = ui.Button({
  label: 'Clear & Redraw',
  style: STYLE_BTN_DANGER
});
leftPanel.add(btnClear);

var aoiStatusLabel = ui.Label('AOI: not drawn yet', {
  fontSize: '10px', color: '#e65100', fontWeight: 'bold',
  margin: '2px 0 6px 0'
});
leftPanel.add(aoiStatusLabel);
leftPanel.add(ui.Label('------------------------------', { color:'#ddd', margin:'0 0 6px 0' }));

leftPanel.add(secLabel('Year (S2 basemap):'));
var yearValueLabel = ui.Label('2024', {
  fontSize:'12px', fontWeight:'bold', color:'#111111', margin:'0 0 2px 0', textAlign:'center'
});
leftPanel.add(yearValueLabel);
var yearSlider = ui.Slider({
  min: APPLY_YEARS[0], max: APPLY_YEARS[APPLY_YEARS.length - 1],
  value: 2024, step: 1, style: { width:'100%', margin:'0 0 8px 0' }
});
leftPanel.add(yearSlider);

leftPanel.add(ui.Label('------------------------------', { color:'#ddd', margin:'0 0 4px 0' }));
leftPanel.add(secLabel('Min Prob Threshold (Seagrass Mask):'));
leftPanel.add(ui.Label(
  'Max prob 2017-2024 >= threshold -> seagrass pixel',
  { fontSize:'9px', color:'#888', fontStyle:'italic', margin:'0 0 3px 0' }
));
var probValueLabel = ui.Label('0.50', {
  fontSize:'13px', fontWeight:'bold', color:'#1565c0',
  margin:'0 0 2px 0', textAlign:'center'
});
leftPanel.add(probValueLabel);
var probSlider = ui.Slider({
  min: 0.00, max: 1.00, value: 0.50, step: 0.05,
  style: { width:'100%', margin:'0 0 4px 0' }
});
leftPanel.add(probSlider);
leftPanel.add(ui.Panel({
  widgets: [
    ui.Label('0.00', { fontSize:'8px', color:'#888', margin:'0' }),
    ui.Label('',     { stretch:'horizontal' }),
    ui.Label('1.00', { fontSize:'8px', color:'#888', margin:'0' })
  ],
  layout: ui.Panel.Layout.flow('horizontal'),
  style: { margin:'0 0 4px 0' }
}));
var areaEstLabel = ui.Label('Seagrass area: -- km2', {
  fontSize:'10px', color:'#333', fontWeight:'bold',
  backgroundColor:'#edf4fb', padding:'4px 6px',
  margin:'2px 0 4px 0'
});
leftPanel.add(areaEstLabel);
leftPanel.add(ui.Label(
  'AGC filter: ' + AGC_MIN + ' - ' + AGC_MAX + ' gC/m2',
  { fontSize:'9px', color:'#555', fontStyle:'italic',
    backgroundColor:'#fff8e1', padding:'3px 6px',
    margin:'0 0 8px 0' }
));
leftPanel.add(ui.Label('------------------------------', { color:'#ddd', margin:'0 0 6px 0' }));

leftPanel.add(secLabel('Optional Layers:'));
var showMPACheck = ui.Checkbox({
  label: 'Show MPA boundaries', value: false,
  style: { fontSize:'11px', color:'#444', margin:'0 0 4px 3px' }
});
leftPanel.add(showMPACheck);
var showTrainingCheck = ui.Checkbox({
  label: 'Show training points', value: false,
  style: { fontSize:'11px', color:'#444', margin:'0 0 4px 3px' }
});
leftPanel.add(showTrainingCheck);
var showAOICheck = ui.Checkbox({
  label: 'Show AOI outline', value: true,
  style: { fontSize:'11px', color:'#444', margin:'0 0 4px 3px' }
});
leftPanel.add(showAOICheck);
var showAGCCheck = ui.Checkbox({
  label: 'Show AGC layer', value: true,
  style: { fontSize:'11px', color:'#444', margin:'0 0 4px 3px' }
});
leftPanel.add(showAGCCheck);
var showMaskCheck = ui.Checkbox({
  label: 'Show Seagrass Mask (>= 0.50)', value: false,
  style: { fontSize:'11px', color:'#444', margin:'0 0 8px 3px' }
});
leftPanel.add(showMaskCheck);

function updateMaskCheckLabel(threshold) {
  showMaskCheck.setLabel('Show Seagrass Mask (>= ' + threshold.toFixed(2) + ')');
}
function updateAGCCheckLabel(year) {
  showAGCCheck.setLabel('Show AGC - ' + year);
}

leftPanel.add(ui.Label('------------------------------', { color:'#ddd', margin:'0 0 6px 0' }));

var statusLabel = ui.Label('Status: ready - zoom & draw AOI', {
  fontSize:'10px', color:'#111111', margin:'0 0 6px 0'
});
leftPanel.add(statusLabel);
leftPanel.add(ui.Label('------------------------------', { color:'#ddd', margin:'0 0 6px 0' }));

leftPanel.add(secLabel('AGC Legend:'));
var legendBox = ui.Panel({
  style:{padding:'0', margin:'0 0 4px 0', backgroundColor:'rgba(0,0,0,0)'}
});
legendBox.add(ui.Label('AGC (gC/m2) - fixed range: ' + AGC_MIN + ' - ' + AGC_MAX, {
  fontSize:'9px', color:'#111111', fontWeight:'bold', margin:'0 0 3px 0'
}));
function hexToRgb(hex) {
  var h = hex.replace('#','');
  return [parseInt(h.slice(0,2),16), parseInt(h.slice(2,4),16), parseInt(h.slice(4,6),16)];
}
function rgbToHex(r,g,b) {
  return '#'+('0'+Math.round(r).toString(16)).slice(-2)
            +('0'+Math.round(g).toString(16)).slice(-2)
            +('0'+Math.round(b).toString(16)).slice(-2);
}
var barRow = ui.Panel({ layout:ui.Panel.Layout.flow('horizontal'), style:{margin:'0 0 2px 0'} });
var N_SLICES = 36;
for (var si = 0; si < N_SLICES; si++) {
  var t    = si / (N_SLICES - 1);
  var pidx = Math.min(Math.floor(t*(AGC_PALETTE.length-1)), AGC_PALETTE.length-2);
  var frac = t*(AGC_PALETTE.length-1) - pidx;
  var c1   = hexToRgb(AGC_PALETTE[pidx]);
  var c2   = hexToRgb(AGC_PALETTE[pidx+1]);
  barRow.add(ui.Label('', {
    backgroundColor: rgbToHex(c1[0]+(c2[0]-c1[0])*frac,
                               c1[1]+(c2[1]-c1[1])*frac,
                               c1[2]+(c2[2]-c1[2])*frac),
    width:'5px', height:'12px', margin:'0', padding:'0'
  }));
}
legendBox.add(barRow);
legendBox.add(ui.Panel({
  widgets: [
    ui.Label(String(AGC_MIN) + ' gC/m2', { fontSize:'9px', color:'#555', margin:'0' }),
    ui.Label('',                           { stretch:'horizontal', margin:'0', padding:'0' }),
    ui.Label(String(AGC_MAX) + ' gC/m2', { fontSize:'9px', color:'#555', margin:'0' })
  ],
  layout: ui.Panel.Layout.flow('horizontal'), style:{ margin:'1px 0 0 0' }
}));
legendBox.add(ui.Label('Consistent scale across all years & AOIs', {
  fontSize:'8px', color:'#aaa', fontStyle:'italic', margin:'2px 0 0 0'
}));
leftPanel.add(legendBox);

leftPanel.add(ui.Label('------------------------------', { color:'#ddd', margin:'6px 0 6px 0' }));
leftPanel.add(ui.Label('MC Bootstrap (n=100)', {
  fontWeight:'bold', fontSize:'12px', color:'#111111', margin:'0 0 4px 0'
}));
var mc  = MC_SUMMARY;
var mcW = ['44px','38px','30px','90px'];
leftPanel.add(ui.Panel({
  widgets: ['Metric','Mean','SD','95% CI'].map(function(h,i){
    return ui.Label(h,{fontSize:'11px',fontWeight:'bold',color:'#fff',
                        width:mcW[i],margin:'1px 1px',
                        backgroundColor:'#2171b5',padding:'3px 3px'});
  }),
  layout: ui.Panel.Layout.flow('horizontal'), style:{margin:'0 0 1px 0'}
}));
[
  ['RMSE', mc.RMSE_mean.toFixed(2), mc.RMSE_sd.toFixed(2), mc.RMSE_lower.toFixed(2)+' - '+mc.RMSE_upper.toFixed(2)],
  ['MAE',  mc.MAE_mean.toFixed(2),  mc.MAE_sd.toFixed(2),  mc.MAE_lower.toFixed(2) +' - '+mc.MAE_upper.toFixed(2)],
  ['R2',   mc.R2_mean.toFixed(2),   mc.R2_sd.toFixed(2),   mc.R2_lower.toFixed(2)  +' - '+mc.R2_upper.toFixed(2)]
].forEach(function(row, ri){
  var bg = (ri%2===0)?'#edf4fb':'#ffffff';
  leftPanel.add(ui.Panel({
    widgets: row.map(function(cell,ci){
      return ui.Label(cell,{fontSize:'11px',color:'#111',fontFamily:'monospace',
                             width:mcW[ci],margin:'1px 1px',
                             backgroundColor:bg,padding:'3px 3px'});
    }),
    layout: ui.Panel.Layout.flow('horizontal'), style:{margin:'0 0 1px 0'}
  }));
});
leftPanel.add(ui.Label('Source: F1_rf_model_AGC.R',{
  fontSize:'8px', color:'#aaa', fontStyle:'italic', margin:'1px 0 6px 0'
}));

leftPanel.add(ui.Label('------------------------------', { color:'#ddd', margin:'0 0 6px 0' }));
leftPanel.add(ui.Label('Model Validation (from R)', {
  fontWeight:'bold', fontSize:'12px', color:'#111111', margin:'0 0 4px 0'
}));
var hW_l = ['36px','40px','44px','40px','34px'];
leftPanel.add(ui.Panel({
  widgets: ['Year','n','RMSE','MAE','R2'].map(function(h,i){
    return ui.Label(h,{fontSize:'11px',fontWeight:'bold',color:'#fff',
                        width:hW_l[i],margin:'1px 1px',
                        backgroundColor:'#2171b5',padding:'3px 3px'});
  }),
  layout: ui.Panel.Layout.flow('horizontal'), style:{margin:'0 0 1px 0'}
}));
[2017,2018,2019,2020,2021,2022,2023].forEach(function(yr){
  var s   = YEARLY_STATS[yr];
  var bg  = (yr%2===0)?'#edf4fb':'#ffffff';
  var r2c = s.R2>=0.6?'#1565c0':(s.R2>=0.4?'#e65100':'#c62828');
  leftPanel.add(ui.Panel({
    widgets: [
      ui.Label(String(yr),        {fontSize:'11px',color:'#111',fontFamily:'monospace',width:hW_l[0],margin:'1px 1px',backgroundColor:bg,padding:'3px 3px'}),
      ui.Label(String(s.n),       {fontSize:'11px',color:'#333',fontFamily:'monospace',width:hW_l[1],margin:'1px 1px',backgroundColor:bg,padding:'3px 3px'}),
      ui.Label(s.RMSE.toFixed(2), {fontSize:'11px',color:'#111',fontFamily:'monospace',width:hW_l[2],margin:'1px 1px',backgroundColor:bg,padding:'3px 3px'}),
      ui.Label(s.MAE.toFixed(2),  {fontSize:'11px',color:'#333',fontFamily:'monospace',width:hW_l[3],margin:'1px 1px',backgroundColor:bg,padding:'3px 3px'}),
      ui.Label(s.R2.toFixed(3),   {fontSize:'11px',color:r2c,fontFamily:'monospace',fontWeight:'bold',width:hW_l[4],margin:'1px 1px',backgroundColor:bg,padding:'3px 3px'})
    ],
    layout: ui.Panel.Layout.flow('horizontal'), style:{margin:'0 0 1px 0'}
  }));
});
leftPanel.add(ui.Label('Unit: gC/m2  |  Source: F3_rf_evalModel_AGC.R',{
  fontSize:'8px', color:'#aaa', fontStyle:'italic', margin:'1px 0 0 0'
}));

// ==========================================================
// 8. RIGHT PANEL
// ==========================================================
rightPanel.add(ui.Label('Analysis Results', {
  fontWeight:'bold', fontSize:'15px', color:'#111111', margin:'4px 0 2px 0'
}));
rightPanel.add(ui.Label('Forecast Method:', {
  fontWeight:'bold', fontSize:'11px', color:'#333', margin:'4px 0 2px 0'
}));
var forecastSelect = ui.Select({
  items: [
    { label:'Linear trend - all years (2017-2024)',        value:'linear'  },
    { label:'Recent trend - last 3 years only (2022-2024)',value:'recent3' },
    { label:'Recent trend - last 4 years only (2021-2024)',value:'recent4' }
  ],
  value:'linear', style:{ width:'100%' }
});
rightPanel.add(ui.Panel({ widgets:[forecastSelect], style:{
  backgroundColor:'rgba(0,0,0,0)', padding:'0', margin:'0 0 8px 0'
}}));
var resultsPanel = ui.Panel({ style:{ margin:'0', padding:'0' } });
rightPanel.add(resultsPanel);
resultsPanel.add(ui.Label('Zoom to target area & draw a polygon AOI', {
  fontSize:'11px', color:'#aaa', fontStyle:'italic', margin:'6px 0'
}));

// ==========================================================
// 9. MAP LAYER HELPERS
// ==========================================================
function loadAGC(year, AOI, dynamicMask) {
  var img = ee.Image(AGC_PREFIX + year)
    .select('AGC_pred')
    .rename('AGC_mean');
  var rangeMask = img.gte(AGC_MIN).and(img.lte(AGC_MAX));
  return img
    .updateMask(rangeMask)
    .updateMask(dynamicMask)
    .updateMask(depthMask)
    .clip(AOI);
}

function buildS2Mosaic(year, AOI) {
  function maskS2(img) {
    var scl = img.select('SCL');
    return img.updateMask(
      scl.neq(1).and(scl.neq(3)).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10))
    );
  }
  return ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(AOI)
    .filterDate(ee.Date.fromYMD(year,1,1),
                ee.Date.fromYMD(ee.Number(year).add(1),1,1))
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
    .map(maskS2)
    .select(['B4','B3','B2'])
    .reduce(ee.Reducer.percentile([60]))
    .rename(['B4','B3','B2'])
    .updateMask(depthMask)
    .clip(AOI);
}

function addS2Layer(year, AOI) {
  mapPanel.addLayer(
    buildS2Mosaic(year, AOI),
    { bands:['B4','B3','B2'], min:0, max:3000, gamma:1.4 },
    'S2 p60 - ' + year, true
  );
}

function addAGCLayer(agcImg, year, visible) {
  mapPanel.addLayer(agcImg,
    { min: AGC_MIN, max: AGC_MAX, palette: AGC_PALETTE },
    'AGC - ' + year, visible
  );
}

function addMaskLayer(dynamicMask, AOI, threshold) {
  mapPanel.addLayer(
    dynamicMask.clip(AOI),
    { palette: MASK_PALETTE, opacity: 0.45 },
    'Seagrass mask (MaxProb >= ' + threshold.toFixed(2) + ')',
    showMaskCheck.getValue()
  );
}

function addMPALayer() {
  mapPanel.addLayer(
    MPA_FC,
    { color:'00dd66', strokeWidth:1, fillColor:'00000000' },
    'MPA Boundaries',
    showMPACheck.getValue()
  );
}

function addTrainingLayer() {
  mapPanel.addLayer(
    trainingAllFC,
    { color:'00cc44', pointSize:3, pointShape:'circle', opacity:0.8 },
    'Training points (2017-2023)',
    showTrainingCheck.getValue()
  );
}

function addAOIOutlineLayer(aoiGeom) {
  mapPanel.addLayer(
    ee.FeatureCollection([ee.Feature(aoiGeom)]),
    { color:'FF6D00', strokeWidth:2, fillColor:'FF6D0020' },
    'AOI (user drawn)',
    showAOICheck.getValue()
  );
}

// ==========================================================
// 10. AGC AGGREGATION - DEFF/ICC method
// ----------------------------------------------------------
// Eq 3: SDtotal = sqrt(RMSE_yr^2 + SDdata^2)
//         RMSE_yr: per-year model error from YEARLY_STATS
//         SDdata : 5.5046 gC/m2 (median field CI SD)
// Eq 7: SDC_adj  = SDtotal x Apx_total / sqrt(N_EFFECTIVE) / 1e6  [ton C]
//         Apx_total = N_pixels x SCALE^2  (total seagrass area, m2)
//         N_EFFECTIVE = 201.3 (DEFF/ICC, pooled -- see CONFIGURATION)
// Eq 8: CI half = 1.96 x SDC_adj
//
// NOTE: RMSE_yr is used as a spatially uniform proxy for
// pixel-level SDmodel, as pixel-level uncertainty assets
// were not pre-computed during national-scale deployment.
// This matches the formula in 08_AGC_rf_modelUncer_indo.js
// (the national total AGC script) exactly, just evaluated
// here per user-drawn AOI instead of the full country.
// ==========================================================
function aggregateYear(year, agcImg, AOI) {

  // Single reduceRegion: sum (-> totalC) + count (-> N_pixels)
  var stats = agcImg.rename('AGC').reduceRegion({
    reducer:    ee.Reducer.sum().combine(ee.Reducer.count(), null, true),
    geometry:   AOI,
    scale:      10,
    maxPixels:  1e13,
    bestEffort: true,
    tileScale:  16
  });

  var totalC_ton = ee.Number(stats.get('AGC_sum')).multiply(SCALE * SCALE).divide(1e6);
  var N_pixels   = ee.Number(stats.get('AGC_count'));
  var area_m2    = N_pixels.multiply(SCALE * SCALE);

  // Eq 3: SDtotal per pixel (gC/m2)
  var RMSE_yr  = YEARLY_STATS[year] ? YEARLY_STATS[year].RMSE : MC_SUMMARY.RMSE_mean;
  var SDtotal  = Math.sqrt(RMSE_yr * RMSE_yr + SDdata * SDdata);

  // Eq 7: DEFF/ICC-adjusted area-level SD (ton C)
  var SDC_adj_ton = area_m2.multiply(SDtotal).divide(Math.sqrt(N_EFFECTIVE)).divide(1e6);

  // Eq 8: 95% CI half-width (ton C)
  var half = SDC_adj_ton.multiply(1.96);

  return ee.Feature(null, {
    year:            year,
    totalCarbon_ton: totalC_ton,
    CI_lower_adj:    totalC_ton.subtract(half).max(0),
    CI_upper_adj:    totalC_ton.add(half),
    CI_error_adj:    half
  });
}

// ==========================================================
// 11. FORECAST FUNCTIONS
// ==========================================================
function clamp0(x) {
  return ee.Number(ee.Algorithms.If(CLAMP_MIN0, ee.Number(x).max(0), x));
}
function forecast_linear(fc, yS, yE) {
  var lr = fc.map(function(f){ return f.set('constant',1); })
             .reduceColumns(ee.Reducer.linearRegression(2,1),
                            ['constant','year','totalCarbon_ton']);
  var B = ee.Array(lr.get('coefficients'));
  return ee.FeatureCollection(ee.List.sequence(yS,yE).map(function(y){
    y = ee.Number(y);
    return ee.Feature(null,{
      year:y, value:clamp0(B.get([0,0]).add(B.get([1,0]).multiply(y))), series:'forecast'
    });
  }));
}
function forecast_recent(fc, k, yS, yE) {
  var cut = ee.Number(fc.aggregate_max('year')).subtract(k-1);
  var lr  = fc.filter(ee.Filter.gte('year',cut))
              .map(function(f){ return f.set('constant',1); })
              .reduceColumns(ee.Reducer.linearRegression(2,1),
                             ['constant','year','totalCarbon_ton']);
  var B = ee.Array(lr.get('coefficients'));
  return ee.FeatureCollection(ee.List.sequence(yS,yE).map(function(y){
    y = ee.Number(y);
    return ee.Feature(null,{
      year:y, value:clamp0(B.get([0,0]).add(B.get([1,0]).multiply(y))), series:'forecast'
    });
  }));
}
function getForecast(method, tsFC) {
  if (method==='recent3') return forecast_recent(tsFC, 3, FUTURE_START, FUTURE_END);
  if (method==='recent4') return forecast_recent(tsFC, 4, FUTURE_START, FUTURE_END);
  return forecast_linear(tsFC, FUTURE_START, FUTURE_END);
}

// ==========================================================
// 12. PIXEL INSPECTOR onClick handler
// ==========================================================
function pixelInspectorClick(coords) {
  if (!_pixelInfoLabel || !_pixelChartPanel || !_dynamicMask) return;

  var pt = ee.Geometry.Point([coords.lon, coords.lat]);
  _pixelInfoLabel.setValue(
    'Lat: ' + coords.lat.toFixed(5) + '   Lon: ' + coords.lon.toFixed(5)
  );

  var layers = mapPanel.layers();
  for (var li = 0; li < layers.length(); li++) {
    if ((layers.get(li).getName() || '') === 'Inspect point') {
      layers.remove(layers.get(li));
      break;
    }
  }
  mapPanel.addLayer(
    ee.FeatureCollection([ee.Feature(pt)]),
    { color: 'FF0000', pointSize: 6, pointShape: 'cross', strokeWidth: 2 },
    'Inspect point', true
  );

  _pixelChartPanel.clear();
  _pixelChartPanel.add(ui.Label('Sampling pixel...', {
    fontSize:'10px', color:'#e65100', margin:'0 0 4px 0'
  }));

  var ptFC = ee.FeatureCollection(
    APPLY_YEARS.map(function(year) {
      var val = ee.Image(AGC_PREFIX + year)
                  .select('AGC_pred')
                  .updateMask(_dynamicMask)
                  .updateMask(depthMask)
                  .reduceRegion({
                    reducer: ee.Reducer.first(),
                    geometry: pt,
                    scale: 10, maxPixels: 1
                  }).values().get(0);
      return ee.Feature(null, { year: year, AGC: val });
    })
  );

  ptFC.evaluate(function(res) {
    _pixelChartPanel.clear();
    if (!res || !res.features) {
      _pixelChartPanel.add(ui.Label('Could not sample pixel.', {
        fontSize:'10px', color:'#c62828'
      }));
      return;
    }
    var hasData = res.features.some(function(f){ return f.properties.AGC !== null; });
    if (!hasData) {
      _pixelChartPanel.add(ui.Label('No seagrass at this pixel (outside mask).', {
        fontSize:'10px', color:'#c62828'
      }));
      return;
    }
    var chartFC = ee.FeatureCollection(
      res.features.map(function(f){
        return ee.Feature(null,{ year:f.properties.year, AGC:f.properties.AGC||0 });
      })
    );
    var pixelChart = ui.Chart.feature.byFeature(chartFC, 'year', ['AGC'])
      .setChartType('ColumnChart')
      .setOptions({
        title: 'Pixel AGC - ' + coords.lat.toFixed(4) + ', ' + coords.lon.toFixed(4),
        backgroundColor:'#ffffff',
        titleTextStyle:{ color:'#222', fontSize:10, bold:true },
        hAxis:{ title:'Year', format:'####', titleTextStyle:{color:'#444'}, textStyle:{color:'#333'} },
        vAxis:{ title:'AGC (gC/m2)', titleTextStyle:{color:'#444'}, textStyle:{color:'#333'}, minValue:0 },
        legend:{ position:'none' },
        colors:['#2171b5'],
        chartArea:{ backgroundColor:'#f8f9fa', width:'82%', height:'65%' }
      });
    pixelChart.style().set({ width:'100%', margin:'0 0 6px 0' });
    _pixelChartPanel.add(pixelChart);
  });
}

// ==========================================================
// 13. TOGGLE MODE
// ==========================================================
function toggleMode(wantInspect) {
  if (wantInspect && !currentAOI) {
    statusLabel.setValue('Draw an AOI first before using Inspect');
    return;
  }

  inspectMode = wantInspect;
  togglePanel.clear();

  var btnDrawNew = ui.Button({
    label: inspectMode ? 'Draw AOI' : 'Draw AOI (ACTIVE)',
    style: inspectMode ? STYLE_BTN_INACTIVE : STYLE_BTN_ACTIVE
  });
  var btnInspectNew = ui.Button({
    label: inspectMode ? 'Inspect Pixel (ACTIVE)' : 'Inspect Pixel',
    style: inspectMode ? STYLE_BTN_ACTIVE : STYLE_BTN_INACTIVE
  });

  btnDrawNew.onClick(function() { toggleMode(false); });
  btnInspectNew.onClick(function() { toggleMode(true); });

  togglePanel.add(btnDrawNew);
  togglePanel.add(btnInspectNew);

  if (inspectMode) {
    mapPanel.drawingTools().setShown(false);
    statusLabel.setValue('Mode: Inspect - click a pixel to plot AGC');
  } else {
    mapPanel.drawingTools().setShown(true);
    mapPanel.drawingTools().setDrawModes(['polygon']);
    mapPanel.drawingTools().setLinked(false);
    clearInspectMarker();
    statusLabel.setValue('Mode: Draw AOI - click vertices on the map');
  }
}

// ==========================================================
// 14. RIGHT PANEL - RESULTS RENDERING
// ==========================================================
function renderResults(method, result, tsFC, dynamicMask, AOI, threshold) {
  resultsPanel.clear();

  _dynamicMask = dynamicMask;

  var fcActual   = tsFC.map(function(f){
    return ee.Feature(null,{year:f.get('year'),value:f.get('totalCarbon_ton'),series:'actual'});
  });
  var fcForecast = getForecast(method, tsFC);
  var xTicks     = APPLY_YEARS.concat([FUTURE_START, FUTURE_START+1, FUTURE_END]);
  var methodLabel = {
    'linear':'Linear trend - all years',
    'recent3':'Recent trend - last 3 yrs',
    'recent4':'Recent trend - last 4 yrs'
  }[method] || method;

  var chart = ui.Chart.feature.groups(
      fcActual.merge(fcForecast).sort('year'), 'year', 'value', 'series'
    )
    .setChartType('LineChart')
    .setOptions({
      title: methodLabel + '  [mask >= ' + threshold.toFixed(2) + ']',
      backgroundColor:'#ffffff',
      titleTextStyle:{ color:'#222', fontSize:11, bold:true },
      hAxis:{ title:'Year', ticks:xTicks, format:'####',
              titleTextStyle:{color:'#444'}, textStyle:{color:'#333'},
              viewWindow:{min:xTicks[0], max:xTicks[xTicks.length-1]} },
      vAxis:{ title:'Total AGC (tC)',
              titleTextStyle:{color:'#444'}, textStyle:{color:'#333'} },
      legend:{ position:'bottom', textStyle:{color:'#333'} },
      pointSize:5, lineWidth:2,
      series:{
        0:{ color:'#1565c0', lineWidth:2 },
        1:{ color:'#c62828', lineDashStyle:[4,4], lineWidth:2 }
      },
      chartArea:{ backgroundColor:'#f8f9fa', width:'82%', height:'65%' }
    });
  chart.style().set({ width:'100%', margin:'0 0 8px 0' });
  resultsPanel.add(chart);

  resultsPanel.add(ui.Label('------------------------------',{color:'#ddd',margin:'2px 0 4px 0'}));
  resultsPanel.add(ui.Label('Total AGC per Year', {
    fontWeight:'bold', fontSize:'13px', color:'#111111', margin:'4px 0 2px 0'
  }));
  resultsPanel.add(ui.Label('Seagrass mask: MaxProb >= ' + threshold.toFixed(2), {
    fontSize:'10px', color:'#1565c0', fontWeight:'bold', margin:'0 0 2px 0'
  }));
  resultsPanel.add(ui.Label('AGC filter: ' + AGC_MIN + ' - ' + AGC_MAX + ' gC/m2', {
    fontSize:'10px', color:'#e65100', fontWeight:'bold', margin:'0 0 4px 0'
  }));
  var tsW = ['40px','72px','160px'];
  resultsPanel.add(ui.Panel({
    widgets:[
      ui.Label('Year',       {fontSize:'12px',fontWeight:'bold',color:'#fff',width:tsW[0],margin:'1px 2px',backgroundColor:'#2171b5',padding:'3px 4px'}),
      ui.Label('tC',         {fontSize:'12px',fontWeight:'bold',color:'#fff',width:tsW[1],margin:'1px 2px',backgroundColor:'#2171b5',padding:'3px 4px'}),
      ui.Label('95% CI (tC)',{fontSize:'12px',fontWeight:'bold',color:'#fff',width:tsW[2],margin:'1px 2px',backgroundColor:'#2171b5',padding:'3px 4px'})
    ],
    layout:ui.Panel.Layout.flow('horizontal'), style:{margin:'0 0 1px 0'}
  }));
  result.features.forEach(function(f) {
    var p  = f.properties;
    var bg = (p.year%2===0)?'#edf4fb':'#ffffff';
    resultsPanel.add(ui.Panel({
      widgets:[
        ui.Label(String(p.year),
          {fontSize:'12px',color:'#111',fontFamily:'monospace',width:tsW[0],margin:'1px 2px',backgroundColor:bg,padding:'3px 4px'}),
        ui.Label(p.totalCarbon_ton.toFixed(1),
          {fontSize:'12px',color:'#111',fontFamily:'monospace',width:tsW[1],margin:'1px 2px',backgroundColor:bg,padding:'3px 4px'}),
        ui.Label('['+p.CI_lower_adj.toFixed(1)+' - '+p.CI_upper_adj.toFixed(1)+']',
          {fontSize:'12px',color:'#333',fontFamily:'monospace',width:tsW[2],margin:'1px 2px',backgroundColor:bg,padding:'3px 4px'})
      ],
      layout:ui.Panel.Layout.flow('horizontal'), style:{margin:'0 0 1px 0'}
    }));
  });
  resultsPanel.add(ui.Label('CI: Eq 3/7/8 DEFF-ICC propagation - SDdata=5.50 gC/m2, N_eff=201.3 - unit: tC',{
    fontSize:'9px', color:'#888', fontStyle:'italic', margin:'1px 0 6px 0'
  }));

  resultsPanel.add(ui.Label('------------------------------',{color:'#ddd',margin:'2px 0 4px 0'}));
  resultsPanel.add(ui.Label('Pixel AGC Inspector', {
    fontWeight:'bold', fontSize:'13px', color:'#111111', margin:'2px 0 2px 0'
  }));
  resultsPanel.add(ui.Label('Switch to Inspect mode, then click a pixel on the map.', {
    fontSize:'10px', color:'#666', fontStyle:'italic', margin:'0 0 4px 0'
  }));

  _pixelInfoLabel = ui.Label('No pixel selected yet.', {
    fontSize:'10px', color:'#aaa', margin:'0 0 4px 0'
  });
  _pixelChartPanel = ui.Panel({ style:{margin:'0', padding:'0'} });

  resultsPanel.add(_pixelInfoLabel);
  resultsPanel.add(_pixelChartPanel);

  resultsPanel.add(ui.Label('------------------------------',{color:'#ddd',margin:'4px 0 6px 0'}));
  var descPanel = ui.Panel({ style:{ backgroundColor:'rgba(0,0,0,0)', padding:'0', margin:'0 0 4px 0' }});
  descPanel.add(ui.Label('About This Map', {
    fontWeight:'bold', fontSize:'12px', color:'#111111', margin:'0 0 6px 0'
  }));
  [
    'These maps are part of PhD research of m.hafizt@uq.edu.au, funded by The Indonesia Endowment Fund (LPDP) and Supervised by: c.roelfsema@uq.edu.au, s.phinn@uq.edu.au, mitchell.lyons@unsw.edu.au, k.mcmahon@ecu.edu.au.',
    '',
    'Data sources:',
    '1. Google Embedding Satellite (2017-2024)',
    '2. Sentinel-2A satellite imagery (2017-2024)',
    '3. https://allencoralatlas.org/',
    '4. Field training data: BRIN, UGM, KKP, ACA',
    '5. https://portaldata.kkp.go.id/ (MPA 50K Rev2022)'
  ].forEach(function(line) {
    descPanel.add(ui.Label(line, {
      fontSize:   line===''?'4px':'11px',
      color:      '#111111',
      fontWeight: line==='Data sources:'?'bold':'normal',
      margin:     line===''?'0':'0 0 2px 0'
    }));
  });
  resultsPanel.add(descPanel);
}

// ==========================================================
// 15. CORE ANALYSIS
// ==========================================================
function runAnalysis(aoiGeom) {
  if (!aoiGeom) return;

  var method      = forecastSelect.getValue() || 'recent4';
  var yearVal     = yearSlider.getValue()     || 2024;
  var threshold   = probSlider.getValue()     || 0.50;
  var AOI         = aoiGeom;
  var dynamicMask = buildDynamicMask(threshold);

  dynamicMask.multiply(pixelArea).reduceRegion({
    reducer:    ee.Reducer.sum(),
    geometry:   AOI,
    scale:      10,
    maxPixels:  1e12,
    bestEffort: true,
    tileScale:  4
  }).evaluate(function(areaResult) {
    if (areaResult) {
      var rawVal = 0;
      for (var k in areaResult) { rawVal = areaResult[k] || 0; break; }
      areaEstLabel.setValue(
        'Seagrass area: ' + (rawVal/1e6).toFixed(2) + ' km2  (mask >= ' + threshold.toFixed(2) + ')'
      );
    }
  });

  mapPanel.layers().reset();
  addS2Layer(yearVal, AOI);
  addMaskLayer(dynamicMask, AOI, threshold);
  updateMaskCheckLabel(threshold);
  addMPALayer();
  addTrainingLayer();
  addAOIOutlineLayer(AOI);
  addAGCLayer(loadAGC(yearVal, AOI, dynamicMask), yearVal, showAGCCheck.getValue());
  updateAGCCheckLabel(yearVal);

  resultsPanel.clear();
  resultsPanel.add(ui.Label('Computing AGC totals...  (mask >= ' + threshold.toFixed(2) + ')', {
    color:'#e65100', fontSize:'11px', margin:'8px 0'
  }));

  // aggregateYear() computes N_pixels and area_m2 internally per year --
  // no shared inflation factor needs to be pre-computed here (DEFF/ICC
  // uses the fixed N_EFFECTIVE constant, not an area-derived ratio).
  var tsFC = ee.FeatureCollection(
    APPLY_YEARS.map(function(year) {
      return aggregateYear(year, loadAGC(year, AOI, dynamicMask), AOI);
    })
  ).sort('year');

  tsFC.evaluate(function(result) {
    if (!result || !result.features) {
      statusLabel.setValue('Evaluation failed - check AGC assets.');
      resultsPanel.clear();
      resultsPanel.add(ui.Label('Failed to load AGC data. Check asset paths.', {
        color:'#c62828', fontSize:'11px'
      }));
      return;
    }
    statusLabel.setValue('Analysis complete - mask >= ' + threshold.toFixed(2));
    renderResults(method, result, tsFC, dynamicMask, AOI, threshold);
  });
}

// ==========================================================
// 16. DRAWING TOOLS - polygon complete
// ==========================================================
mapPanel.drawingTools().setShown(true);
mapPanel.drawingTools().setDrawModes(['polygon']);
mapPanel.drawingTools().setLinked(false);

mapPanel.drawingTools().onDraw(function(geometry) {
  mapPanel.drawingTools().layers().reset();

  ee.Geometry(geometry).area({ maxError: 1000 }).evaluate(function(areaM2) {
    var areaKm2 = areaM2 / 1e6;

    if (areaKm2 > MAX_AOI_KM2) {
      aoiStatusLabel.setValue(
        'AOI ' + areaKm2.toFixed(1) + ' km2 exceeds ' + MAX_AOI_KM2 + ' km2 - please redraw smaller'
      );
      statusLabel.setValue('AOI too large - redraw');
      currentAOI = null;
      mapPanel.drawingTools().layers().reset();
      return;
    }

    currentAOI = ee.Geometry(geometry);
    aoiStatusLabel.setValue('AOI: ' + areaKm2.toFixed(1) + ' km2');
    statusLabel.setValue('Running analysis...');
    runAnalysis(currentAOI);
  });
});

function clearInspectMarker() {
  var layers = mapPanel.layers();
  for (var li = 0; li < layers.length(); li++) {
    if ((layers.get(li).getName() || '') === 'Inspect point') {
      layers.remove(layers.get(li));
      break;
    }
  }
}

// ==========================================================
// 17. CLEAR & REDRAW BUTTON
// ==========================================================
btnClear.onClick(function() {
  currentAOI       = null;
  _pixelInfoLabel  = null;
  _pixelChartPanel = null;
  _dynamicMask     = null;

  mapPanel.drawingTools().layers().reset();
  mapPanel.layers().reset();

  aoiStatusLabel.setValue('AOI: not drawn yet');
  areaEstLabel.setValue('Seagrass area: -- km2');
  statusLabel.setValue('Status: ready - zoom & draw AOI');

  resultsPanel.clear();
  resultsPanel.add(ui.Label('Zoom to target area & draw a polygon AOI', {
    fontSize:'11px', color:'#aaa', fontStyle:'italic', margin:'6px 0'
  }));

  toggleMode(false);
});

// ==========================================================
// 18. SLIDER CALLBACKS
// ==========================================================
var probPending = false;
probSlider.onChange(function(val) {
  probValueLabel.setValue(val.toFixed(2));
  areaEstLabel.setValue('Seagrass area: computing...');
  if (!currentAOI) return;
  statusLabel.setValue('Waiting... ' + val.toFixed(2));
  probPending = true;
  var capturedVal = val;
  ui.util.setTimeout(function() {
    if (probSlider.getValue() !== capturedVal) return;
    probPending = false;
    statusLabel.setValue('Recomputing mask >= ' + capturedVal.toFixed(2) + '...');
    runAnalysis(currentAOI);
  }, 600);
});

var yearPending = false;
yearSlider.onChange(function(val) {
  yearValueLabel.setValue(String(val));
  if (!currentAOI) return;
  statusLabel.setValue('Waiting... ' + val);
  yearPending = true;
  var capturedYear = val;
  ui.util.setTimeout(function() {
    if (yearSlider.getValue() !== capturedYear) return;
    yearPending = false;
    var threshold   = probSlider.getValue() || 0.50;
    var dynamicMask = buildDynamicMask(threshold);
    statusLabel.setValue('Updating basemap -> ' + capturedYear + '...');
    mapPanel.layers().reset();
    addS2Layer(capturedYear, currentAOI);
    addMaskLayer(dynamicMask, currentAOI, threshold);
    updateMaskCheckLabel(threshold);
    addMPALayer();
    addTrainingLayer();
    addAOIOutlineLayer(currentAOI);
    addAGCLayer(loadAGC(capturedYear, currentAOI, dynamicMask), capturedYear, showAGCCheck.getValue());
    updateAGCCheckLabel(capturedYear);
    statusLabel.setValue('S2 p60 - ' + capturedYear);
  }, 600);
});

forecastSelect.onChange(function() {
  if (!currentAOI) return;
  runAnalysis(currentAOI);
});

// ==========================================================
// 19. LAYER CHECKBOX CALLBACKS
// ==========================================================
showMPACheck.onChange(function(checked) {
  mapPanel.layers().forEach(function(layer) {
    if ((layer.getName() || '') === 'MPA Boundaries') layer.setShown(checked);
  });
});

showTrainingCheck.onChange(function(checked) {
  mapPanel.layers().forEach(function(layer) {
    if ((layer.getName() || '') === 'Training points (2017-2023)') layer.setShown(checked);
  });
});

showAOICheck.onChange(function(checked) {
  mapPanel.layers().forEach(function(layer) {
    if ((layer.getName() || '') === 'AOI (user drawn)') layer.setShown(checked);
  });
});

showAGCCheck.onChange(function(checked) {
  var yr = yearSlider.getValue() || 2024;
  mapPanel.layers().forEach(function(layer) {
    if ((layer.getName() || '') === 'AGC - ' + yr) layer.setShown(checked);
  });
});

showMaskCheck.onChange(function(checked) {
  mapPanel.layers().forEach(function(layer) {
    if ((layer.getName() || '').indexOf('Seagrass mask') === 0) layer.setShown(checked);
  });
});

// ==========================================================
// 20. INITIALISE
// ==========================================================
mapPanel.onClick(function(coords) {
  if (!inspectMode) return;
  pixelInspectorClick(coords);
});

toggleMode(false);
