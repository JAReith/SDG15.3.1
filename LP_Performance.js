var area = ee.FeatureCollection('users/jonathanreith/KK')
var empty = ee.Image().byte();
var outline = empty.paint({
  featureCollection: area,
  color: 1,
  width: 3
});
Map.addLayer(outline, {palette: '000000'}, 'Kiteto & Kongwa');
Map.centerObject(area,9);

// ################################################################
// ### FUNCTIONS ###
// ################################################################
// Define coefficients supplied by Roy et al. (2016) for translating ETM+
// surface reflectance to OLI surface reflectance.
var coefficients = {
  itcps: ee.Image.constant([0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172]).multiply(10000),
  slopes: ee.Image.constant([0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071])
};
// Define function to get and rename bands of interest from OLI.
function renameOLI(img) {
  return img.select(
		['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'pixel_qa'],
		['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'pixel_qa']
	);
}
// Define function to get and rename bands of interest from ETM+.
function renameETM(img) {
  return img.select(
		['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'pixel_qa'],
		['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'pixel_qa']
  );
}
// Define function to apply harmonization transformation.
function etm2oli(img) {
  return img.select(['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'])
    .multiply(coefficients.slopes)
    .add(coefficients.itcps)
    .round()
    .toShort()
    .addBands(img.select('pixel_qa')
  );
}
// Define function to mask out clouds and cloud shadows.
function fmask(img) {
  var cloudShadowBitMask = 1 << 3;
  var cloudsBitMask = 1 << 5;
  var qa = img.select('pixel_qa');
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
    .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return img.updateMask(mask);
}
// Define function to calculate NDVI
function calcnd(img) {
  return img.normalizedDifference(['NIR', 'Red']).rename('nd');
}
// Define function to prepare OLI images.
function prepOLI(img) {
  var orig = img;
  img = renameOLI(img);
  img = fmask(img);
  img = calcnd(img);
  return ee.Image(img.copyProperties(orig, orig.propertyNames()));
}
// Define function to prepare ETM+ images.
function prepETM(img) {
  var orig = img;
  img = renameETM(img);
  img = fmask(img);
  img = etm2oli(img);
  img = calcnd(img);
  return ee.Image(img.copyProperties(orig, orig.propertyNames()));
}
// ################################################################
// ### APPLICATION ###
// ################################################################

// Get Landsat surface reflectance collections for OLI, ETM+ and TM sensors.
var oliCol = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
                    .filterDate('1999-11-01', '2019-10-30')
                    .filter(ee.Filter.calendarRange(11,6,'month')) ;
var etmCol = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
                    .filterDate('1999-11-01', '2019-10-30')
                    .filter(ee.Filter.calendarRange(11,6,'month')) ;
var tmCol  = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
                    .filterDate('1999-11-01', '2015-10-30')
                    .filter(ee.Filter.calendarRange(11,6,'month')) ;

// Define a collection filter.
var colFilter = ee.Filter.and(
  ee.Filter.bounds(area),
ee.Filter.lt('CLOUD_COVER', 80),         
  ee.Filter.lt('GEOMETRIC_RMSE_MODEL', 10),
ee.Filter.or(
    ee.Filter.eq('IMAGE_QUALITY', 9), //9 = best quality
    ee.Filter.eq('IMAGE_QUALITY_OLI', 9)
  ));

// Filter collections and prepare them for merging.
var l8 = oliCol.filter(colFilter).map(prepOLI);
var l7 = etmCol.filter(colFilter).map(prepETM);
var l5 = tmCol.filter(colFilter).map(prepETM);
print(l8,'Landsat 8');
print(l7,'Landsat 7');
print(l5,'Landsat 5');

// Merge the collections.
var collMerged = l8
  .merge(l7)
  .merge(l5);
var Baseline = collMerged.filterDate('1999-11-01', '2015-10-30')
                    .filter(ee.Filter.calendarRange(11,6,'month')) ;
var Comparison = collMerged.filterDate('2015-11-01', '2019-10-30')
                    .filter(ee.Filter.calendarRange(11,6,'month')) ;
///////////////////////////////////////////////////////////////////////////////////////////
//Define the analysis period, and use the time series of NDVI to compute mean the NDVI for each pixel.
var nd_meanb = Baseline.mean()
var nd_meanc = Comparison.mean()
//print(nd_mean,'mean')

/*
For each unit, extract  the mean NDVI values computed in step 1, and create a frequency
distribution. From this distribution determine the value which represents the 90th percentile (we
don’t recommend using the absolute maximum NDVI value to avoid possible errors due to the
presence of outliers). The value representing the 90th percentile will be considered the maximum
productivity for that unit.*/
//Prepare the LC
var image18 = ee.Image('users/jonathanreith/RCMDR_18')
var buildColorPalette = function(colors) {
  var c = []
  for (var k in colors) {
    c.push(colors[k])
  }
  return c.join(',')
}
var rcmrdLandCovers = {
  'NoData'      : 'gray',
  'Forestland'  : 'darkgreen',
  'Grassland'   : 'lime',
  'Cropland'    : 'yellow',
  'Wetland'     : 'blue',
  'Settlement'  : 'red',
  'Otherland'   : 'orange',
  'Cloud'      : 'black',
  'Shadow'      : 'lightgray'
}
var paletteS1 = buildColorPalette(rcmrdLandCovers)
var visS1 = {min: 0, max: 8, palette: paletteS1}

// Weighted smoothing using a 3x3 window euclidean distance weighting from corners
// create a 3x3 kernel with the weights
// kernel W and H must equal weight H and H
// apply mode on neightborhood using weights
// and force operations to be done at native scale
var weights = [[1,2,1],
               [2,3,2],
               [1,2,1]];
var kernel = ee.Kernel.fixed(3,3,weights);
var SCALE = 30;

var lc = image18.reduceNeighborhood({
  reducer: ee.Reducer.mode(),
  kernel: kernel
}).reproject('EPSG:4326', null, SCALE);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////LOAD DATA//////////////////////////////////////////////////////////
//Datasets for ##land capability units## (Wessels 2008)
//Define similar ecologically similar units as the unique intersection of land cover and soil type.
var soil    = ee.Image("OpenLandMap/SOL/SOL_GRTGROUP_USDA-SOILTAX_C/v01")
              .select(['grtgroup'],['soil']);

//Create ecounits
var lc  = lc.multiply(1000); 
var ecounit = soil.add(lc);

//Map.addLayer(soil.clip(area), getVisParams, 'soil');
//Map.addLayer(ecounit.clip(area),{"min":1019,"max":8430}, 'ecounit');

//merge the NDVI with the ecounit. 
var econdvib = ecounit.add(nd_meanb);
var econdvib = ee.Image.cat([ecounit, nd_meanb]);

var econdvic = ecounit.add(nd_meanc);
var econdvic = ee.Image.cat([ecounit, nd_meanc]);
//print(econdvi,'ecounit with NDVI')

//Convert the ecounit raster to polygon
var eeList = ee.List([90]);
var eco90b = econdvib.reduceToVectors({
  geometry: area,
  geometryType: 'polygon',
  eightConnected: false,
  labelProperty: 'ecounit',
  reducer: ee.Reducer.percentile(eeList)
});
//print(eco90,'ecounit 90th percentile')
var eco90c = econdvic.reduceToVectors({
  geometry: area,
  geometryType: 'polygon',
  eightConnected: false,
  labelProperty: 'ecounit',
  reducer: ee.Reducer.percentile(eeList)
});
var ecomaxb = eco90b
  .filter(ee.Filter.notNull(['p90']))
  .reduceToImage({
    properties: ['p90'],
    reducer: ee.Reducer.max()
});
var ecomaxc = eco90c
  .filter(ee.Filter.notNull(['p90']))
  .reduceToImage({
    properties: ['p90'],
    reducer: ee.Reducer.max()
});
//print(ecomax,'ecounit 90th percentile RASTER')

//Compute the ratio of mean NDVI and maximum productivity. Values between 0-1.
var performanceb = nd_meanb.divide(ecomaxb)
var performancec = nd_meanc.divide(ecomaxc)
//print(performance)
var palettes = require('users/gena/packages:palettes');
var palette = palettes.colorbrewer.RdYlGn[11];
Map.addLayer(performanceb.clip(area),{palette: palette}, 'Performance Baseline');
Map.addLayer(performancec.clip(area),{palette: palette}, 'Performance Comparison');
//Find the values that are degraded--- That are less than half of it...so less than 0.5
var degradedb = performanceb.lt(0.5);
var degradedc = performancec.lt(0.5);

Map.addLayer(degradedb.clip(area),{"palette":['#ffffff',"ff0000"]}, 'Performance Degradation Baseline');
Map.addLayer(degradedc.clip(area),{"palette":['#ffffff',"ff0000"]}, 'Performance Degradation Comparison');
///////////////////////////////////////////////////////////////
///////////////// UNCCD Reporting /////////////////////////////
///////////////////////////////////////////////////////////////
// 1 and 2 are degraded. 2 indicates degradation in both years, 1 just in one of there
var performance = degradedb.add(degradedc)
Map.addLayer(performance.clip(area),{"palette":['#ffffff',"ff0000","ff0000","ff0000"]}, 'Performance Degradation 2019');


var options = {
  title: 'Productivity Performance Histogram',
  fontSize: 20,
  hAxis: {title: 'Performance Ratio'},
  vAxis: {title: 'Count of Pixels'},
  }
  
var histogram = ui.Chart.image.histogram(performanceb, area, 1000)
                   .setSeriesNames(['Ratio of observed and maximum NDVI in the same ecounit'])
                   .setOptions(options);
                   
print(histogram);
  /*
Export.image.toDrive({
  image: performance.clip(area),
  description: "ProductivityPerformance",
  scale: 30,
  region: area,
  maxPixels:3000000000
  });
  
  Export.image.toDrive({
  image: degradedb.clip(area),
  description: "ProductivityPerformanceB",
  scale: 30,
  region: area,
  maxPixels:3000000000
  });
  
  Export.image.toDrive({
  image: degradedc.clip(area),
  description: "ProductivityPerformanceC",
  scale: 30,
  region: area,
  maxPixels:3000000000
  });

Export.image.toAsset({
  image: performance.clip(area),
  description: "ProductivityPerformance",
  scale: 30,
  region: area,
  pyramidingPolicy: {
    'default': 'sample'
  }
  });
  */
  
