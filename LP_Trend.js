//import study area
var area = ee.FeatureCollection('users/jonathanreith/KK')
// Set study area as map center.
Map.centerObject(area,9);
//Map.setCenter (36.691, -5.6174, 12)
// Create an empty image into which to paint the features, cast to byte. 
var empty = ee.Image().byte();
// Paint all the polygon edges with the same number and width, display.
var outline = empty.paint({
  featureCollection: area,
  color: 1,
  width: 3
});
Map.addLayer(outline, {palette: '000000'}, 'Kiteto & Kongwa');

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
                    .filterDate('1999-11-01', '2019-10-30')
                    .filter(ee.Filter.calendarRange(11,6,'month')) ;

// Define a collection filter.
var colFilter = ee.Filter.and(
  ee.Filter.bounds(area),
ee.Filter.lt('CLOUD_COVER', 80),       
  ee.Filter.lt('GEOMETRIC_RMSE_MODEL', 10),
ee.Filter.or(
    ee.Filter.eq('IMAGE_QUALITY', 9),
    ee.Filter.eq('IMAGE_QUALITY_OLI', 9)
  ));

// Filter collections and prepare them for merging.
var l8 = oliCol.filter(colFilter).map(prepOLI);
var l7 = etmCol.filter(colFilter).map(prepETM);
var l5 = tmCol.filter(colFilter).map(prepETM);
//print(l8,'Landsat 8');
//print(l7,'Landsat 7');
//print(l5,'Landsat 5');

// Merge the collections.
var collMerged = l8.merge(l7).merge(l5);
/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////PRECIPITATION//////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
var dataset = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
              .filter(ee.Filter
              .date('1999-11-01', '2019-10-30')) 
              .filter(ee.Filter.calendarRange(11,7,'month')) ;
var precipitation = dataset.select('precipitation');
var palettes = require('users/gena/packages:palettes');
var palette = palettes.colorbrewer.YlGnBu[7];
var precipitationVis = {
  max: 3.1257293224334717,
  min: 1.9586076736450195,
  palette: palette,
};
var prec = precipitation.reduce(ee.Reducer.mean());
//Map.addLayer(prec.clip(area), precipitationVis, 'Precipitation');
///////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
//MINMAXNORM

var byears = ee.List.sequence(1999,2014); //always one year earlier because rainy season starts earlier
var cyears = ee.List.sequence(2010,2018); //eight year long Trend minimum
print
var Baseline = ee.ImageCollection.fromImages(
  byears.map(function (y) {
  var start   = ee.Date.fromYMD(y, 11, 1); // From November on
  var stop    = start.advance(7, 'month'); // till june
  var sumVeg  = collMerged.filterDate(start, stop).select('nd').mean() 
  var sumprec = precipitation.filterDate(start, stop).select('precipitation').sum()
  var RUE     = sumVeg.divide(sumprec).rename('RUE')
return ee.Image.cat(sumVeg, sumprec, RUE).set('year', y);
}));

//print(Baseline,'BL RUE')
var Comparison = ee.ImageCollection.fromImages(
  cyears.map(function (y) {
  var start   = ee.Date.fromYMD(y, 11, 1); // From November on
  var stop    = start.advance(7, 'month'); // till june
  var sumVeg  = collMerged.filterDate(start, stop).select('nd').mean()
  var sumprec = precipitation.filterDate(start, stop).select('precipitation').sum()
  var RUE     = sumVeg.divide(sumprec).rename('RUE')
return ee.Image.cat(sumVeg, sumprec, RUE).set('year', y);
}));
//print(Comparison,'CP RUE')
////////////////////////////////////////////
var chart = ui.Chart.image.series({
  imageCollection: Baseline.select('RUE'),
  region: area,
  reducer: ee.Reducer.mean(),
  scale: 1000,
  xProperty: 'year'
}).setChartType('ScatterChart')
.setOptions({
  title: 'RUE curve over Baseline',
  hAxis: {title: 'Date'},
  vAxis: {title: 'RUE'},
 // lineWidth: 1,
//  pointSize: 6,
//dataOpacity: 0.8
  trendlines: {
        0: {
          type: 'linear',
          color: 'lightblue',
          lineWidth: 3,
          opacity: 0.7,
          showR2: true,
          visibleInLegend: true
        }}
})
//.setSeriesNames(['Landsat 5', 'Landsat 7'])
;
//print(chart,'chart');

var chart = ui.Chart.image.series({
  imageCollection: Comparison.select('RUE'),
  region: area,
  reducer: ee.Reducer.mean(),
  scale: 1000,
  xProperty: 'year'
}).setChartType('ScatterChart')
.setOptions({
  title: 'RUE curve over Comparison Period',
  hAxis: {title: 'Date'},
  vAxis: {title: 'RUE'},
 // lineWidth: 1,
//  pointSize: 6,
//dataOpacity: 0.8
  trendlines: {
        0: {
          type: 'linear',
          color: 'lightblue',
          lineWidth: 3,
          opacity: 0.7,
          showR2: true,
          visibleInLegend: true
        }}
})
;
//print(chart,'chart');
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
//Creates a reducer that computes the Kendall's Tau-b rank correlation and p-value on a two-sided test of H0:
//x and y are independent. A positive tau value indicates an increasing trend; 
//negative value indicates a decreasing trend
var RUE_kendall = Baseline.select('RUE').reduce(ee.Reducer.kendallsCorrelation()).clip(area);
//print(RUE_kendall,'RUE_kendall')

var Nb = ee.Image(Baseline.select('RUE').count().clip(area));
//print(Nb,"Nb")
// get the z scode  
//Positive Z scores indicate a trend of increasing productivity
var taub = RUE_kendall.select('RUE_tau');
var num = ee.Image(3).multiply(taub.multiply((Nb.multiply( Nb.subtract(ee.Image(1)))).sqrt()));
var den = (ee.Image(2).multiply((ee.Image(2).multiply(Nb)).add(ee.Image(5)))).sqrt();
var zb = (num.divide(den));
//print(zb,"zb")

//z value lower than -1,96 or higher than 1,96
var bimprovement = zb.gt(1.96);
var bdegradation = zb.lt(-1.96);
var bmask = bdegradation.add(bimprovement);
//Just the relevant change remains!
var bmaskedImage = taub.updateMask(bmask);
            var palettes = require('users/gena/packages:palettes');
            var palette = palettes.colorbrewer.RdYlGn[11];
//Map.addLayer(zb,{palette: palette},'zb value'); //All values
//Map.addLayer(bmaskedImage.clip(area),{ min: -0.6, max: 0.6,'palette': 'FF0000,FFFF00,00FF00'},'Baseline Trend sign');
//just significant values

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
var RUE_kendall = Comparison.select('RUE').reduce(ee.Reducer.kendallsCorrelation()).clip(area);
print(RUE_kendall,'RUE_kendall')
//Number of images per pixel--->should be 16
var Nc = ee.Image(Comparison.select('RUE').count().clip(area));
// get the z scode  
//Positive Z scores indicate a trend of increasing productivity
var tauc = RUE_kendall.select('RUE_tau');
var num = ee.Image(3).multiply(tauc.multiply((Nc.multiply( Nc.subtract(ee.Image(1)))).sqrt()));
var den = (ee.Image(2).multiply((ee.Image(2).multiply(Nc)).add(ee.Image(5)))).sqrt();
var zc = (num.divide(den));
print(zc,"zc")

//z value lower than -1,96 or higher than 1,96
var cimprovement = zc.gt(1.96);
var cdegradation = zc.lt(-1.96);
var cmask = cdegradation.add(cimprovement);
//Just the relevant change remains!
var cmaskedImage = tauc.updateMask(cmask);
            var palettes = require('users/gena/packages:palettes');
            var palette = palettes.colorbrewer.RdYlGn[11];
//Map.addLayer(zc,{palette: palette},'z value'); //All values
//Map.addLayer(cmaskedImage.clip(area),{ min: -0.6, max: 0.6,'palette': 'FF0000,FFFF00,00FF00'},'Comparison Trend sign');
//just significant values

///////////////////////////////////////////////////////////////
///////////////// UNCCD Reporting /////////////////////////////
///////////////////////////////////////////////////////////////
var bimprovement  = bimprovement.remap([0,1],[0,2])   //So that improvement is 2 and will not confused with degradation
var bdegr   = bimprovement.add(bdegradation);         // 0 is ND, 1 =D, 2= IMP
Map.addLayer(bdegr,{ max: 2 ,palette: ['white','red','green']},'Baseline'); //All values

var cimprovement  = cimprovement.remap([0,1],[0,2])
var CPdegr  = cimprovement.add(cdegradation); //0=stable, 1 = degradation, 2 =improvement
Map.addLayer(CPdegr,{ max: 2 ,palette: ['white','red','green']},'Comparison'); //All values

var Degradation19 = bdegr.multiply(10).add(CPdegr);
var Degradation19 = Degradation19
              .remap([0,	1,	2,
                      10,	11,	12,
                      20,	21,	22
                      ],
                     [0,	2,	1,
                      2,	2,	0,
                      1,	2,	1
                      ]); 
                     //0 Not Degraded
                     //2 Degraded= Was degraded and stayed degraded OR became degraded in 2019!! 1 Improvment
Map.addLayer(Degradation19.clip(area),{max: 2,"palette":['white','green','red']},'Degradation binär 2019')



 /* 
Export.image.toAsset({
  image: Degradation19.clip(area),
  description: "ProductivityTrendSDG",
  scale: 30,
  region: area,
  pyramidingPolicy: {
    'default': 'sample'
  }
  });


Export.image.toDrive({
  image: Degradation19.clip(area),
  description: "Productivity TREND SDG",
  scale: 30,
  region: area,
  maxPixels:3000000000
  });
  
Export.image.toDrive({
  image: bdegr.clip(area),
  description: "Productivity TREND 15",
  scale: 30,
  region: area,
  maxPixels:3000000000
  });
  
  Export.image.toDrive({
  image: CPdegr.clip(area),
  description: "Productivity TREND19",
  scale: 30,
  region: area,
  maxPixels:3000000000
  });
  
  */
