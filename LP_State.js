//import study area
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
  )
);

// Filter collections and prepare them for merging.
var l8 = oliCol.filter(colFilter).map(prepOLI);
var l7 = etmCol.filter(colFilter).map(prepETM);
var l5 = tmCol.filter(colFilter).map(prepETM);
//print(l8,'Landsat 8');
//print(l7,'Landsat 7');
//print(l5,'Landsat 5');

// Merge the collections.
var collMerged = l8
  .merge(l7)
  .merge(l5);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////STATE///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
For each pixel, use the annual mean of NDVI for the baseline period to compute a frequency
distribution. In case the baseline period missed some extreme values in NDVI, add 5% on both
extremes of the distribution. 
*/
//Baseline till 2012
var yearsb     = ee.List.sequence(1999,2011); //because november, one year earlier
var baseline  = ee.ImageCollection.fromImages(
  yearsb.map(function (y) {
  var start = ee.Date.fromYMD(y, 11, 1); // From November on
  var stop = start.advance(7, 'month');
  var sumVeg = collMerged.filterDate(start, stop).select('nd').mean(); // 
return stop, sumVeg.set('year', y);
}));

//Baseline Comparison Period 2013-2015
var yearsbc    = ee.List.sequence(2012,2014); //because november, one year earlier
var bcompare  = ee.ImageCollection.fromImages(
  yearsbc.map(function (y) {
  var start = ee.Date.fromYMD(y, 11, 1); // From November on
  var stop = start.advance(7, 'month');
  var sumVeg = collMerged.filterDate(start, stop).select('nd').mean();
return sumVeg.set('year', y);
}));
//Comparison Period 2016-2019
var yearsc = ee.List.sequence(2015,2018); //because november, one year earlier
var compare = ee.ImageCollection.fromImages(
  yearsc.map(function (y) {
  var start = ee.Date.fromYMD(y, 11, 1); // From November on
  var stop = start.advance(7, 'month');
  var sumVeg = collMerged.filterDate(start, stop).select('nd').mean(); 
return sumVeg.set('year', y);
}));

//THE BASELINE TO COMPARE TO
// add 5% on both extremes of the frequency  distribution
var max     = baseline.reduce(ee.Reducer.max());
var max105  = max.add(max.multiply(0.05));
var min     = baseline.reduce(ee.Reducer.min());
var min05   = min.subtract(max.multiply(0.05));

// Determine the mean NDVI for the Baseline and comparison period
var meanbase  = baseline.mean();
var meanbcomp = bcompare.mean();
var meancomp  = compare.mean();

//Assign to the mean NDVI for the baseline period the number corresponding to that percentile class (1-10)
// USe ceil to assign it to a full integer, instead of float!
//E.g. 7.2 belongs to class 8, as 9.9 belongs to 10 and 0.2 belongs to the first class
var classbase = meanbase.expression(
    '((mean - min) / (max - min) * 10)', {
      'mean': meanbase,
      'min': min05,
      'max': max105
})
.ceil().rename('class');
//print(classbase,'classbase norm')
var palettes = require('users/gena/packages:palettes');
var palette = palettes.cmocean.Speed[7];
//Map.addLayer(classbase.clip(area),{ "min":0,"max":10, palette: palette},'Baseline period');
/*
// Baseline Comparison Min Max
var maxc = bcompare.reduce(ee.Reducer.max());
var minc = bcompare.reduce(ee.Reducer.min());
// Comparison MinMax
var maxc = compare.reduce(ee.Reducer.max());
var minc = compare.reduce(ee.Reducer.min());
*/
var classbcomp = meanbcomp.expression(
    '((mean - min) / (max - min) * 10)', {
      'mean': meanbcomp,
      'min': min05,    
      'max': max105
    })
    .ceil().rename('class');

var classcomp = meancomp.expression(
    '((mean - min) / (max - min) * 10)', {
      'mean': meancomp,
      'min': min05,    
      'max': max105
    })
    .ceil().rename('class');

//print(classcomp,'classbase norm')
//Map.addLayer(classbcomp.clip(area),{ "min":0,"max":10, palette: palette},'Baseline Comparison Period 13-15');
//Map.addLayer(classcomp.clip(area),{ "min":0,"max":10, palette: palette},'Comparison Period 16-19');

///////////////////Histrogram
var joined = ee.Image.cat([classbase, classbcomp,  classcomp]);
//print(joined,'joined')
var joined = joined.select(
    ['class', 'class_1', 'class_2'], // old names
    ['Baseline', 'Baseline Comparison', 'Comparison']  // new names
    );
var options = {
  title: 'Productivity State Histogram',
  fontSize: 20,
  hAxis: {title: 'Productivity Class'},
  vAxis: {title: 'Count of Pixels'},
  series: {
    0: {color: 'blue'},
    1: {color: 'green'},
    2: {color: 'red'},
    }};
var histogram = ui.Chart.image.histogram(joined, area, 1000)
          .setSeriesNames(['2000-2012', '2013-2015',['2016-2019']])
          .setOptions(options);
//print(histogram);

////////////Computation  of the Baseline Degradation classes /////////////////////7
var Bstate = classbase.subtract(classbcomp).rename('degradationstate')
//print(Bstate)
//Values equal or higher than 2 are degradation (7-4=3)
//Values lower than  -2 on the contrary to be seen as improvement
var deg_stab_imp = Bstate.gt(-2).add(Bstate.gte(2));
//0= better
//1= stable
//2=degradation
//Map.addLayer(deg_stab_imp.clip(area),{"max":2,"palette":['02ff0a','ffffff','ff0000']},'Baseline Degradation and Improvement');
////////////Computation  of the Baseline Degradation classes /////////////////////7
var Cstate = classbase.subtract(classcomp).rename('degradationstate')
//print(Cstate)
//Values equal or higher than 2 are degradation (7-4=3)
//Values lower than are -2 the contrary to be seen as improvement
var deg_stab_imp = Cstate.gt(-2).add(Cstate.gte(2));
//0= better
//1= stable
//2=degradation
//Map.addLayer(deg_stab_imp.clip(area),{"max":2,"palette":['02ff0a','ffffff','ff0000']},'Comparison Degradation and Improvement');


///////////////////////////////////////////////////////////////
///////////////// UNCCD Reporting /////////////////////////////
//Computation  of the Baseline Degradation classes
//Values equal or higher than 2 are degradation (7-4=3)
var BLdegradation = Bstate.gt(-2).add(Bstate.gte(2));
Map.addLayer(BLdegradation.clip(area),{"max":2,"palette":['02ff0a','ffffff','ff0000']},'BL DEGRAD');
//0= better ....below minus 2 Or rather the rest
//1= stable  ...All the values between -2 and 2
//2=degradation ... All the values greater than 2
var CPdegra       = Cstate.gt(-2).add(Cstate.gte(2));
Map.addLayer(CPdegra.clip(area),{"max":2,"palette":['02ff0a','ffffff','ff0000']},'Comparison Degradation');
//0= better ....below minus 2 Or rather the rest
//1= stable  ...All the values between -2 and 2
//2=degradation ... All the values greater than 2
/// Degradation for 2019
var BLdegradation_  = BLdegradation.multiply(10); 
var Degradation19 = BLdegradation_.add(CPdegra);
var Degradation19 = Degradation19
              .remap([0,	1,	2,
                      10,	11,	12,
                      20,	21,	22
                      ],
                     [1,	1,	2,
                      1,	0,	2,
                      0,	2,	2
                      ]);
                     //0 Not Degraded
                     //2 Degraded= Was degraded and stayed degraded OR bacame degraded in 2019!! 1 Improvement
Map.addLayer(Degradation19.clip(area),{"palette":['ffffff','green','red']},'Degradation 2019')
///////////////////////////EXPORT/////////////////////////////////////////////////////////////
/*
Export.image.toAsset({
  image: BLdegradation.clip(area),
  description: "ProductivityState15",
  scale: 30,
  region: area,
  pyramidingPolicy: {
    'default': 'sample'
  }
  });

Export.image.toDrive({
  image: Degradation19.clip(area),
  description: "Productivity State",
  scale: 30,
  region: area,
  maxPixels:3000000000
  });
/*
Export.image.toDrive({
  image: classbase.clip(area),
  description: "Productivity State Baseline",
  scale: 30,
  region: area,
  maxPixels:3000000000
  });
  
Export.image.toDrive({
  image: classcomp.clip(area),
  description: "Productivity State Comparison",
  scale: 30,
  region: area,
  maxPixels:3000000000
  });
/*
////////////////////////////////////////////////////////////////////////////////////////////////////////7
var chart = ui.Chart.image.series({
  imageCollection: collMerged,
  region: area,
  reducer: ee.Reducer.mean(),
  scale: 30
}).setOptions({title: 'NDVI over time',trendlines: {0: {
        color: 'CC0000'
      }},
      lineWidth: 1,
      pointSize: 1,
    });
print(chart);
//////////////////////////////////////
var yearchart = ui.Chart.image.series(byYear, area,ee.Reducer.mean(),30,'year')
    .setChartType('ScatterChart')
    .setOptions({
      title: 'yearly NDVI mean time series',
      trendlines: {0: {
        color: 'CC0000'
      }},
      lineWidth: 1,
      pointSize: 3,
    });
print(yearchart);
*/
