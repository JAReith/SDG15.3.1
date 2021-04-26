//Upload area of interest and relevant Layers
var area = ee.FeatureCollection('users/jonathanreith/KK')
var LC  = ee.Image("users/jonathanreith/LD/LandCover19"),
    LPP = ee.Image("users/jonathanreith/LD/ProductivityPerformance19"), // 0ND/ 1D
    LPS = ee.Image("users/jonathanreith/LD/ProductivityState19"),      //0ND / 1 IMP / 2D
    LPT = ee.Image("users/jonathanreith/LD/ProductivityTrend19"),       //0ND / 1 IMP / 2D
    SOC = ee.Image("users/jonathanreith/SOC18")                 //1 nodata /2 deg / 3stable /4 imp
//Visualize the sub-indicators of Land Productivity and compute it

//Map.addLayer(LPS.clip(area),{max: 2,"palette":['ffffff','green','red']},'STATE')

var LPS = LPS       //imp, stab, deg
              .remap([0,	1,	2 ],
                     [1,	0,	2 ]);
var LPT = LPT
              .remap([0,	1,	2
                      ],
                    [0,	2,	1  ]);

//Map.addLayer(LPP.clip(area),{max: 1,"palette":['ffffff','red'        ]},'PERFORMANCE')
Map.addLayer(LPS.clip(area),{max: 2,"palette":['ffffff','green','red']},'STATE')
Map.addLayer(LPT.clip(area),{max: 2,"palette":['ffffff','green','red']},'TREND')
//First Trend and State
var LPTS = (LPT.multiply(10)).add(LPS)
var LPTS = LPTS
              .remap([0,	1,	2,
                      10,	11,	12,
                      20,	21,	22],
                    [5,	0,	15,
                      10,	10,	8,
                      20,	20,	20
                      ]);
// Then combine it with Performance
var LP = LPTS.add(LPP)
var LP = LP
              .remap([0,	1,
                        5,	6,
                        8,	9,
                        10,	11,
                        15,	16,
                        20,	21,],
                    [1,	1,
                      1,	2,
                      0,	1,
                      0,	0,
                      3,	4,
                      4,	4
                      ]);
//Land Producitvity 0 Improvement, 1 stable, 2 stablebutstressed, 3 early signs, 4 degraded
Map.addLayer(LP.clip(area),{"palette":['green','White','Gold','Orange','red']} , 'Land Productivity');
Map.addLayer(LC.clip(area),{'max': 2, "palette":['PapayaWhip','green','red']} , 'Land Cover');
Map.addLayer(SOC.clip(area),{'min':1,'max': 4, "palette":['PapayaWhip','red','PapayaWhip','green']} , 'Soil Organic Carbon');
// Combine Land PRoducitivty with Land Cover and Soil Organic Carbon. Method: One out, all out
var LD = (LP.multiply(10)).add(LC)
var LD = LD
              .remap([0,	1,	2,
                      10,	11,	12,
                      20,	21,	22,
                      30,	31,	32,
                      40,	41,	42],
                     [0,	0,	2,
                      1,	0,	2,
                      1,	0,	2,
                      1,	0,	2,
                      2,	2,	2

                      ]);
var LD = (LD.multiply(10)).add(SOC)
var LD = LD
              .remap([1,	2,	3,	4,
                      11,	12,	13,	14,
                      21,	22,	23,	24 ],
                      [1,	3,	1,	1,
                      2,	3,	2,	1,
                      3,	3,	3,	3         ]);
Map.addLayer(LD.clip(area),{'max': 3, "palette":['green','PapayaWhip','red']} , 'Land Degradation');


Export.image.toDrive({
  image: LD.clip(area),
  description: "LandDegradation2019",
  scale: 30,
  region: area,
  maxPixels:3000000000
  });
/*  
Export.image.toDrive({
  image: SOC.clip(area),
  description: "SOCDEG2020",
  scale: 30,
  region: area,
  maxPixels:3000000000
  });
  */
