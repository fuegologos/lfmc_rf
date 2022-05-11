/*
LFMC - A remote sensing approach to predict LFMC at large scale
Extract the monthly average value and standard deviation of NDVI from Landsat 5 images
for each site from the included pixels to an extent equivalent to the MODIS resolution.
Code adapted from https://doi.org/10.1038/s41597-019-0164-9
*/

// (0) Load plots and bounding box
var roi = ee.FeatureCollection("users/acunill86/roi_med");
var plots = ee.FeatureCollection("users/acunill86/plots_med");

// (1) Function to mask clouds, clouds shadows and snow
var maskCloudsShadowsSnow = function(image) {
  var bitCloudsShadow = (1 << 3);
  var bitSnow = (1 << 4);
  var bitClouds = (1 << 5);
  var pixel_qa = image.select('pixel_qa');
  var mask = pixel_qa.bitwiseAnd(bitCloudsShadow).eq(0)
                     .and(pixel_qa.bitwiseAnd(bitSnow).eq(0))
                     .and(pixel_qa.bitwiseAnd(bitClouds).eq(0));
  return image.updateMask(mask);
};

// (2) Function to add NDVI band to Landsat 5 image
var addNDVI = function(image) {
  var ndvi = image.normalizedDifference(['B4', 'B3']).rename('NDVI');
  return image.addBands(ndvi);
};

// (3) Add NDVI band to Landsat image collection and mask clouds, clouds shadows and snow
var landsatCollectionWithNDVI = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
  .filterBounds(roi.geometry())
  .map(addNDVI)
  .map(maskCloudsShadowsSnow);

// (4) Function to add features to "plots" and to export to Google Drive as a "CSV" file (one file for each month and year)
var addFeatures = function(n) {
  // Get one object from the list by index 'n' and convert to GEE image
  var image = imagesByMonth.get(n);
  image = ee.Image(image);
  // Compute and add features to plots
  var finalTable = plots.map(function(feature) {
    var buffer500 = feature.buffer(232).bounds();  // create a square buffer matching MODIS resolution (~464 m)
    return feature.set('ndviSD',
                        image.reduceRegion({geometry: buffer500.geometry(),
                                            reducer: ee.Reducer.stdDev(),
                                            scale: 30,
                                            })
                                            .get('NDVI'))
                  .set('ndviMean',
                        image.reduceRegion({geometry: buffer500.geometry(),
                                            reducer: ee.Reducer.mean(),
                                            scale: 30,
                                            }).get('NDVI'))
                  .set('unmaskedPixels',
                        image.reduceRegion({geometry: buffer500.geometry(),
                                            reducer: ee.Reducer.count(),
                                            scale: 30,
                                            }).get('NDVI'))
                  .set('totPixels',
                        image.unmask()
                             .reduceRegion({geometry: buffer500.geometry(),
                                            reducer: ee.Reducer.count(),
                                            scale: 30,
                                            }).get('NDVI'));
  });
  return Export.table.toDrive({collection:finalTable,
                              description: 'L5_month_' + (n + monthStart).toString() + '_year_' + year.toString(),
                              fileFormat: 'CSV'
                              });
  };

// (5) Function to create a list of image collections of monthly NDVI inside year range
var createMonthlyNDVI = function(m){
  m = ee.Number(m);
  var yr = ee.Number(year);
  var startDate = ee.Date.fromYMD(yr, m, 1);
  var endDate = startDate.advance(1, 'month');
  var monthMeanNDVI = landsatCollectionWithNDVI.select('NDVI')
    .filterDate(startDate, endDate)
    .mean();
  return monthMeanNDVI;
};

// (6) Loop through target dates: 01-02-2000 to 31-10-2011
for (var year = 2000; year <= 2011; year += 1) {
  // Define sequence of month to be included
  if (year == 2000) {
    var monthStart = 2;
    var monthEnd = 12;
  } else if (year == 2011) {
    var monthStart = 1;
    var monthEnd = 10;
  } else {
    var monthStart = 1;
    var monthEnd = 12;
  }
  
  // Create list of image collections containing NDVI image of each month inside years range
  var months = ee.List.sequence(monthStart, monthEnd);
  var imagesByMonth = months.map(createMonthlyNDVI).flatten();
  
  // Create list of integers corresponding to the months
  var imageIndicesArray = [];
  var endIndex = monthEnd - monthStart;
  for (var i = 0; i <= endIndex; i++) {imageIndicesArray.push(i)}

  // Make the final computation
  imageIndicesArray.map(addFeatures);  
}

// END ---
