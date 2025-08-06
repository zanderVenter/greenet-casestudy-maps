/*******************************
 Grassland Relative Yield Estimation
 using Sentinel-2 imagery (NDVI-based)
 for any Area of Interest (AOI)

 Developed as part of GreeNet project
 Contact: zander.venter@nina.no
********************************/

// ----- USER INPUTS -----

// Import or draw your area of interest as a geometry asset
var AOI = geometry; // Replace 'geometry' with your imported asset if needed

// Define analysis year and growing season months
var startYear = 2021;
var endYear = 2021;
var startMonth = 5;
var endMonth = 9;

// Cloud masking parameters
var cloudFilterThreshold = 30;  // Maximum % cloudiness per scene
var clearScoreThreshold = 0.60; // Cloud Score Plus threshold

// Output projection and scale
var outputCRS = 'EPSG:3035';
var outputScale = 10; // meters

// Output file name
var exportLabel = 'relative_yield_output';

// ----- GRASSLAND MASKING -----

// Load CLC+ 2018 raster (10 m) and mask for herbaceous classes
var clcplus = ee.Image('projects/nina/Europe_misc/CLMS_CLCplus_RASTER_2018_010m_eu_03035_V1_1');
var grass = clcplus.eq(6).or(clcplus.eq(7)); // 6: Permanent herbaceous, 7: Periodically herbaceous

// Exclude moss & lichen areas using ESA WorldCover
var moss = ee.ImageCollection("ESA/WorldCover/v200").mosaic().eq(100); // 100: Moss & Lichen
grass = grass.where(moss, 0).unmask(0); // Set moss to 0, others unmasked

// ----- SENTINEL-2 PROCESSING FUNCTION -----

function getS2(aoi, startYear, endYear, startMonth, endMonth) {
  var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

  var s2 = ee.ImageCollection("COPERNICUS/S2_HARMONIZED")
    .filterBounds(aoi)
    .filter(ee.Filter.calendarRange(startYear, endYear, 'year'))
    .filter(ee.Filter.calendarRange(startMonth, endMonth, 'month'))
    .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', cloudFilterThreshold);

  // Link with cloud score plus and mask cloudy pixels
  s2 = s2.linkCollection(csPlus, ['cs_cdf']).map(function(img) {
    return img.updateMask(img.select('cs_cdf').gte(clearScoreThreshold));
  });

  // Calculate NDVI and return image collection
  s2 = s2.select(['B4', 'B8'], ['red', 'nir'])
         .map(function(img) {
            return img.addBands(img.normalizedDifference(['nir', 'red']).rename('ndvi'));
         });

  return s2;
}

// ----- NDVI COMPOSITE AND NORMALIZATION -----

function processNDVI(aoi) {
  var s2 = getS2(aoi, startYear, endYear, startMonth, endMonth);
  
  var ndviMean = s2.select('ndvi').mean().clip(aoi);
  
  // Apply grassland mask
  ndviMean = ndviMean.updateMask(grass);

  // Compute 2nd and 98th percentiles for normalization
  var perc = ndviMean.reduceRegion({
    reducer: ee.Reducer.percentile([2, 98]),
    geometry: aoi,
    scale: outputScale,
    bestEffort: true
  });

  var min = ee.Number(perc.get('ndvi_p2'));
  var max = ee.Number(perc.get('ndvi_p98'));

  // Clamp and scale to 0–100
  var relYield = ndviMean
    .clamp(min, max)
    .unitScale(min, max)
    .multiply(100)
    .round()
    .unmask(999) // Nodata
    .int();

  return relYield;
}

// ----- RUN AND EXPORT -----

var result = processNDVI(AOI);
Map.centerObject(AOI);
Map.addLayer(result.updateMask(result.neq(999)), {min: 0, max: 100, palette: ['white', 'green']}, 'Relative Yield');

// Optional export to Drive
Export.image.toDrive({
  image: result,
  region: AOI,
  description: exportLabel,
  scale: outputScale,
  crs: outputCRS,
  maxPixels: 1e11
});
