#----------------------------------------------------------------------------------------------------------------------------------# 
# Kalyani Ramanan 30/04/24
#
# Script based on Plume_Search.ipynb, designed to retrieve L2 methane concentration enhancmement maps from Sentinel-2 data. Using 
# GEE Python API. First part calculates fractional TOA reflectance using Sentinel-2 L1C TOA reflectance data bands 11 and 12., u
# sing methods from Varon 2021. Then. these fractional values are converted to absolute values of concentration in mol/m2.
# This is first attempted over a small area, which will then be run in batches in parallel.
#--------------------------------------------------------------------------------------

# Import libraries and initialize GEE
import ee
ee.Authenticate()
ee.Initialize(project = 'ee-krgeophys')
import os
import geemap
from geemap import *
import numpy as np
import xarray as xr # May need to install using %pip install xarray
import datetime
import matplotlib
from matplotlib import pyplot as plt
print('hi')
#=-------------------------------------------------------------------------------------

# Define date and location of Hassi Messaoud plume
source_lat = 31.6585
source_lon = 5.9053
print('Plume source at (', source_lat, ',',source_lon, ')')

# Define dates to filter image collection by (date of plume detection)

plume_date_start = '2019-11-20'
plume_date_end = '2019-11-21'
print('Plume detection date:' , plume_date_start)

# Define geometry as rectangle with coordinates at 2km distance from source (found using google maps)
lonmin = 5.88419
lonmax = 5.92609
latmin = 31.64011
latmax = 31.67640

# Define GEE geometry shape over this area
geometry_area = ee.Geometry.Rectangle([lonmin, latmin, lonmax, latmax], None, False)

# Load TOA collection
collection = (
    ee.ImageCollection("COPERNICUS/S2_HARMONIZED")
    .filterBounds(geometry_area) # Filter by geographical area using geometry rectangle defined above
    .filterDate(plume_date_start, plume_date_end) # Filter by date of plume detection
)

# Get image as first image for filtered collection
img = collection.first()
# Print image ID and spacecraft name
print('Image ID:', img.id().getInfo())
print('Spacecraft name:', img.get('SPACECRAFT_NAME').getInfo())

## Plot Sentinel-2 B12 TOA reflectance
# Plot only B12 to see plume as show in Varon (2021) Figure 2b
mapB12 = geemap.Map()
# Define visualisation parameters
vis = {
    'min': 0.35,
    'max': 0.75,   
    'bands': ['B12']
   
}

mapB12.setCenter(source_lon, source_lat, 13.5) # Define map centre
mapB12.addLayer(img.divide(10000), vis, 'Sentinel-2') # Add bands of Sentinel-2 data
# Map.addLayer(geometry_area, {}, 'geometry_area') # Overlay rectangle used in Varon (2021) figure 2e
# Map.remove_colorbar(clear)
mapB12.add_colorbar(vis, label='TOA Refelctance')

display(mapB12)

# Plot only B11 to see no plume as show in Varon (2021) Figure 2b
mapB11 = geemap.Map()
# Define visualisation parameters
vis = {
    'min': 0.35,
    'max': 0.75,   
    'bands': ['B11']
   
}

mapB11.setCenter(source_lon, source_lat, 13.5) # Define map centre
mapB11.addLayer(img.divide(10000), vis, 'Sentinel-2') # Add bands of Sentinel-2 data
# Map.addLayer(geometry_area, {}, 'geometry_area') # Overlay rectangle used in Varon (2021) figure 2e
# Map.remove_colorbar(clear)
mapB11.add_colorbar(vis, label='TOA Refelctance')

display(mapB11)
#-----------------------------------------------------------------------------------------



# Calculate fractional difference in two bands and plot

# Calculate linear regression of B12 over B11 with intercept forced through 0
# Import image using image name
img = ee.Image('COPERNICUS/S2_HARMONIZED/20191120T101321_20191120T101408_T31SGR')
# Create new image that is the concatenation of three images: a constant, the SWIR1 band and the SWIR2 band
B11 = img.select('B11').divide(10000)
B12 = img.select('B12').divide(10000)
imRegress = ee.Image.cat(B11, B12)

linearRegression = imRegress.reduceRegion(
    reducer = ee.Reducer.linearRegression(numX = 1, numY = 1),
    geometry = geometry_area,
    scale = 20,
)

# Print scaling factor value, should be nearish to 1
print('scaling factor:',linearRegression.get('coefficients').getInfo()[0])

# Calculate fractional difference usign MBMP method

# Use image.expression() to write equation for R_MBSP
R_MBSP = img.expression(
    '(c*R12 - R11)/R11 ', # Where this number comes from least squares difference scale above
    {
        
        'c': linearRegression.get('coefficients').getInfo()[0],
        # 'c': c, 
        'R11': img.select('B11').divide(10000),
        'R12': img.select('B12').divide(10000),
    },
)

# Fractional methane column enhancement
F_MBSP = R_MBSP.subtract(-0.029)

# Get max and min R_MBSP values to compare with Fei's (min - -0.16 and max = -0.0150)
minval = F_MBSP.reduceRegion(ee.Reducer.min(), geometry = geometry_area).get('constant').getInfo()
maxval = F_MBSP.reduceRegion(ee.Reducer.max(), geometry = geometry_area).get('constant').getInfo()
print('min_val:',F_MBSP.reduceRegion(ee.Reducer.min(), geometry = geometry_area).get('constant').getInfo())
print('max_val:',F_MBSP.reduceRegion(ee.Reducer.max(), geometry = geometry_area).get('constant').getInfo())

#  Display F_MBSP on map

Map = geemap.Map()

# Define visualisation parameters for plot
# vis = { 'min': -0.5, 'max': 0.05, "palette": ['FF0000',	'FFFFFF', '0000FF']}
# vis = { 'min': -2, 'max': 2, "palette": ['FF0000',	'FFFFFF', '0000FF']}
vis = { 'min': minval, 'max': maxval, "palette": ['FF0000',	'FFFFFF', '0000FF']}

Map.setCenter(source_lon, source_lat, 13.75) # Define map centre
Map.addLayer(F_MBSP, vis, 'Sentinel-2') # Add bands of Sentinel-2 data
# Map.addLayer(geometry_area, {}, 'geometry_area') # Overlay rectangle used in Varon (2021) figure 2e
Map.add_colorbar(vis, label='Methane enhancement ')
display(Map)

#--------------------------------------------------------------------------------------------------

# Now we will use the radiative transfer model Radtran to convert from these fractional values of methane
# enhancement to absolute values. Code provided by D. Varon, licensed by GHGSat

# BEcause the model takes a long time to run, we will define many small areas which can be run in parallel
# to get this result. First, we are testing this with a small 3x3 pixel area

# Define date and location of Hassi Messaoud plume
source_lat = 31.6585
source_lon = 5.9053
print('Plume source at (', source_lat, ',',source_lon, ')')

# Define dates to filter image collection by (date of plume detection)
plume_date_start = '2019-11-20'
plume_date_end = '2019-11-21'
print('Plume detection date:' , plume_date_start)


#### Make a rectangle over area over which Varon 2021 images the Hassi Messaoud plume

# Define geometry as rectangle 
# Define geometry as rectangle with coordinates at 2km distance from source (found using google maps)
lonmin = source_lon - 0.0002
lonmax = source_lon + 0.0002
latmin = source_lat - 0.0002
latmax = source_lat + 0.0002

# Set up GEE geometry shape
tiny_area = ee.Geometry.Rectangle([lonmin, latmin, lonmax, latmax], None, False)

# Turn list into 2d array and check its shape
F_tiny_area = F_MBSP.sampleRectangle(tiny_area)
F_array = np.array((F_tiny_area.get('constant').getInfo()))
print(F_array.shape)

