#----------------------------------------------------------------------------------------------------------------------------------# 
# Kalyani Ramanan 01/05/2024
#
# Function based on first part of Plume_Search.ipynb, designed to return F_MBSP using Varon 2021 method for input date and source coordinates.
# Using GEE Python API
#--------------------------------------------------------------------------------------


def retrieve_F_MBSP(source_lon, source_lat, image_date, image_next_day):

    # Variables:
    # source_lon: coordinates of source longitude
    # source_lat: coordinates of source latitude
    # image_date: date of image to use
    #
    print('Plume source at (', source_lat, ',',source_lon, ')')

# Make 10x10km roi around the source to define the scene !WARNING! THIS IS SET UP FOR NORTHERN HEMISPHERE

    # Calculate coordinates in meters and back to coordinates to make the bounding box
    source_utm = utm.from_latlon(source_lon, source_lat)
    source_lon_utm = source_utm[0]
    source_lat_utm = source_utm[1]

    minlon = utm.to_latlon(source_lon_utm - 10000)
    maxlon = utm.to_latlon(source_lon_utm + 10000)
    minlat = utm.to_latlon(source_lat_utm - 10000)
    maxlat = utm.to_latlon(source_lat_utm + 10000)


    print(minlon, minlat)








#utm.from_latlon(lon, lat)




# BELOW IS ONLY FOR REFERENCE - TO BE DELETED 







# Import libraries and initialize GEE - WILL I NEED THIS IF THEY ARE IMPORTED IN THE ORIGINAL BATCHES SCRIPT?
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


# PLOTTING B12 and B11 TOA, COMMENTED OUT=====================================================================
# ## Plot Sentinel-2 B12 TOA reflectance
# # Plot only B12 to see plume as show in Varon (2021) Figure 2b
# mapB12 = geemap.Map()
# # Define visualisation parameters
# vis = {
#     'min': 0.35,
#     'max': 0.75,   
#     'bands': ['B12']
   
# }

# mapB12.setCenter(source_lon, source_lat, 13.5) # Define map centre
# mapB12.addLayer(img.divide(10000), vis, 'Sentinel-2') # Add bands of Sentinel-2 data
# # Map.addLayer(geometry_area, {}, 'geometry_area') # Overlay rectangle used in Varon (2021) figure 2e
# # Map.remove_colorbar(clear)
# mapB12.add_colorbar(vis, label='TOA Refelctance')

# display(mapB12)

# # Plot only B11 to see no plume as show in Varon (2021) Figure 2b
# mapB11 = geemap.Map()
# # Define visualisation parameters
# vis = {
#     'min': 0.35,
#     'max': 0.75,   
#     'bands': ['B11']
   
# }

# mapB11.setCenter(source_lon, source_lat, 13.5) # Define map centre
# mapB11.addLayer(img.divide(10000), vis, 'Sentinel-2') # Add bands of Sentinel-2 data
# # Map.addLayer(geometry_area, {}, 'geometry_area') # Overlay rectangle used in Varon (2021) figure 2e
# # Map.remove_colorbar(clear)
# mapB11.add_colorbar(vis, label='TOA Refelctance')

# display(mapB11)
#==================================================================================================================
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