#----------------------------------------------------------------------------------------------------------------------------------# 
# Kalyani Ramanan 01/05/2024
#
# Function based on first part of Plume_Search.ipynb, designed to return F_MBSP using Varon 2021 method for input date and source coordinates.
# Using GEE Python API

# Inputs: source coordinates, plume date, next day
#--------------------------------------------------------------------------------------
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
import utm

## Hassi Messaoud plume source
source_lat = 31.6585
source_lon = 5.9053

plume_date = '2019-11-20'
plume_next_day = '2019-11-21'



def retrieve_F_MBSP(source_lon, source_lat, image_date, image_next_day):

    ## Get coordinates of 10km surrounding area

    ##comm out for now
    # (s_lon, s_lat, zone, northing) = utm.from_latlon(source_lon,source_lat)
    # maxlon, maxlat = utm.to_latlon(s_lon +5000, s_lat+5000, zone ,northing)
    # minlon, minlat = utm.to_latlon(s_lon -5000, s_lat-5000, zone ,northing)
    # surrounding10km = ee.Geometry.Rectangle([minlon, minlat, maxlon, maxlat], None, False)

    # Create square geometry area centered on source_lon, source_lat with 10km diameter
    size_add = 5000  # 10km in meters
    minlon = source_lon - size_add / 111000
    maxlon = source_lon + size_add / 111000
    minlat = source_lat - size_add / (111000 * np.cos(np.deg2rad(source_lat)))
    maxlat = source_lat + size_add / (111000 * np.cos(np.deg2rad(source_lat)))
    surrounding10km = ee.Geometry.Rectangle([minlon, minlat, maxlon, maxlat], None, False)

    
    ## Get image from image collection and clip to 10x10km source surrounding area
    Collection = (
        ee.ImageCollection("COPERNICUS/S2_HARMONIZED")
        .filterBounds(surrounding10km)
        .filterDate(image_date, image_next_day) # filter to only one image by date)
        )
    first = Collection.first()
    img = first.clip(surrounding10km)
    
    # Select SWIR bands (S2)
    B11 = img.select('B11').divide(10000)
    B12 = img.select('B12').divide(10000)

    ## Calculate linear regression of B12 over B11
    imRegress = ee.Image.cat(B11, B12)

    linearRegression = imRegress.reduceRegion(
    reducer = ee.Reducer.linearRegression(numX = 1, numY = 1),
    geometry = surrounding10km,
    scale = 20)

    scale_factor = linearRegression.get('coefficients').getInfo()[0]


    ## Calculate R_MBSP
    R_MBSP = img.expression(
    '(c*R12 - R11)/R11 ', # Where this number comes from least squares difference scale above
    {
        
        'c': scale_factor,
        'R11': img.select('B11').divide(10000),
        'R12': img.select('B12').divide(10000),
    },
    )

    ## For now, ignoring the model thing you subtract, so F_MBSP = R_MBSP
    F_MBSP = R_MBSP

    ## Get min and max vals for colorbar in plot
    minval = F_MBSP.reduceRegion(ee.Reducer.min(), geometry = surrounding10km).get('constant').getInfo()
    maxval = F_MBSP.reduceRegion(ee.Reducer.max(), geometry = surrounding10km).get('constant').getInfo()

    return F_MBSP, minval, maxval

# F_MBSP, minval, maxval = retrieve_F_MBSP(source_lon, source_lat, plume_date_start, plume_date_end)
