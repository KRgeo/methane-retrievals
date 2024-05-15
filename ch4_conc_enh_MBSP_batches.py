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
import F_MBSP
from F_MBSP import retrieve_F_MBSP

#=-------------------------------------------------------------------------------------
# Define date and location of Hassi Messaoud plume
source_lat = 31.6585
source_lon = 5.9053
plume_date_start = '2019-11-20'  # Define plume_date_start variable

# Get MBSP retrieval for day with plume by calling function in F_MBSP.py
plume_date = plume_date_start  # Use plume_date_start variable
plume_next_day = '2019-11-21'
MBSP, minval, maxval = retrieve_F_MBSP(source_lon, source_lat, plume_date, plume_next_day)

# Get MBSP retrieval for day without plume
ref_img_date = '2019-10-06'
ref_next_day = '2019-10-07'
MBSP_np, minval_np, maxval_np = retrieve_F_MBSP(source_lon, source_lat, ref_img_date, ref_next_day)

## Calculate and plot MBMP retrieval using these two
MBMP = MBSP.subtract(MBSP_np)

## Plot MBMP image
Map = geemap.Map() 
# Define visualisation parameters for plot
vis = {'min': -0.2, 'max' : 0.05, "palette": ['FF0000',	'FFFFFF', '0000FF']}
Map.setCenter(source_lon, source_lat, 13.75) # Define map centre
Map.addLayer(MBMP, vis, 'Sentinel-2') # Add bands of Sentinel-2 data
# Map.addLayer(geometry_area, {}, 'geometry_area') # Overlay rectangle used in Varon (2021) figure 2e
Map.add_colorbar(vis, label='Methane enhancement ')
display(Map)
#--------------------------------------------------------------------------------------------------

# Now we will use the radiative transfer model Radtran to convert from these fractional values of methane
# enhancement to absolute values. Code provided by D. Varon, licensed by GHGSat

# BEcause the model takes a long time to run, we will define many small areas which can be run in parallel
# to get this result. First, we are testing this with a small 3x3 pixel area

## Run retreival using radtran and Varon setup

# Import libraries
import sys

# sys.path.append('Users\s1709837\plumes\my_code\s2_ch4_code_share' )
# from s2_ch4_code_share import hapi.py
import setup
import radtran as rt
import hapi


## Make a rectangle over area around source to calculate absolute concentration enhancement values
size_add = 0.0002
lonmin = source_lon - size_add
lonmax = source_lon + size_add
latmin = source_lat - size_add
latmax = source_lat + size_add
tiny_area = ee.Geometry.Rectangle([lonmin, latmin, lonmax, latmax], None, False)

# Turn list into 2d array and check its shape
F_tiny_area = MBSP.sampleRectangle(tiny_area)
F_array = np.array((F_tiny_area.get('constant').getInfo()))
print('Array shape:',F_array.shape)

# Configuration for radtran
num_layers = 100
targheight = 0
obsheight = 100
solarangle = 40
obsangle = 0
instrument = 'S2A'
method = 'MBSP'

# Invented 2D array of fractional reflectance data
frac_refl_data = F_array


from tqdm import tqdm
import time

with tqdm(range(F_array.size)) as pbar:
    for i in tqdm(range(F_array.size), desc = 'tqdm() Progress Bar'):
        plume_retrieval = rt.retrieve(frac_refl_data, instrument, method, targheight, obsheight, solarangle, obsangle, num_layers=num_layers)
        # time.sleep(0.5)
        pbar.update()

print(plume_retrieval)


#------------------------------------------------------------------------------

# Same for "background" image

# Configuration for radtran
num_layers = 100
targheight = 0
obsheight = 100
solarangle = 40
obsangle = 0
instrument = 'S2A'
method = 'MBSP'

# Invented 2D array of fractional reflectance data
# Turn list into 2d array and check its shape
F_tiny_area = MBSP_np.sampleRectangle(tiny_area)
F_array = np.array((F_tiny_area.get('constant').getInfo()))
print('Array shape:',F_array.shape)
frac_refl_data = F_array


from tqdm import tqdm
import time

with tqdm(range(F_array.size)) as pbar:
    for i in tqdm(range(F_array.size), desc = 'tqdm() Progress Bar'):
        np_retrieval = rt.retrieve(frac_refl_data, instrument, method, targheight, obsheight, solarangle, obsangle, num_layers=num_layers)
        time.sleep(0.5)
        pbar.update()

print(np_retrieval)


# Calculate and plot using these two
MBMP_retrieval = plume_retrieval - np_retrieval
print(MBMP_retrieval)
