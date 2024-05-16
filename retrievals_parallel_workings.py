#----------------------------------------------------------------------------------------------------------------------------------# 
# Kalyani Ramanan 15/05/24

# Based on ch4_conc_enh_MBSP. Get methane concentration enhancement values in mol/m2 from S2 B11 and B12 images. Use F_MBSP function to claculate MBSP transmittance maps (as per Varon 2021) to input to radtran model.
#
# Split to parallel jobs using MPI (?)
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
import itertools
import utm
import webbrowser
# from mpi4py import MPI

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

# Create square geometry area centered on source_lon, source_lat with 10km diameter
size_add = 4000  # 10km in meters
minlon = source_lon - size_add / 111000
maxlon = source_lon + size_add / 111000
minlat = source_lat - size_add / (111000 * np.cos(np.deg2rad(source_lat)))
maxlat = source_lat + size_add / (111000 * np.cos(np.deg2rad(source_lat)))
surrounding8km = ee.Geometry.Rectangle([minlon, minlat, maxlon, maxlat], None, False)

#com out for now
## Plot MBSP image
Map = geemap.Map() 
# Define visualisation parameters for plot
vis = {'min': -0.5, 'max' : 0.05, "palette": ['FF0000',	'FFFFFF', '0000FF']}
Map.setCenter(source_lon, source_lat, 13.75) # Define map centre
Map.addLayer(MBSP, vis, 'Sentinel-2') # Add bands of Sentinel-2 data
Map.addLayer(surrounding8km, {}, 'geometry_area') # Overlay rectangle used in Varon (2021) figure 2e
Map.add_colorbar(vis, label='Methane enhancement ')
# display(Map)
# Save the map to an HTML file
output_file = 'output_map.html'
Map.to_html(output_file)

# Open the map in a web browser
webbrowser.open(output_file)

# MBSP_array = MBSP.sampleRectangle(surrounding10km, defaultValue=0).get('constant').getInfo()
# print(MBSP_array)



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
import itertools


##----------------------------------------------------
# Define input data for radtran

# Configuration for radtran
num_layers = 100
targheight = 0
obsheight = 100
solarangle = 40
obsangle = 0
instrument = 'S2A'
method = 'MBSP'

# Invented 2D array of fractional reflectance data
MBSP_array = MBSP.sampleRectangle(surrounding8km).get('constant').getInfo()
frac_refl_data = MBSP_array


# np_retrieval = rt.retrieve(frac_refl_data, instrument, method, targheight, obsheight, solarangle, obsangle, num_layers=num_layers)


class data_container:
        """
                Put these data into a class to contain them as one object as they get passed to a 
                function.
        """
        def __init__(self, frac_refl_data, instrument, method, targheight, obsheight, solarangle, obsangle, num_layers):
                self.frac_refl_data = frac_refl_data
                self.instrument = instrument
                self.method = method
                self.targheight = targheight
                self.obsheight = obsheight
                self.solarangle = solarangle
                self.obsangle = obsangle
                self.num_layers = num_layers






# Create object containing all relevant input data to be passed to radtran
data_input =  data_container(frac_refl_data, instrument, method, targheight, obsheight, solarangle, obsangle, num_layers)
index_row = list(range(0, np.shape(frac_refl_data)[0]))
index_column = list(range(0, np.shape(frac_refl_data)[1]))

# Create iterator that is the equivalent of a nested loop list
index_iterator = list(itertools.product(index_row, index_column))
index_iterator = [(int(i), int(j)) for i, j in index_iterator]  # Convert tuple objects to integers

# Find out which MPI rank with process is
comm = MPI.COMM_WORLD
processes = comm.Get_size() #MPI size
mpi_rank = comm.Get_rank() #MPI rank
print("mpi size: ", processes, "mpi rank: ", mpi_rank)


def compute_trials_ind(num_trials_total, mpi_size, mpi_rank):
    """Boilerplate even(ish) division of workload logic"""
    num_trials_this_rank = num_trials_total // mpi_size
    remainder = num_trials_total % mpi_size
    start_ind = num_trials_total // mpi_size * (mpi_rank)
    end_ind = num_trials_total // mpi_size * (mpi_rank + 1)
    if (mpi_rank) == mpi_size:
        end_ind = -1
    return start_ind, end_ind

start_ind, end_ind = compute_trials_ind(len(index_iterator), processes, mpi_rank)














# from tqdm import tqdm
# import time

# with tqdm(range(F_array.size)) as pbar:
#     for i in tqdm(range(F_array.size), desc = 'tqdm() Progress Bar'):
#         plume_retrieval = rt.retrieve(frac_refl_data, instrument, method, targheight, obsheight, solarangle, obsangle, num_layers=num_layers)
#         # time.sleep(0.5)
#         pbar.update()

# print(plume_retrieval)




# ##-----------------------------------------------------------------------------------------------------------------------

