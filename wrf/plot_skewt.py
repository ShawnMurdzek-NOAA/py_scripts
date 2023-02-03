"""
Plot Skew-Ts from a WRF Simulation

This script uses raw wrfout files rather than UPP output.

shawn.s.murdzek@noaa.gov
Date Created: 8 December 2022
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import plot_model_data as pmd
import numpy as np
import os
import datetime as dt
import netCDF4 as nc


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

#path = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_3km/WRF/run/'
path = '/scratch1/BMC/wrfruc/murdzek/nature_run_3km/WRF/run/'
ntimes = 2
start = dt.datetime(2021, 7, 24, 23)

# Timestep between output files (s)
step = 900

# Location to plot Skew-Ts
lat = 41.7658
lon = -72.6734

save_fname = '/scratch1/BMC/wrfruc/murdzek/py_scripts/wrf/skewts.gif'

#---------------------------------------------------------------------------------------------------
# Create Plots
#---------------------------------------------------------------------------------------------------

# Notes:
#     For CONUS, thin = 80 is good for the wind barbs

for n in range(ntimes):
    fname = path + (start + dt.timedelta(seconds=(step*n))).strftime('wrfout_d01_%Y-%m-%d_%H_%M_%S')
    print(fname)

    fig = plt.figure(figsize=(8, 8))

    fptr = nc.Dataset(fname)
    out = pmd.PlotOutput([fptr], 'wrf', fig, 1, 1, 1)
    out.skewt(lon, lat, barbs=True)

    plt.savefig('tmp%02d.png' % n)
    plt.close()

os.system('convert -delay 30 -loop 0 tmp*.png %s' % save_fname)
os.system('rm *.png')
