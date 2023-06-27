"""
Plot Underlying Terrain For a WRF Simulation

This script uses raw wrfout files rather than UPP output.

shawn.s.murdzek@noaa.gov
Date Created: 4 October 2022
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import os
import datetime as dt
import netCDF4 as nc

import pyDA_utils.plot_model_data as pmd


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

path = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_1km_CAN_rockies_smooth_terrain/WRF/run/'
start = dt.datetime(2021, 7, 24, 23)

# Additional points to plot
pts={'lat':[50.33566],
     'lon':[-122.63756],
     'c':['b'],
     'm':['o'],
     'size':[10]}

save_fname = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_1km_CAN_rockies_smooth_terrain/wrf_graphics/terrain.png'

#---------------------------------------------------------------------------------------------------
# Create Plots
#---------------------------------------------------------------------------------------------------

fname = path + start.strftime('wrfout_d01_%Y-%m-%d_%H_%M_%S')

fig = plt.figure(figsize=(8, 8))

fptr = nc.Dataset(fname)
out = pmd.PlotOutput([fptr], 'wrf', fig, 1, 1, 1)
out.contourf('HGT', cntf_kw={'cmap':'inferno', 'levels':np.arange(0, 2700, 100), 'extend':'max'})
for i in range(len(pts['lat'])):
    out.plot(pts['lon'][i], pts['lat'][i], plt_kw={'c':pts['c'][i], 'marker':pts['m'][i],
                                                   'markersize':pts['size'][i]})
out.config_ax()
out.ax_title()

plt.savefig(save_fname)
plt.close()
