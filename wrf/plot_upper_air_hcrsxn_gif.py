"""
Plot Horizontal Cross Section from a WRF Simulation

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

path = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_1km_Detroit/WRF/run/'
ntimes = 4
start = dt.datetime(2021, 7, 24, 23, 15)

# Timestep between output files (s)
step = 4

save_fname = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_1km_Detroit/wrf_graphics/ref_1km.gif'


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
    out.contourf('dbz', ingest_kw={'interp_field':'height_agl', 'interp_lvl':1000},
                 cntf_kw={'cmap':'gist_ncar', 'levels':np.arange(5, 75, 5), 'extend':'max'})
    #out.barbs('ua', 'va', thin=10, ingest_kw={'interp_field':'pressure', 'interp_lvl':500}, barb_kw={'length':6})
    out.config_ax()
    out.ax_title()

    plt.savefig('tmp%02d.png' % n)
    plt.close()

os.system('convert -delay 30 -loop 0 tmp*.png %s' % save_fname)
os.system('rm *.png')


"""
End plot_upper_air_hcrsxn_gif.py
"""
