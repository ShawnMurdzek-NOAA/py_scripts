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
import plot_model_data as pmd
import numpy as np
import os
import datetime as dt


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

path = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_1km_CAN_rockies_smooth_terrain/WRF/run/'
ntimes = 19
start = dt.datetime(2021, 7, 24, 23)

# Timestep between output files (s)
step = 5

# Plotting bounds (set to None to use defaults)
lat_lims = [50.2, 50.4]
lon_lims = [-122.7, -122.5]

save_fname = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_1km_CAN_rockies_smooth_terrain/wrf_graphics/psfc_wnd10_zoomed.gif'

#---------------------------------------------------------------------------------------------------
# Create Plots
#---------------------------------------------------------------------------------------------------

# Notes:
#     For CONUS, thin = 80 is good for the wind barbs

for n in range(ntimes):
    fname = path + (start + dt.timedelta(seconds=(step*n))).strftime('wrfout_d01_%Y-%m-%d_%H_%M_%S')
    print(fname)

    fig = plt.figure(figsize=(8, 8))

    out = pmd.PlotOutput(fname, 'wrf', fig, 1, 1, 1)
    out.contourf('PSFC', ingest_kw={},
                 cntf_kw={'cmap':'inferno', 'levels':np.arange(68000, 104000, 1500), 'extend':'both'})
    out.barbs('U10', 'V10', thin=1, ingest_kw={}, barb_kw={'length':6})
    out.config_ax()
    out.ax_title()

    if lat_lims != None:
        out.set_lim(lat_lims[0], lat_lims[1], lon_lims[0], lon_lims[1])

    plt.savefig('tmp%02d.png' % n)
    plt.close()

os.system('convert -delay 30 -loop 0 tmp*.png %s' % save_fname)
os.system('rm *.png')
