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

path = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_3km/WRF/run/'
ntimes = 5
start = dt.datetime(2021, 7, 24, 23)

# Timestep between output files (min)
step = 15

save_fname = '/mnt/lfs4/BMC/wrfruc/murdzek/figs/nature_run_3km/500mb_z_winds.gif'

#---------------------------------------------------------------------------------------------------
# Create Plots
#---------------------------------------------------------------------------------------------------

for n in range(ntimes):
    fname = path + (start + dt.timedelta(minutes=(step*n))).strftime('wrfout_d01_%Y-%m-%d_%H_%M_%S')
    print(fname)

    fig = plt.figure(figsize=(8, 8))

    out = pmd.PlotOutput(fname, 'wrf', fig, 1, 1, 1)
    out.contourf('z', ingest_kw={'interp_field':'pressure', 'interp_lvl':500},
                 cntf_kw={'cmap':'inferno', 'levels':np.arange(5500, 6100, 25), 'extend':'both'})
    out.barbs('ua', 'va', thin=80, ingest_kw={'interp_field':'pressure', 'interp_lvl':500}, barb_kw={'length':6})
    out.config_ax()
    out.ax_title()

    plt.savefig('tmp%02d.png' % n)

os.system('convert -delay 30 -loop 0 tmp*.png %s' % save_fname)
os.system('rm *.png')
