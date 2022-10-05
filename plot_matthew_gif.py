"""
Plot Horizontal Cross Section from a WRF Simulation

shawn.s.murdzek@noaa.gov
Date Created: 4 October 2022
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import plot_realdata as pr
import numpy as np
import os
import datetime as dt


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

path = '/mnt/lfs4/BMC/wrfruc/murdzek/hurr_matthew_test/WRF/test/em_real/'
ntimes = 17
start = dt.datetime(2016, 10, 6)

for n in range(ntimes):
    fname = path + (start + dt.timedelta(hours=(3*n))).strftime('wrfout_d01_%Y-%m-%d_%H:%M:%S')
    print(fname)

    fig = plt.figure(figsize=(8, 8))

    out = pr.PlotOutput(fname, 'wrf', fig, 1, 1, 1)
    out.contourf('slp', cntf_kw={'cmap':'inferno', 'levels':np.arange(970, 1035, 5), 'extend':'both'})
    out.barbs('ua', 'va', thin=8, ingest_kw={'interp_field':'z', 'interp_lvl':500}, barb_kw={'length':6})
    out.config_ax()
    out.ax_title()

    plt.savefig('tmp%02d.png' % n)

os.system('convert -delay 30 -loop 0 tmp*.png matthew.gif')
os.system('rm *.png')
