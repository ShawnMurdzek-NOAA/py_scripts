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


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

fname = '/mnt/lfs4/BMC/wrfruc/murdzek/hurr_matthew_test/WRF/test/em_real/wrfout_d01_2016-10-06_03:00:00'

fig = plt.figure(figsize=(8, 8))

out = pr.PlotOutput(fname, 'wrf', fig, 1, 1, 1)

#out.contourf('slp')
out.contourf('z', units='dm', interp_field='pressure', interp_lvl=500)

out.config_ax()
out.ax_title()

plt.savefig('test.png')
