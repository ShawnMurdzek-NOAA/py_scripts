"""
Plot WRF Postage Stamp Plots

Parameters
----------
sys.argv[1] : Path containing ensemble output
sys.argv[2] : YYYYMMDDHHMMSS timestamp of ensemble output
sys.argv[3] : Output image file name

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import sys
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import pyart.graph as art_cm

import pyDA_utils.plot_postage_stamp as pps


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Path containing ensemble output
path = sys.argv[1]

# Names of ensemble members
names = ['ctrl'] + [f"skeb{n}" for n in range(1,18)]

# WRF output file name
wrf_time = dt.datetime.strptime(sys.argv[2], '%Y%m%d%H%M%S')
wrf_fname = dt.datetime.strftime(wrf_time, 'wrfout_d01_%Y-%m-%d_%H:%M:%S')

# File names containing ensemble output
fnames = [[f"{path}/ctrl/WRF/run/{wrf_fname}"]] + [[f"{path}/{n}/run/{wrf_fname}"] for n in names[1:]]

# Output filename
out_fname = sys.argv[3]

# Text size
size = 16


#---------------------------------------------------------------------------------------------------
# Main Program
#---------------------------------------------------------------------------------------------------

# Create postage stamp plot object
ps_obj = pps.PlotPostageStamp(fnames, names, 'wrf', 4, 5, fig_kw={'figsize':(12, 10)}, proj=None)

# Plot composite reflectivity
lvls = np.arange(0, 75, 5)
ps_obj.contourf('mdbz', label_kw={'size':size}, 
                kw={'cntf_kw':{'levels':lvls, 'cmap':'HomeyerRainbow'}})

# Add figure annotations
ps_obj.ax_titles(size=size)
ps_obj.sup_title(size=(size+4))

# Eliminate some whitespace
plt.subplots_adjust(left=0.04, bottom=0.03, right=0.90, top=0.93, wspace=0.1, hspace=0.15)

# Save figure
plt.savefig(out_fname, dpi=500)


"""
End plot_hcrsxn_postage_stamp.py
"""
