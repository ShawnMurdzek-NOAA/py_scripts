"""
Plot IMS Snow Cover Fraction

shawn.s.murdzek@noaa.gov
Date Created: 2 May 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import xarray as xr
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
import plot_model_data as pmd


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

#template = '/work2/noaa/wrfruc/murdzek/RRFS_input_data/winter/snow/ims96/grib2/22%03d2200000000.grib2'
template = '/work2/noaa/wrfruc/murdzek/nature_run_spring/synthetic_ims/22%03d2200000000.grib2'
fnames = [template % i for i in range(119, 127)]
#fnames = [template % i for i in range(32, 39)]

minlon = -130.
maxlon = -60.
minlat = 20.
maxlat = 55.

save_dir = './'


#---------------------------------------------------------------------------------------------------
# Plot Data
#---------------------------------------------------------------------------------------------------

for f in fnames:
    print('Plotting for %s' % f)

    ds = xr.open_dataset(f, engine='pynio')
    fig = plt.figure(figsize=(8, 4))
    plt.subplots_adjust(left=0.02, bottom=0.04, right=0.98, top=0.97)
    out = pmd.PlotOutput([ds], 'upp', fig, 1, 1, 1)
    out.pcolormesh('SNOWC_P0_L1_GST0')
    out.config_ax(grid=False)
    out.ax_title(size=14)
    out.set_lim(minlat, maxlat, minlon, maxlon)

    plt.savefig('%s/ims_snow_NR_%s.png' % (save_dir, f.split('/')[-1].split('.')[0]))
    plt.close()


"""
End plot_ims_snow.py
""" 
