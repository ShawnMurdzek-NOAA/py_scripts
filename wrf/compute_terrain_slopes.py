"""
Compute Maximum Terrain Slope for Several Simulations

This script uses raw wrfout files rather than UPP output.

shawn.s.murdzek@noaa.gov
Date Created: 12 December 2022
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
import numpy as np
import os
import datetime as dt


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Input files, sim names, and terrain height fields
sims = {}
for i in range(24):
    sims['prc%d' % i] = {'file': '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_1km_CONUS/WPS/geo_em.d01.nc_%04d' % i,
                         'field': 'HGT_M'}

# Horizontal grid spacing (m)
dx = 1000.


#---------------------------------------------------------------------------------------------------
# Compute Terrain Slope
#---------------------------------------------------------------------------------------------------

for n in sims.keys():
    ds = xr.open_dataset(sims[n]['file'])

    # Determine the maximum gradient
    g = np.gradient(ds[sims[n]['field']][0, :, :], dx)
    mag = np.sqrt((g[0]*g[0]) + (g[1]*g[1]))
    gmax = np.amax(mag)
    gmaxloc = np.unravel_index(np.argmax(mag), ds[sims[n]['field']][0, :, :].shape)

    print('%s, max slope = %.3f at (%d, %d)' % (n, gmax, gmaxloc[0], gmaxloc[1]))


"""
End compute_terrain_slopes.py
"""
