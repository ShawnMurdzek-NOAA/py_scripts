"""
Smooth Terrain Field from Geogrid

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import numpy as np
from scipy.ndimage import convolve


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Input geo_em file
in_fname = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_1km_CAN_rockies_py_smooth/WPS/geo_em.d01.nc.original'

# Output geo_em file
out_fname = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_1km_CAN_rockies_py_smooth/WPS/geo_em.d01.nc.1'

# Smoothing option
opt = '1-2-1'
npass = 1

# Gridspacing (m)
dx = 1000.


#---------------------------------------------------------------------------------------------------
# Apply Smoothing
#---------------------------------------------------------------------------------------------------

ds = xr.open_dataset(in_fname)
field = 'HGT_M'

for i in range(npass):
    print('pass = %d' % (i+1))
    if opt == '1-2-1':
        ds[field][0, :, :] = convolve(ds[field][0, :, :], 
                                      np.array([[0, 1, 0], [1, 2, 1], [0, 1, 0]])) / 6.

# Compute maximum gradient
g = np.gradient(ds[field][0, :, :], dx)
mag = np.sqrt((g[0]*g[0]) + (g[1]*g[1]))
gmax = np.amax(mag)
loc_gmax = np.unravel_index(np.argmax(mag), ds[field][0, :, :].shape)
print('max gradient = %.3f at (%d, %d)' % (gmax, loc_gmax[0], loc_gmax[1]))

ds.to_netcdf(out_fname)


"""
End smooth_terrain.py
""" 
