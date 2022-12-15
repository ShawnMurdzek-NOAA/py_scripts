"""
Plot a Horizontal Cross Section Using XArray

Quick script to plot a 2D field directly from a wrfout file. Designed to not be memory intensive so
that fields can be plotted from really big netCDF files

shawn.s.murdzek@noaa.gov
Date Created: 15 December 2022
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Input netCDF file
infile = '/mnt/lfs4/BMC/wrfruc/murdzek/JOINER/wrfout_d01_2021-07-24_23:00:00'

# Plot parameters
field = 'HGT'
proj = ccrs.PlateCarree()

# Output file
save_fname = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_1km_CONUS/wrf_graphics/hgt.png'


#---------------------------------------------------------------------------------------------------
# Plot Field
#---------------------------------------------------------------------------------------------------

ds = xr.open_dataset(infile)

# Create figure
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, proj=proj)

cax = ax.contourf(ds['XLONG'][0, :, :], ds['XLAT'][0, :, :], ds[field][0, :, :], 
                  np.arange(0, 3000, 100), cmap='inferno', extend='max')

cbar = plt.colorbar(cax, ax=ax, orientation='horizontal')
cbar.set_label('%s (%s)' % (ds[field].attrs['description'], ds[field].attrs['units']), size=12)

# Add coastlines and state borders
ax.coastlines('50m')
borders = cfeature.NaturalEarthFeature(category='cultural',
                                       scale='50m',
                                       facecolor='none',
                                       name='admin_1_states_provinces')
ax.add_feature(borders, linewidth=0.5, edgecolor='k')

plt.suptitle(str(ds['XTIME'][0].values).split('.')[0], size=16)
plt.savefig(save_fname)
plt.close()


"""
End plot_hcrsxn_xarray.py 
"""
