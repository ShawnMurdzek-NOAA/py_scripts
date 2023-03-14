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
import cartopy.feature as cfeature
import matplotlib.pyplot as plt


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Input netCDF file
infile = '/scratch1/BMC/wrfruc/murdzek/nature_run_winter/output/202202021300/JOIN/wrfout_d01_2022-02-02_14:00:00'

# Plot parameters
field = 'LANDMASK'
#proj = ccrs.LambertConformal()  # For some reason, the Lambert Conformal grid doesn't work
proj = ccrs.PlateCarree()

# Domain limits
lon = [-74.4, -73.6]
lat = [40.4, 41]

# Output file
save_fname = './landmask_1km.png'


#---------------------------------------------------------------------------------------------------
# Plot Field
#---------------------------------------------------------------------------------------------------

ds = xr.open_dataset(infile)

# Create figure
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, projection=proj)

#cax = ax.contourf(ds['XLONG'][0, :, :], ds['XLAT'][0, :, :], ds[field][0, :, :], 
#                  np.arange(0, 1.1), cmap='winter')
cax = ax.pcolormesh(ds['XLONG'][0, :, :], ds['XLAT'][0, :, :], ds[field][0, :, :], 
                    vmin=0, vmax=1, cmap='winter', transform=proj)

#cbar = plt.colorbar(cax, ax=ax, orientation='horizontal')
#cbar.set_label('%s (%s)' % (ds[field].attrs['description'], ds[field].attrs['units']), size=12)

# Add coastlines and state borders
ax.coastlines('10m', linewidth=0.75)
borders = cfeature.NaturalEarthFeature(category='cultural',
                                       scale='10m',
                                       facecolor='none',
                                       name='admin_1_states_provinces')
ax.add_feature(borders, linewidth=0.75, edgecolor='k')

ax.set_extent([lat[0], lat[1], lon[0], lon[1]])

plt.suptitle(str(ds['XTIME'][0].values).split('.')[0], size=16)
plt.savefig(save_fname)
plt.close()


"""
End plot_hcrsxn_xarray.py 
"""
