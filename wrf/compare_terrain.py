"""
Plot Underlying Terrain For Two WRF Simulations and Take the Difference

This script uses raw wrfout files rather than UPP output.

shawn.s.murdzek@noaa.gov
Date Created: 9 December 2022
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

# wrfout files and titles for each simulation
sims = ['/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_1km_CAN_rockies/WPS/geo_em.d01.nc.original',  
        '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_1km_CAN_rockies/WPS/geo_em.d01.nc']  
ttls = ['Original',
        'New']

# Cartopy map projection
proj = ccrs.PlateCarree()

# Grid spacing (m)
dx = 1000.

save_fname = '/mnt/lfs4/BMC/wrfruc/murdzek/figs/terrain_compare_original_smooth.png'

#---------------------------------------------------------------------------------------------------
# Create Plots
#---------------------------------------------------------------------------------------------------

lat = 'XLAT_M'
lon = 'XLONG_M'
hgt = 'HGT_M'

fig = plt.figure(figsize=(8, 10))

ds = {}
axes = []
for i, (s, t) in enumerate(zip(sims, ttls)):
    ds[t] = xr.open_dataset(s)

    # Plot terrain height
    axes.append(fig.add_subplot(3, 1, i+1, projection=proj))
    cax = axes[i].contourf(ds[t][lon][0, :, :], ds[t][lat][0, :, :], ds[t][hgt][0, :, :], 
                           transform=proj, levels=np.arange(0, 2700, 100), cmap='inferno', 
                           extend='max')

    # Determine the maximum gradient
    g = np.gradient(ds[t][hgt][0, :, :], dx)
    gmax = np.amax(np.sqrt((g[0]*g[0]) + (g[1]*g[1])))

    axes[i].set_title('%s\n(max gradient = %.3f)' % (t, gmax), size=12)

# Add colorbar
cbar = plt.colorbar(cax, ax=axes, orientation='vertical')
cbar.set_label('terrain height (m)', size=12)

# Plot difference
ax = fig.add_subplot(3, 1, 3, projection=proj)
cax = ax.contourf(ds[t][lon][0, :, :], ds[t][lat][0, :, :], 
                  ds[ttls[0]][hgt][0, :, :] - ds[ttls[1]][hgt][0, :, :], 
                  transform=proj, levels=np.arange(-1000, 1050, 50), cmap='bwr', extend='both')
ax.set_title('%s $-$ %s' % (ttls[0], ttls[1]), size=12)
cbar2 = plt.colorbar(cax, ax=ax, orientation='vertical')
cbar2.set_label('height difference (m)', size=12)

plt.savefig(save_fname)
plt.close()


"""
End compare_terrain.py
"""
