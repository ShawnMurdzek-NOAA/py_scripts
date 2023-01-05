"""
Plot a 2-D Field For Two WRF Simulations and Take the Difference

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
sims = ['/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_10km_hrrr/WRF/run/wrfout_d01_2021-07-25_00_00_00',  
        '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_10km_hrrr2/WRF/run/wrfout_d01_2021-07-25_00_00_00']  
ttls = ['Deleted Fields',
        'Original']

# Cartopy map projection
proj = ccrs.PlateCarree()

# Grid spacing (m)
field = 'T2'

save_fname = '/mnt/lfs4/BMC/wrfruc/murdzek/figs/compare_T2m.png'

#---------------------------------------------------------------------------------------------------
# Create Plots
#---------------------------------------------------------------------------------------------------

lat = 'XLAT'
lon = 'XLONG'

fig = plt.figure(figsize=(8, 10))

ds = {}
axes = []
for i, (s, t) in enumerate(zip(sims, ttls)):
    ds[t] = xr.open_dataset(s)

    # Plot Field
    axes.append(fig.add_subplot(3, 1, i+1, projection=proj))
    cax = axes[i].contourf(ds[t][lon][0, :, :], ds[t][lat][0, :, :], ds[t][field][0, :, :], 
                           transform=proj, cmap='inferno', extend='both')

    axes[i].set_title(t, size=12)

# Add colorbar
cbar = plt.colorbar(cax, ax=axes, orientation='vertical')
cbar.set_label('%s (%s)' % (ds[t][field].attrs['description'], ds[t][field].attrs['units']), size=12)

# Compute RMSD
rmsd = np.sqrt(np.mean(ds[ttls[0]][field][0, :, :] - ds[ttls[1]][field][0, :, :]))

# Plot difference
ax = fig.add_subplot(3, 1, 3, projection=proj)
cax = ax.contourf(ds[t][lon][0, :, :], ds[t][lat][0, :, :], 
                  ds[ttls[0]][field][0, :, :] - ds[ttls[1]][field][0, :, :], 
                  transform=proj, levels=np.arange(-0.525, 0.526, 0.05), cmap='bwr', extend='both')
ax.set_title('%s $-$ %s (RMSD = %.2e)' % (ttls[0], ttls[1], rmsd), size=12)
cbar2 = plt.colorbar(cax, ax=ax, orientation='vertical')
cbar2.set_label('difference', size=12)

plt.savefig(save_fname)
plt.close()


"""
End compare_field.py
"""
