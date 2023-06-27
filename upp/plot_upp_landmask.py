"""
Plot The Model Landmask from Two Different Simulations

shawn.s.murdzek@noaa.gov
Date Created: 3 February 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
import xarray as xr
import cartopy.feature as cfeature

import pyDA_utils.plot_model_data as pmd


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# UPP output files
fname1 = '/scratch1/BMC/wrfruc/murdzek/nature_run_spring/wrfnat_202205041900.grib2'
fname2 = '/scratch1/BMC/wrfruc/murdzek/HRRR_data/2211911000001.grib2'

# Field and level to plot (set level to np.nan for a 2D field)
field = 'LAND_P0_L1_GLC0'
zind = np.nan

# Colorbar and contour levels
cmap = 'winter'
vmin = 0
vmax = 1
extend = 'neither'

# Titles
name1 = '1-km WRF'
name2 = '3-km HRRR'

# Domain limits
lon = [-74.4, -73.6]
lat = [40.4, 41]

# Single points to plot
lon_pts = [-73.8740, -73.7781, -74.1745]
lat_pts = [40.7769, 40.6413, 40.6895]

# Output file name
save_fname = './landmask_1km_3km.png'


#---------------------------------------------------------------------------------------------------
# Plot Data
#---------------------------------------------------------------------------------------------------

print('Plotting for %s' % fname1)
fig = plt.figure(figsize=(10, 6))
ds = []
for j, (f, n) in enumerate(zip([fname1, fname2], [name1, name2])):
    ds.append(xr.open_dataset(f, engine='pynio'))
    out = pmd.PlotOutput([ds[-1]], 'upp', fig, 1, 2, j+1)
    out.pcolormesh(field, pcm_kw={'cmap':cmap, 'vmin':vmin, 'vmax':vmax}, 
                   ingest_kw={'zind':zind}, cbar_kw={'orientation':'horizontal', 'aspect':20})
    for x, y in zip(lon_pts, lat_pts):
        out.plot(x, y, plt_kw={'markersize':10, 'marker':'*', 'color':'k'})
    out.ax.coastlines('10m', edgecolor='w', linewidth=0.75)
    borders = cfeature.NaturalEarthFeature(category='cultural',
                                           scale='10m',
                                           facecolor='none',
                                           name='admin_1_states_provinces')
    out.ax.add_feature(borders, linewidth=0.75, edgecolor='w')
    out.ax_title(txt=n, size=14)
    out.set_lim(lat[0], lat[1], lon[0], lon[1])

plt.subplots_adjust(left=0.02, bottom=0.04, right=0.98, top=0.93, wspace=0.05)
plt.savefig(save_fname)
plt.close()


"""
End plot_upp_landmask.py
""" 
