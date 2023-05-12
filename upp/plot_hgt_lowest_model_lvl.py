"""
Plot Height of Lowest Model Level AGL

shawn.s.murdzek@noaa.gov
Date Created: 12 May 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.feature as cfeature
import cartopy.crs as ccrs


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

upp_files = ['/work2/noaa/wrfruc/murdzek/nature_run_winter/UPP/20220202/wrfnat_202202022000.grib2']
save_fnames = ['upp_lml_202202022000.png']

# Domain limits
lon = [-135, -60]
lat = [20, 53]

# Single points to plot
lon_pts = []
lat_pts = []


#---------------------------------------------------------------------------------------------------
# Plot Data
#---------------------------------------------------------------------------------------------------

for fname, save_fname in zip(upp_files, save_fnames):
    print('Plotting for %s' % fname)

    ds = xr.open_dataset(fname, engine='pynio')

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    cax = ax.pcolormesh(ds['gridlon_0'], ds['gridlat_0'], 
                        ds['HGT_P0_L105_GLC0'][0, :, :] - ds['HGT_P0_L1_GLC0'], 
                        vmin=2, vmax=8, transform=ccrs.PlateCarree())
    cbar = plt.colorbar(cax, ax=ax, orientation='horizontal')
    cbar.set_label('height of LML (m AGL)', size=14)

    for x, y in zip(lon_pts, lat_pts):
        out.plot(x, y, plt_kw={'markersize':10, 'marker':'*', 'color':'k'})

    ax.set_extent([lon[0], lon[1], lat[0], lat[1]])
    ax.coastlines('10m', edgecolor='k', linewidth=0.75)
    borders = cfeature.NaturalEarthFeature(category='cultural',
                                           scale='10m',
                                           facecolor='none',
                                           name='admin_1_states_provinces')
    ax.add_feature(borders, linewidth=0.75, edgecolor='k')

    plt.suptitle(fname.split('/')[-1], size=16)
    plt.subplots_adjust(left=0.02, bottom=0.04, right=0.98, top=0.97)
    plt.savefig(save_fname)
    plt.close()


"""
End plot_upp.py
""" 
