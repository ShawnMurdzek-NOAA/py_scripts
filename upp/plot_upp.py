"""
Plot UPP Output Fields

shawn.s.murdzek@noaa.gov
Date Created: 27 January 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
import plot_model_data as pmd
import xarray as xr
import cartopy.feature as cfeature


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

upp_files = []
save_fnames = []
for day in ['0201', '0202', '0203', '0204', '0205', '0206', '0207']:
    upp_files = upp_files + ['/work2/noaa/wrfruc/murdzek/nature_run_winter/UPP/2022%s/wrfnat_2022%s2200.grib2' 
                             % (day, day)]   
    save_fnames = save_fnames + ['upp_snow_2022%s2200.png' % day]

#upp_files = ['/work2/noaa/wrfruc/murdzek/nature_run_spring/UPP/20220503/wrfnat_202205032200.grib2']
#save_fnames = ['upp_ceil_202205032200.png']

# UPP field to plot
field = 'SNOD_P0_L1_GLC0'

# Domain limits
lon = [-130, -60]
lat = [20, 55]

# Single points to plot
lon_pts = []
lat_pts = []

name = 'Cycled Snow Depth'


#---------------------------------------------------------------------------------------------------
# Plot Data
#---------------------------------------------------------------------------------------------------

# Sample colors from plasma colorbar and define plotting levels
#plevels = np.array([0.01, 0.1, 0.5, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]) * 1e-3
#nlvl = len(plevels)
#colors = [mcm.plasma(i / (nlvl-1)) for i in range(nlvl)]

for fname, save_fname in zip(upp_files, save_fnames):
    print('Plotting for %s' % fname)

    fig = plt.figure(figsize=(8, 8))
    plt.subplots_adjust(left=0.02, bottom=0.04, right=0.98, top=0.97)
    ds = xr.open_dataset(fname, engine='pynio')
    out = pmd.PlotOutput([ds], 'upp', fig, 1, 1, 1)
    out.pcolormesh(field, pcm_kw={'cmap':'plasma', 'vmin':0, 'vmax':2})
    #out.barbs('UGRD_P0_L103_GLC0', 'VGRD_P0_L103_GLC0', thin=5, ingest_kw={'zind':0}, barb_kw={})

    for x, y in zip(lon_pts, lat_pts):
            out.plot(x, y, plt_kw={'markersize':10, 'marker':'*', 'color':'k'})

    out.set_lim(lat[0], lat[1], lon[0], lon[1])
    #out.config_ax(grid=False)
    out.ax.coastlines('10m', edgecolor='k', linewidth=0.75)
    borders = cfeature.NaturalEarthFeature(category='cultural',
                                           scale='10m',
                                           facecolor='none',
                                           name='admin_1_states_provinces')
    out.ax.add_feature(borders, linewidth=0.75, edgecolor='k')
    out.ax_title(txt=name, size=14)

    plt.savefig(save_fname)
    plt.close()


"""
End plot_upp.py
""" 
