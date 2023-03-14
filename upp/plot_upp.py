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
for mn in range(15, 46, 15):
    upp_files.append('/scratch1/BMC/wrfruc/murdzek/nature_run_spring/wrfnat_2022050419%02d.grib2' % mn)
    save_fnames.append('./upp_REFC_2022050419%02d.png' % mn)

# UPP field to plot
field = 'REFC_P0_L200_GLC0'

# Domain limits
lon = [-74.4, -73.6]
lat = [40.4, 41]

# Single points to plot
lon_pts = [-73.8740, -73.7781, -74.1745]
lat_pts = [40.7769, 40.6413, 40.6895]

name = '1-km WRF'


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
    out.contourf(field, cntf_kw={'cmap':'gist_ncar', 'extend':'max', 'levels':np.arange(5, 76, 5)})
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
