"""
Plot A Potential OSSE Case Study: Line of Convection in S Texas

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import pyart.graph.cm_colorblind as art_cm
import datetime as dt

import pyDA_utils.plot_model_data as pmd


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Figure size and files to include
figsize = (8, 6)
#files = np.array([['/work2/noaa/wrfruc/murdzek/nature_run_spring/UPP/20220504/wrfprs_202205041600_er.grib2',
#                   '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data_rrfs-workflow_orion/spring/NCO_dirs/ptmp/prod/rrfs.20220504/16/rrfs.t16z.prslev.f000.conus_3km.grib2',
#                   '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data_rrfs-workflow_orion/spring_uas_35km/NCO_dirs/ptmp/prod/rrfs.20220504/16/rrfs.t16z.prslev.f000.conus_3km.grib2'],
#                  ['/work2/noaa/wrfruc/murdzek/nature_run_spring/UPP/20220504/wrfprs_202205041900_er.grib2',
#                   '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data_rrfs-workflow_orion/spring/NCO_dirs/ptmp/prod/rrfs.20220504/16/rrfs.t16z.prslev.f003.conus_3km.grib2',
#                   '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data_rrfs-workflow_orion/spring_uas_35km/NCO_dirs/ptmp/prod/rrfs.20220504/16/rrfs.t16z.prslev.f003.conus_3km.grib2']])
#titles = np.array([['Nature Run 16:00 UTC', 'Ctrl f00 hr', '35-km UAS f00 hr'],
#                   ['Nature Run 19:00 UTC', 'Ctrl f03 hr', '35-km UAS f03 hr']])
files = np.array([['/work2/noaa/wrfruc/murdzek/nature_run_spring/UPP/20220504/wrfprs_202205041500_er.grib2',
                   '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data_rrfs-workflow_orion/spring/NCO_dirs/ptmp/prod/rrfs.20220504/15/rrfs.t15z.prslev.f000.conus_3km.grib2',
                   '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data_rrfs-workflow_orion/spring_uas_35km/NCO_dirs/ptmp/prod/rrfs.20220504/15/rrfs.t15z.prslev.f000.conus_3km.grib2'],
                  ['/work2/noaa/wrfruc/murdzek/nature_run_spring/UPP/20220504/wrfprs_202205041900_er.grib2',
                   '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data_rrfs-workflow_orion/spring/NCO_dirs/ptmp/prod/rrfs.20220504/15/rrfs.t15z.prslev.f004.conus_3km.grib2',
                   '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data_rrfs-workflow_orion/spring_uas_35km/NCO_dirs/ptmp/prod/rrfs.20220504/15/rrfs.t15z.prslev.f004.conus_3km.grib2']])
titles = np.array([['Nature Run 15:00 UTC', 'Ctrl f00 hr', '35-km UAS f00 hr'],
                   ['Nature Run 19:00 UTC', 'Ctrl f04 hr', '35-km UAS f04 hr']])

# Field plotting parameters
field = 'REFC_P0_L200_GLC0'
field_name = 'composite reflectivity (dB$Z$)'
cmap = art_cm.HomeyerRainbow
clevels = np.arange(5, 80, 5)

# Subdomain (lat, lon) coordinates
lat = [28, 36]
lon = [-105, -96]

# Output file
out_fname = 'S_TX_2022050419_f04.png'

# Option for debugging (only use one GRIB2 file to make the code run faster)
debug = False


#---------------------------------------------------------------------------------------------------
# Make Plot
#---------------------------------------------------------------------------------------------------

start = dt.datetime.now()
print("Start time =", start)

fig = plt.figure(figsize=figsize)
nrows = files.shape[0]
ncols = files.shape[1]
out_obj = []

if debug:
    fname = files[1, 2]
    ds = xr.open_dataset(fname, engine='pynio')

for i in range(nrows):
    for j in range(ncols):
        fname = files[i, j]
        print()
        print(f"opening {fname}")
        if not debug:
            ds = xr.open_dataset(fname, engine='pynio')

        print("creating pmd object")
        out = pmd.PlotOutput([ds], 'upp', fig, nrows, ncols, (i*3) + j + 1, 
                             proj=ccrs.PlateCarree())

        print("making plot")
        out.contourf(field, cbar=False, cntf_kw={'cmap':cmap, 'levels':clevels})

        print("adding annotations")
        out.set_lim(lat[0], lat[1], lon[0], lon[1])
        out.ax.coastlines('10m', edgecolor='k', linewidth=0.75)
        borders = cfeature.NaturalEarthFeature(category='cultural',
                                               scale='10m',
                                               facecolor='none',
                                               name='admin_1_states_provinces')
        out.ax.add_feature(borders, linewidth=0.75, edgecolor='k')
        out.ax.set_title(titles[i][j], size=14)
        out_obj.append(out)

# Add colorbar
cb_ax = fig.add_axes([0.05, 0.1, 0.9, 0.03])
cbar = plt.colorbar(out_obj[0].cax, cax=cb_ax, orientation='horizontal', aspect=35, pad=0.05)
cbar.set_label(field_name, size=14)
cbar.ax.tick_params(labelsize=10)

plt.subplots_adjust(left=0.01, bottom=0.15, right=0.99, top=0.95, hspace=0.2, wspace=0.02)
plt.savefig(out_fname, dpi=300)
plt.close()

print()
print(f"Elapsed time = {(dt.datetime.now() - start).total_seconds()} s")


"""
End plot_OSSE_case_study.py
"""
