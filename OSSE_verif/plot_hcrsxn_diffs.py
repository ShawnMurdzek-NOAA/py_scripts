"""
Plot Difference Between NR and RRFS Forecasts

Odd CartoPy Bug: Sometimes, CartoPy fails to plot all the data from RRFS when using a Lambert
Conformal projection (this has been observed to happen with 700-mb T). To fix this problem, simply
use a different map projection, like Plate Carree. I'm not sure why this bug exists, but using
changing the map projection is a good enough workaround.

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
import xarray as xr
import cartopy.crs as ccrs

import pyDA_utils.plot_model_data as pmd


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# UPP files (use pressure-level output, otherwise vertical levels will not match)
NR_files = ['/work2/noaa/wrfruc/murdzek/nature_run_spring/UPP/20220430/wrfprs_202204300100.grib2']
rrfs_files = ['/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data/spring_uas_35km/NCO_dirs/ptmp/prod/rrfs.20220430/00/rrfs.t00z.prslev.f001.conus_3km.grib2']

# Field and level to plot (prslev is in Pa. Set to np.nan for 2D fields)
field = 'TMP_P0_L100_GLC0'
prslev = 70000

# Colorbar and contour levels
cmap = 'plasma'
lvls = np.arange(250, 295, 5)
extend = 'both'

# Domain (set to None to use defaults)
lat_lim = None
lon_lim = None

# Projection
proj = ccrs.PlateCarree()

# Output file name tag
save_tag = ''


#---------------------------------------------------------------------------------------------------
# Plot Data
#---------------------------------------------------------------------------------------------------

# Create plots
for fnameNR, fnameRRFS in zip(NR_files, rrfs_files):
    print('Plotting for %s' % fnameRRFS.split('/')[-1])
    fig = plt.figure(figsize=(8, 11))
    ds = []
    zind = []
    out = []
    for j, (f, n) in enumerate(zip([fnameRRFS, fnameNR], ['RRFS', 'NR'])):
        ds.append(xr.open_dataset(f, engine='pynio'))
        if not np.isnan(prslev):
            zind.append(np.where(ds[-1]['lv_ISBL0'] == prslev)[0][0])
        else:
            zind.append(np.nan)
        out.append(pmd.PlotOutput([ds[-1]], 'upp', fig, 3, 1, j+1, proj=proj))
        out[-1].contourf(field, cntf_kw={'cmap':cmap, 'extend':extend, 'levels':lvls}, 
                         ingest_kw={'zind':[zind[-1]]})
        out[-1].config_ax(grid=False)
        if lat_lim != None:
            out[-1].set_lim(lat_lim[0], lat_lim[1], lon_lim[0], lon_lim[1])
        if not np.isnan(prslev):
            out[-1].ax_title(txt='%d-hPa, %s' % (prslev/100., n), size=14)
        else:
            out[-1].ax_title(txt=n, size=14)
    
    diff = pmd.PlotOutput(ds, 'upp', fig, 3, 1, 3, proj=proj)
    diff.plot_diff(field, auto=False, cntf_kw={'cmap':'bwr', 'extend':'both', 
                                               'levels':np.linspace(-5, 5, 20)},
                                      ingest_kw={'zind':zind, 'indices':[None, [[2, -1, 3], [2, -1, 3]]]})
    diff.config_ax(grid=False)
    if lat_lim != None:
        diff.set_lim(lat_lim[0], lat_lim[1], lon_lim[0], lon_lim[1])
    diff.ax_title(txt='difference (RRFS $-$ NR)', size=14)

    plt.subplots_adjust(left=0.02, bottom=0.04, right=0.98, top=0.97, hspace=0.1)
    
    # Create output file name
    fhr = ds[0][field].attrs['forecast_time'][0]
    init_time = dt.datetime.strptime(ds[0][field].attrs['initial_time'], '%m/%d/%Y (%H:%M)')
    if np.isnan(prslev):
        plt.savefig('RRFSvsNR_%s_%s_f%03d_%s.png' % 
                    (save_tag, init_time.strftime('%Y%m%d%H'), fhr, field.split('_')[0]))
    else:
        plt.savefig('RRFSvsNR_%s_%s_f%03d_%dmb_%s.png' % 
                    (save_tag, init_time.strftime('%Y%m%d%H'), fhr, prslev/1e2, field.split('_')[0]))
    plt.close()


"""
End plot_hcrsxn_diffs.py
""" 
