"""
Plot Difference Between Two Modeling Systems

Useful for comparing the NR to RRFS forecasts

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
# Also specify indices to be passed to _ingest_data for difference plot

msu_home = '/work2/noaa/wrfruc/murdzek/'
model_info = {'RRFS':{'files':[msu_home + 'RRFS_OSSE/syn_data/winter_updated/NCO_dirs/ptmp/prod/rrfs.20220203/12/rrfs.t12z.prslev.f000.conus_3km.grib2',
                               msu_home + 'RRFS_OSSE/syn_data/winter_updated/NCO_dirs/ptmp/prod/rrfs.20220203/11/rrfs.t11z.prslev.f001.conus_3km.grib2',
                               msu_home + 'RRFS_OSSE/syn_data/winter_updated/NCO_dirs/ptmp/prod/rrfs.20220203/09/rrfs.t09z.prslev.f003.conus_3km.grib2'],
                      'indices':None},
              'NR':{'files':3*[msu_home + 'nature_run_winter/UPP/20220203/wrfprs_202202031200.grib2'],
                    'indices':[[2, -1, 3], [2, -1, 3]]}}

model_info = {'RRFS_real':{'files':[msu_home + 'RRFS_OSSE/real_red_data/spring/NCO_dirs/ptmp/prod/rrfs.20220505/12/rrfs.t12z.prslev.f000.conus_3km.grib2',
                                    msu_home + 'RRFS_OSSE/real_red_data/spring/NCO_dirs/ptmp/prod/rrfs.20220505/11/rrfs.t11z.prslev.f001.conus_3km.grib2',
                                    msu_home + 'RRFS_OSSE/real_red_data/spring/NCO_dirs/ptmp/prod/rrfs.20220505/09/rrfs.t09z.prslev.f003.conus_3km.grib2'],
                      'indices':None},
              'HRRR':{'files':3*[msu_home + 'HRRR_data/20220505/hrrr.t12z.wrfprsf00.grib2'],
                      'indices':None}}

model_info = {'RRFS':{'files':[msu_home + 'RRFS_data/20240104/rrfs.t12z.prslev.f000.conus_3km.grib2',
                               msu_home + 'RRFS_data/20240104/rrfs.t11z.prslev.f001.conus_3km.grib2',
                               msu_home + 'RRFS_data/20240104/rrfs.t09z.prslev.f003.conus_3km.grib2'],
                      'indices':None},
              'HRRR':{'files':3*[msu_home + 'HRRR_data/20240104/hrrr.t12z.wrfprsf00.grib2'],
                      'indices':None}}

# Field and level to plot (prslev is in Pa. Set to np.nan for 2D fields)
# In addition to UPP output fields, 'wspd' and 'wdir' are options
field = 'gridlat'
prslev = 25000

# Domain (set to None to use defaults)
lat_lim = None
lon_lim = None

# Projection
proj = ccrs.PlateCarree()

# Output file name tag
save_tag = ''


#---------------------------------------------------------------------------------------------------
# Saved Colorbar and Contour Level Configurations
#---------------------------------------------------------------------------------------------------

# Defaults
cmap = 'plasma'
lvls = np.arange(-75, 75.1, 5)
diff_lvls = np.arange(-17, 17.1, 2)
extend = 'both'

# Specific settings for various fields
if (field == 'gridlat_0'):
    lvls = np.arange(20, 60, 5)
    diff_lvls = np.arange(-17, 17.1, 2)*0.0001
    prslev = np.nan
elif (field == 'gridlon_0'):
    lvls = np.arange(-130, -60, 5)
    diff_lvls = np.arange(-17, 17.1, 2)*0.0001
    prslev = np.nan
elif (field in ['wspd', 'UGRD_P0_L100_GLC0', 'VGRD_P0_L100_GLC0']) and (prslev == 25000):
    cmap = 'PiYG'
    lvls = np.arange(-75, 75.1, 5)
    diff_lvls = np.arange(-17, 17.1, 2)
elif (field == 'wdir') and (prslev == 25000):
    lvls = np.arange(-90, 271, 20)
    diff_lvls = np.arange(-17, 17.1, 2)*2.5


#---------------------------------------------------------------------------------------------------
# Plot Data
#---------------------------------------------------------------------------------------------------

# Create plots
mnames = list(model_info.keys())
for fname1, fname2 in zip(model_info[mnames[0]]['files'], model_info[mnames[1]]['files']):
    print('Plotting for %s' % fname1.split('/')[-1])
    fig = plt.figure(figsize=(8, 11))
    ds = []
    zind = []
    out = []
    for j, (f, n) in enumerate(zip([fname1, fname2], mnames)):
        ds.append(xr.open_dataset(f, engine='pynio'))

        # Compute wind speed and direction
        if field in ['wspd', 'wdir']:
            ds[-1]['wspd'] = np.sqrt(ds[-1]['UGRD_P0_L100_GLC0']**2 + ds[-1]['VGRD_P0_L100_GLC0']**2)
            ds[-1]['wspd'].attrs = ds[-1]['UGRD_P0_L100_GLC0'].attrs
            ds[-1]['wspd'].attrs['long_name'] = 'wind speed'
            ds[-1]['wdir'] = 90. - np.rad2deg(np.arctan2(-ds[-1]['VGRD_P0_L100_GLC0'], -ds[-1]['UGRD_P0_L100_GLC0']))
            ds[-1]['wdir'].attrs = ds[-1]['UGRD_P0_L100_GLC0'].attrs
            ds[-1]['wdir'].attrs['long_name'] = 'wind direction'
            ds[-1]['wdir'].attrs['units'] = 'deg'

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
                                               'levels':diff_lvls},
                                      ingest_kw={'zind':zind, 
                                                 'indices':[model_info[mnames[0]]['indices'], 
                                                            model_info[mnames[1]]['indices']]})
    diff.config_ax(grid=False)
    if lat_lim != None:
        diff.set_lim(lat_lim[0], lat_lim[1], lon_lim[0], lon_lim[1])
    diff.ax_title(txt='difference ({n1} $-$ {n2})'.format(n1=mnames[0], n2=mnames[1]), 
                  size=14)

    plt.subplots_adjust(left=0.02, bottom=0.04, right=0.98, top=0.97, hspace=0.1)
    
    # Create output file name
    try:
        fhr = ds[0][field].attrs['forecast_time'][0]
        init_time = dt.datetime.strptime(ds[0][field].attrs['initial_time'], '%m/%d/%Y (%H:%M)')
    except KeyError:
        print('Warning: Missing fhr and init_time')
        fhr = 0
        init_time = dt.datetime(2000, 1, 1, 0, 0)
    if np.isnan(prslev):
        plt.savefig('{n1}vs{n2}_{tag}_{t}_f{fhr:03d}_{v}.png'.format(n1=mnames[0], n2=mnames[1], 
                    tag=save_tag, t=init_time.strftime('%Y%m%d%H'), fhr=fhr, v=field.split('_')[0])) 
    else:
        plt.savefig('{n1}vs{n2}_{tag}_{t}_f{fhr:03d}_{p}mb_{v}.png'.format(n1=mnames[0], n2=mnames[1], 
                    tag=save_tag, t=init_time.strftime('%Y%m%d%H'), fhr=fhr, v=field.split('_')[0],
                    p=prslev/1e2)) 
    plt.close()


"""
End plot_hcrsxn_diffs.py
""" 
