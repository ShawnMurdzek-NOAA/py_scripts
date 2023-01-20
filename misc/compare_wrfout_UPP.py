"""
Compare UPP and WRF Output Files

shawn.s.murdzek@noaa.gov
Date Created: 20 January 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# File names (use WRF output on native grid)
wrfout = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_10km_hrrr/WRF/run/wrfout_d01_2021-07-24_23_15_00'
upp = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_10km_hrrr/UPP/out/wrfnat_202107242315.grib2'

# Fields to compare (the _P0_L105_GLC0 suffix is added to the UPP fields later)
wrf_fields = ['QVAPOR', 'QCLOUD', 'QRAIN',  'QICE',  'QSNOW',  'QGRAUP',  
                        'QNDROP',                    'QNSNOW', 'QNGRAUPEL']
upp_fields = ['SPFH',   'CLWMR',  'RWMR',   'ICMR',  'SNMR',   'GRLE',     
                        'NCONCD',                    'SPNCS',  'SPNCG']

#wrf_fields = ['QRAIN']
#upp_fields = ['RWMR']

# Vertical level for difference plots
zlvl = 8

# Output filename for difference plots. Include %s placeholder for field name
save_fname = 'upp_wrf_diffs_%s.png'


#---------------------------------------------------------------------------------------------------
# Compare UPP and WRF Fields
#---------------------------------------------------------------------------------------------------

wrf_ds = xr.open_dataset(wrfout)
upp_ds = xr.open_dataset(upp, engine='pynio')

# Extract lat and lon arrays for plotting
wrflat = wrf_ds['XLAT'][0, :, :].values
wrflon = wrf_ds['XLONG'][0, :, :].values
upplat = upp_ds['gridlat_0'].values
upplon = upp_ds['gridlon_0'].values

for wf, up in zip(wrf_fields, upp_fields):
    wrf_data = np.squeeze(wrf_ds[wf].values)
    try:
        upp_data = upp_ds[up + '_P0_L105_GLC0'].values
    except KeyError:
        upp_data = upp_ds[up].values
    diff = wrf_data - upp_data
    diff[np.isclose(wrf_data, 0)] = np.nan
    diff[np.isclose(upp_data, 0)] = np.nan
    rdiff = diff / wrf_data
    print()
    print(wf)
    print('UPP maxval        = %.4e' % upp_data.max())
    print('WRF maxval        = %.4e' % wrf_data.max())
    print('max absolute diff = %.4e' % np.nanmax(np.abs(diff)))
    print('RMSD              = %.4e' % np.sqrt(np.nanmean(diff*diff)))
    print('max rel diff      = %.4f' % np.nanmax(np.abs(rdiff)))
    print('RMSD (rel diff)   = %.4f' % np.sqrt(np.nanmean(rdiff*rdiff)))

    # Plot differences
    fig = plt.figure(figsize=(8, 10))
    for j, (data, lat, lon, ttl) in enumerate(zip([wrf_data, upp_data], [wrflat, upplat], 
                                                  [wrflon, upplon], ['WRF', 'UPP'])):
        ax = fig.add_subplot(3, 1, j+1, projection=ccrs.PlateCarree())
        cax = ax.contourf(lon, lat, data[zlvl, :, :])
        cbar = plt.colorbar(cax, ax=ax, orientation='vertical')
        cbar.set_label(wf)
        ax.set_title(ttl, size=14)
        ax.coastlines('50m')

    ax = fig.add_subplot(3, 1, 3, projection=ccrs.PlateCarree()) 
    mxd = np.nanmax(np.abs(diff[zlvl, :, :]))
    cax = ax.contourf(wrflon, wrflat, diff[zlvl, :, :], np.linspace(-1*mxd, mxd, 12), 
                      cmap='bwr')
    cbar = plt.colorbar(cax, ax=ax, orientation='vertical')
    cbar.set_label('diff (WRF $-$ UPP)', size=14)
    ax.coastlines('50m')

    plt.suptitle('z lvl = %d' % zlvl, size=16)
    plt.savefig(save_fname % wf)
    plt.close()


"""
End compare_wrfout_UPP.py
"""
