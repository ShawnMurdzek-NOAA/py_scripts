"""
Compute 6-hr Precipitation Totals from 1-hr Precipitation Totals

This script is designed to operate on netCDF files containing 1-hr precip amounts that are output by 
NR_eval/extract_wrfnat_fields.py. 

It is assumed that these 1-hr precip amounts are in 1-hr increments.

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import datetime as dt
import numpy as np


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

start_date = dt.datetime(2022, 2, 1)
in_files = ['/work2/noaa/wrfruc/murdzek/nature_run_winter/UPP/%s/precip1hr_all_%s.nc' % 
            ((start_date + dt.timedelta(days=i)).strftime('%Y%m%d'), 
             (start_date + dt.timedelta(days=i)).strftime('%Y%m%d')) for i in range(7)]

# Files containing 1-hr precip amounts from the previous day. Set to None if a file does not exist
prev_files = [None] + in_files[:-1]


#---------------------------------------------------------------------------------------------------
# Compute 6-hr Precip Amounts
#---------------------------------------------------------------------------------------------------

field = 'APCP_P8_L1_GLC0_acc'

for f, pf in zip(in_files, prev_files):
    print('processing %s' % f)
    new_ds = xr.open_dataset(f).copy()
    new_ds[field][:, :, :] = np.nan
    if pf != None:
        ds = xr.open_mfdataset([f, pf])
    else:
        ds = xr.open_dataset(f)
    for j, t in enumerate(new_ds['time'].values):
        ind = []
        for k in range(6):
            try: 
                ind.append(np.where(ds['time'] == (t - np.timedelta64(1, 'h')))[0][0])
            except IndexError:
                # Not enough data to compute 6-hr precip amounts
                break
        if len(ind) == 6:
            new_ds[field][j, :, :] = np.sum(ds[field][ind, :, :], axis=0)
    new_ds[field].attrs['statistical_process_duration'] = '6 hrs'
    new_ds.to_netcdf('%s/precip6hr_%s' % ('/'.join(f.split('/')[:-1]), f.split('/')[-1][-11:]))
        

"""
End compute_precip6hr.py
"""
