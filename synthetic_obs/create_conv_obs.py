"""
Create Synthetic Conventional Observations from WRF Nature Run Output

Conventional obs are written into a new CSV file. UAS obs can then be added later.

Prepbufr CSV files for real observations can be created using GSI-utils/bin/prepbufr_decode_csv.x.
To use this utility, do the following:

1. Copy the prepbufr file you wish to decode to GSI-utils/bin and rename it 'prepbufr'
2. Run ./prepbufr_decode_csv.x
3. Move resulting CSV file to wherever you choose

This code uses UPP output on WRF native levels (wrfnat). This is likely more efficient than using 
the wrfred files, which can take some time (and a lot of memory) to create.

shawn.s.murdzek@noaa.gov
Date Created: 22 November 2022
Environment: adb_graphics
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import numpy as np
import datetime as dt
import json
import math

import create_ob_utils as cou 
import bufr


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Directory containing wrfnat output from UPP
wrf_dir = '/scratch1/BMC/wrfruc/murdzek/nature_run_spring_v2/output/202204291200/UPP/'

# Directory containing real prepbufr CSV output
bufr_dir = '/scratch1/BMC/wrfruc/murdzek/sample_real_obs/obs_rap/'

# File containing conventional observation error characteristics (use GSI errtable file)
error_fname = '/scratch1/BMC/wrfruc/murdzek/sample_real_obs/errtable.rrfs'

# Observation platforms to use (aka subsets, same ones used by BUFR)
ob_platforms = ['ADPUPA', 'AIRCAR', 'AIRCFT', 'PROFLR', 'ADPSFC', 'SFCSHP', 'MSONET', 'GPSIPW']

# Output directory for synthetic prepbufr CSV output
fake_bufr_dir = '/scratch1/wrfruc/murdzek/nature_run_spring_v2/output/202204291200/synthetic_obs/'

# Start and end times for prepbufrs. Step is in min
bufr_start = dt.datetime(2022, 4, 29, 12)
bufr_end = dt.datetime(2022, 4, 29, 13)
bufr_step = 120

# Start and end times for wrfnat UPP output. Step is in min
wrf_start = dt.datetime(2022, 4, 29, 12, 0)
wrf_end = dt.datetime(2022, 4, 29, 13, 0)
wrf_step = 15


#---------------------------------------------------------------------------------------------------
# Create Synthetic Observations
#---------------------------------------------------------------------------------------------------

# Open up obs error file
errors = cou.read_ob_errors(error_fname)

ntimes = int((bufr_end - bufr_start) / dt.timedelta(minutes=bufr_step) + 1)
for i in range(ntimes):
    t = bufr_start + dt.timedelta(minutes=(i*bufr_step))
    print()
    print('t = %s' % t.strftime('%Y-%m-%d %H%M%S'))
    bufr_fname = bufr_dir + t.strftime('/%Y%m%d%H%M.rap.prepbufr.csv')
    bufr_csv = bufr.bufrCSV(bufr_fname)

    # Only keep platforms if we are creating synthetic obs for them
    obs = bufr_csv.df['subset'].unique()
    for s in obs:
        if s not in ob_platforms:
            bufr_csv.df.drop(index=np.where(bufr_csv.df['subset'] == s)[0], inplace=True)
            bufr_csv.df.reset_index(drop=True, inplace=True)

    # Remove obs outside of the range of wrfnat files
    hr_min = ((wrf_start - t).days * 24) + ((wrf_start - t).seconds / 3600.)
    hr_max = ((wrf_end - t).days * 24) + ((wrf_end - t).seconds / 3600.)
    bufr_csv.df.drop(index=np.where(np.logical_or(bufr_csv.df['DHR'] < hr_min, 
                                                  bufr_csv.df['DHR'] > hr_max))[0], inplace=True)
    bufr_csv.df.reset_index(drop=True, inplace=True)

    # Create output DataFrame
    out_df = bufr_csv.df.copy()

    # Open wrfnat files
    # Ideally this would be done using xr.open_mfdataset(), but the wrfnat files lack a time 
    # dimension, which makes Xarray unable to combine them (and using concatenate requires a lot
    # of memory, so that is not a viable option)
    hr_start = math.floor(out_df['DHR'].min()*4) / 4
    hr_end = math.ceil(out_df['DHR'].max()*4) / 4 + (wrf_step / 60.)
    wrf_ds = {}
    wrf_hr = []
    for hr in np.arange(hr_start, hr_end, wrf_step / 60.):
        wrf_t = t + dt.timedelta(hours=hr)
        print(wrf_dir + wrf_t.strftime('wrfnat_%Y%m%d%H%M.grib2'))
        wrf_ds[hr] = xr.open_dataset(wrf_dir + wrf_t.strftime('wrfnat_%Y%m%d%H%M.grib2'),
                                     engine='pynio')
        wrf_hr.append(hr)
 
    # Loop over each observation
    wrf_hr = np.array(wrf_hr)
    for j in range(len(out_df))[:3]:
        
        subset = out_df.iloc[j]

        # Determine WRF hour right before observation and weight for temporal interpolation
        ihr = np.where((wrf_hr - subset['DHR']) <= 0)[0][-1]
        twgt = wrf_hr[ihr] / (wrf_hr[ihr] - wrf_hr[ihr+1]) 

        # Interpolate horizontally
        xi, yi = cou.wrf_coords(subset['YOB'], subset['XOB'] - 360., wrf_ds[wrf_hr[ihr]], ftype='UPP')        

        # Interpolate vertically
        # Should probably add a conditional here to not interpolate vertically for 2-m and 10-m obs.
        # Also need to add code to interpolate surface pressure (PRSS) and precipitable water (PWO)
        zi = np.zeros(2)
        for k in [ihr, ihr+1]:
            z1d = wrf_ds[wrf_hr[k]]['HGT_P0_L105_GLC0'].interp(ygrid_0=yi, xgrid_0=xi)
            z0i = np.where((z1d - subset['ZOB']) <= 0)[0][-1]
            zi[k] = z0i + (subset['ZOB'] - z1d[z0i]) / (z1d[z0i+1] - z1d[z0i]) 

        # Interpolate temporally
        # Double check that these units are right!
        obs_name = ['POB', 'QOB', 'TOB', 'UOB', 'VOB']
        wrf_name = ['PRES_P0_L105_GLC0']
        for o, m in zip(obs_name, wrf_name):
            if not np.isnan(subset[o]):
                out_df['POB'].iloc[j] =  (wgt * wrf_ds[wrf_hr[ihr]][m].interp(lv_HYBL1=zi, ygrid_0=yi, xgrid_0=xi) +
                                         (1-wgt) * wrf_ds[wrf_hr[ihr+1]][m].interp(lv_HYBL1=zi, ygrid_0=yi, xgrid_0=xi))

        # Add errors

    # Write output DataFrame to a CSV file
    #bufr.df_to_csv(out_df, fake_bufr_dir + t.strftime('/%Y%m%d%H%M.fake.prepbufr.csv'))


"""
End create_conv_obs.py
"""
