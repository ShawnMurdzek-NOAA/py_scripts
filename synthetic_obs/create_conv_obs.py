"""
Create Synthetic Conventional Observations from WRF Nature Run Output

Conventional obs are written into a new CSV file. UAS obs can then be added later.

Prepbufr CSV files for real observations can be created using GSI-utils/bin/prepbufr_decode_csv.x

This code uses UPP output on WRF native levels (wrfnat). This is likely more efficient than using 
the wrfred files, which can take some time (and a lot of memory) to create.

This code requires the Dask module, which is not part of the adb_graphics environment, so it must
be installed separately and the location must be added to the PYTHONPATH variable

shawn.s.murdzek@noaa.gov
Date Created: 22 November 2022
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
wrf_dir = '/lfs4/BMC/wrfruc/murdzek/nature_run_3km/WRF/run/'

# Directory containing real prepbufr CSV output
bufr_dir = '/lfs4/BMC/wrfruc/murdzek/sample_real_obs/obs_rap/'

# File containing conventional observation error characteristics (use GSI errtable file)
error_fname = '/lfs4/BMC/wrfruc/murdzek/sample_real_obs/errtable.rrfs'

# Observation platforms to use (aka subsets, same ones used by BUFR)
ob_platforms = ['ADPUPA', 'AIRCAR', 'AIRCFT', 'PROFLR', 'ADPSFC', 'SFCSHP', 'MSONET', 'GPSIPW']

# Output directory for synthetic prepbufr CSV output
fake_bufr_dir = 'lfs4/BMC/wrfruc/murdzek/nature_run_3km/synthetic_obs/'

# Start and end times for prepbufrs. Step is in min
bufr_start = dt.datetime(2021, 7, 24, 23)
bufr_end = dt.datetime(2021, 7, 25, 0)
bufr_step = 120

# Start and end times for wrfnat UPP output. Step is in min
wrf_start = dt.datetime(2021, 7, 24, 23, 0)
wrf_end = dt.datetime(2021, 7, 25, 0, 0)
wrf_step = 15


#---------------------------------------------------------------------------------------------------
# Create Synthetic Observations
#---------------------------------------------------------------------------------------------------

# Open up obs error file
errors = cou.read_ob_errors(error_fname)

ntimes = int((bufr_end - bufr_start) / dt.timedelta(minutes=bufr_step) + 1)
for i in range(ntimes):
    t = bufr_start + dt.timedelta(minutes=(i*bufr_step))
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
    for hr in np.arange(hr_start, hr_end, wrf_step / 60.):
        wrf_t = t + dt.timedelta(hours=hr)
        wrf_ds[hr] = xr.open_dataset(wrf_dir + wrf_t.strftime('/wrfnat_hrconus_%Y%m%d%H%M.grib2'),
                                     engine='pynio')
 
    # Loop over each station ID (SID)
    # This approach is preferred over looping over each ob time b/c it will make it easier to add
    # dynamic errors (owing to nonzero response times) and instrument bias
    #for s in out_df['SID'].unique():
        
        # Interpolate horizontally

        # Interpolate vertically

        # Interpolate temporally (between the two datasets)

        # Add errors

    # Write output DataFrame to a CSV file
    #bufr.df_to_csv(out_df, fake_bufr_dir + t.strftime('/%Y%m%d%H%M.fake.prepbufr.csv'))


"""
End create_conv_obs.py
"""
