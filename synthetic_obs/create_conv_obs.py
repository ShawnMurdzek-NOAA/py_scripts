"""
Create Synthetic Conventional Observations from WRF Nature Run Output

Conventional obs are written into a new CSV file. UAS obs can then be added later.

Prepbufr CSV files for real observations can be created using GSI-utils/bin/prepbufr_decode_csv.x

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

import create_ob_utils as cou 
import bufr


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Directory containing wrfred output
wrf_dir = '/lfs4/BMC/wrfruc/murdzek/nature_run_3km/WRF/run/'

# Directory containing real prepbufr CSV output
bufr_dir = '/lfs4/BMC/wrfruc/murdzek/sample_real_obs/obs_rap/'

# File containing conventional observation error characteristics
# This file also contains the synthetic obs that are actually created
error_fname = '/lfs4/BMC/wrfruc/murdzek/sample_real_obs/obs_errors.json'

# Output directory for synthetic prepbufr CSV output
fake_bufr_dir = 'lfs4/BMC/wrfruc/murdzek/nature_run_3km/synthetic_obs/'

# Start and end times for prepbufrs. Step is in min
tstart = dt.datetime(2021, 7, 24, 23)
tend = dt.datetime(2021, 7, 25, 0)
step = 60


#---------------------------------------------------------------------------------------------------
# Create Synthetic Observations
#---------------------------------------------------------------------------------------------------

# Open up obs error file
errors = json.load(open(error_fname, 'r'))

ntimes = int((tend - tstart) / dt.timedelta(minutes=step) + 1)
for i in range(ntimes):
    t = tstart + dt.timedelta(minutes=(i*step))
    bufr_fname = bufr_dir + t.strftime('%Y%m%d%H%M.rap.prepbufr.csv')
    bufr_csv = bufr.bufrCSV(bufr_fname)

    # Only keep platforms if we are creating synthetic obs for them
    obs = bufr_csv.df['subset'].unique()
    for s in obs:
        if s not in errors.keys():
            bufr_csv.df.drop(index=np.where(bufr_csv.df['subset'] == s)[0], inplace=True)
            bufr_csv.df.reset_index(drop=True, inplace=True)


"""
End create_conv_obs.py
"""
