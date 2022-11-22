"""
Create Synthetic UAS Observations from WRF Nature Run Output

UAS observations are written to a preexisting prepbufr CSV output file

shawn.s.murdzek@noaa.gov
Date Created: 22 November 2022
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import numpy as np
import datetime as dt

import create_ob_utils as cou 


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Directory containing wrfred output
wrf_dir = '/lfs4/BMC/wrfruc/murdzek/nature_run_3km/WRF/run/'

# File containing UAS ob locations
uas_locs = '/lfs4/BMC/wrfruc/murdzek/nature_run_3km/synthetic_obs/uas_locs.csv' 

# Output directory for synthetic prepbufr CSV output
fake_bufr_dir = 'lfs4/BMC/wrfruc/murdzek/nature_run_3km/synthetic_obs/'

# Start and end times for prepbufrs
tstart = dt.datetime(2021, 7, 24, 23)
tend = dt.datetime(2021, 7, 25, 0)

# Parameters for creating UAS obs
# Units:
#     ascent_rate:   m/s
#     sampling_rate: s
#     max_z:         m
#     response:      s
#     temperature:   K
#     pressure:      hPa
#     rel_humidity:  %
#     wind_spd:      m/s
#     wind_dir:      deg
ascent_rate = 3.
sampling_rate = 5.
max_z = 2000.
fields = {'temperature':  {'error':    0.5,
                           'response': 2.},
          'pressure':     {'error':    1.5,
                           'response': 0.},
          'rel_humidity': {'error':    3.,
                           'response': 5.},
          'wind_spd':     {'error':    1.,
                           'response': 5.},
          'wind_dir':     {'error':    5.,
                           'response': 5.}}


#---------------------------------------------------------------------------------------------------
# Perform Interpolation
#---------------------------------------------------------------------------------------------------

ds = xr.open_dataset(wrf_fname)


"""
End create_obs.py
"""
