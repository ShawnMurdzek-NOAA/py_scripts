"""
Check How Many Data Points from the Surface Station Obs are Missing

Real surface station obs come from the Iowa Environmental Mesonet database:
https://mesonet.agron.iastate.edu/request/download.phtml

shawn.s.murdzek@noaa.gov
Date Created: 3 October 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
import pickle


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Parameters for real obs
real_obs_dir = './spring/'
with open('station_list.txt', 'r') as fptr:
    station_ids = fptr.readlines()
    for i in range(len(station_ids)):
        station_ids[i] = station_ids[i].strip()
years = np.arange(1993, 2023)
startdate = '04290000'
enddate = '05070000'
#startdate = '02010000'
#enddate = '02080000'

# Parameters for fake obs
# analysis_times are dt.timedelta objects relative to 0000
analysis_days = [dt.datetime(2022, 4, 29) + dt.timedelta(days=i) for i in range(8)]
#analysis_days = [dt.datetime(2022, 2, 1) + dt.timedelta(days=i) for i in range(7)]
analysis_times = [dt.timedelta(hours=i) for i in range(24)]

# Maximum time allowed between analysis_times and either the real or fake ob (sec)
max_time_allowed = 450.

# Maximum fraction of analysis times that are allowed to be NaN
nan_thres = 0.2


#---------------------------------------------------------------------------------------------------
# Check Surface Station Ob Files
#---------------------------------------------------------------------------------------------------

# Extract some values that will be used a lot
ndays = len(analysis_days)
ntimes = len(analysis_times)
analysis_year = analysis_days[0].year
min_hr = -max_time_allowed / 3600.
max_hr = max_time_allowed / 3600.

# Loop over each station
real_varnames = ['tmpf', 'dwpf', 'drct', 'sknt', 'alti']
for ID in station_ids:
    print()
    print('Extracting real obs at %s' % ID)
    n_nan = {}
    for v in real_varnames:
        n_nan[v] = 0
    for j, yr in enumerate(years):
        tmp_df = pd.read_csv('%s/%s_%d%s_%d%s.txt' % (real_obs_dir, ID, yr, startdate, yr, enddate), 
                             skiprows=5)

        if len(tmp_df) == 0:
            print('No data for %s %d (entire file is empty)' % (ID, yr))
            for v in real_varnames:
                n_nan[v] = n_nan[v] + (ntimes * ndays)
            continue

        # Only retain times with thermodynamic and kinematic obs (this eliminates special obs for 
        # gusts, etc.)
        ss_df = tmp_df.loc[(tmp_df['tmpf'] != 'M') & (tmp_df['drct'] != 'M')].copy()
        ss_df.reset_index(inplace=True, drop=True)

        station_times = [dt.datetime.strptime(s, '%Y-%m-%d %H:%M').replace(year=analysis_year) 
                         for s in ss_df['valid'].values]
        station_times_tot = np.array([(s - dt.datetime(analysis_year, 1, 1)).total_seconds() 
                                      for s in station_times])
      
        for k, d in enumerate(analysis_days):
            for l, time in enumerate(analysis_times):
                diff = np.abs(station_times_tot - ((d + time) - dt.datetime(analysis_year, 1, 1)).total_seconds())
                if len(diff) == 0:
                    for v in ['tmpf', 'dwpf', 'drct', 'alti']:
                        n_nan[v] = n_nan[v] + 1
                    continue
                idx = np.argmin(diff)
                for v in ['tmpf', 'dwpf', 'drct', 'alti']:
                    if (diff.min() > max_time_allowed) or (ss_df.loc[idx, v] == 'M'):
                        n_nan[v] = n_nan[v] + 1

    for v in real_varnames:
        if (n_nan[v] / (ndays*ntimes*len(years))) > nan_thres:
            print('> %s data is NaN for %s %s' % (nan_thres, v, ID))



"""
End check_sfc_station_data.py
"""
