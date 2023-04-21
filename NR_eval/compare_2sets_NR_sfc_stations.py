"""
Compare Two Sets of Simulated Surface Station Observations

Real surface station obs come from the Iowa Environmental Mesonet database:
https://mesonet.agron.iastate.edu/request/download.phtml

This script uses "perfect" station matching. This means that the NR output WAS intentionally
interpolated to the real surface stations used for comparison here.

shawn.s.murdzek@noaa.gov
Date Created: 11 April 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import bufr
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
import metpy.calc as mc
from metpy.units import units


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Parameters for real obs
station_ids = ['ABR', 'ALB', 'BHM', 'CRP', 'DTW', 'GJT', 'MIA', 'OAK', 'SLE', 'TUS']
station_ids = ['OAK', 'TUS']

# Parameters for fake obs
# analysis_times are dt.timedelta objects relative to 0000
fake_obs_dir1 = '/work2/noaa/wrfruc/murdzek/nature_run_spring/sfc_stat_obs_csv/perfect'
fake_obs_dir2 = '/work2/noaa/wrfruc/murdzek/nature_run_spring/sfc_stat_obs_csv/old_obs/perfect'
analysis_days = [dt.datetime(2022, 4, 29) + dt.timedelta(days=i) for i in range(8)]
analysis_times = [dt.timedelta(hours=i) for i in range(24)]

# Title
title = 'new $-$ old'

# Maximum time allowed between analysis_times and either the real or fake ob (sec)
max_time_allowed = 450.

# Output file name (include %s placeholder for station ID)
out_fname = '%s_sfc_station_compare_sim_obs.png'


#---------------------------------------------------------------------------------------------------
# Compare NR to Surface Station Obs
#---------------------------------------------------------------------------------------------------

# Variables to extract
fake_varnames = ['TOB', 'QOB', 'POB', 'UOB', 'VOB', 'WSPD', 'WDIR']
fake_varnames_noplot = []

# Extract some values that will be used a lot
ndays = len(analysis_days)
ntimes = len(analysis_times)
analysis_year = analysis_days[0].year

# Extract fake observations
print()
fake_stations = {}
for ob_dir in [fake_obs_dir1, fake_obs_dir2]:
    fake_stations[ob_dir] = {}
    for ID in station_ids:
        fake_stations[ob_dir][ID] = {}
        for v in fake_varnames + fake_varnames_noplot:
            fake_stations[ob_dir][ID][v] = np.ones([ndays, ntimes]) * np.nan
    for j, d in enumerate(analysis_days):
        for k, t in enumerate(analysis_times):
            time = d + t
            print(time.strftime('%Y%m%d %H:%M'))
            try:
                full_bufr_csv = bufr.bufrCSV('%s/%s.sfc.fake.prepbufr.csv' % 
                                             (ob_dir, time.strftime('%Y%m%d%H%M')), use_all_col=True)
                full_bufr_csv.df = bufr.compute_wspd_wdir(full_bufr_csv.df)
            except FileNotFoundError:
                print('file not found, continuing to next time')
                continue
        
            # Remove rows that do not have thermodynamic AND kinematic data
            bufr_csv_no_nan = full_bufr_csv.df.loc[np.logical_not(np.isnan(full_bufr_csv.df['TOB']) |
                                                                  np.isnan(full_bufr_csv.df['UOB']))].copy()

            for ID in station_ids:
                red_csv = bufr_csv_no_nan.loc[bufr_csv_no_nan['SID'] == ID].copy()
                if len(red_csv) == 0:
                    # No obs, switching to next station
                    print('no simulated obs for %s' % ID)
                    continue
                red_csv.reset_index(inplace=True, drop=True)
                idx = np.argmin(np.abs(red_csv['DHR']))
                if (3600*np.abs(red_csv['DHR'].loc[idx])) < max_time_allowed:
                    for v in fake_varnames + fake_varnames_noplot:
                        fake_stations[ob_dir][ID][v][j, k] = red_csv.loc[idx, v]

# Plot results
plot_hr = np.array([t.total_seconds() / 3600. for t in analysis_times])
for ID in station_ids:
    ncols = 4
    fig, axes = plt.subplots(nrows=2, ncols=ncols, figsize=(12, 8))
    plt.subplots_adjust(left=0.07, bottom=0.08, right=0.99, top=0.9, wspace=0.48)
    for j, v in enumerate(fake_varnames):
        ax = axes[int(j/ncols), j%ncols]

        for k in range(ndays):
            ax.plot(plot_hr, fake_stations[fake_obs_dir1][ID][v][k, :] - 
                             fake_stations[fake_obs_dir2][ID][v][k, :], 'k-', lw=1)

        ax.grid()
        ax.axhline(0, color='b', lw=1)
        ax.set_ylabel('%s (%s)' % (v, full_bufr_csv.meta[v]['units']), size=14)
        ax.set_xlim([plot_hr.min(), plot_hr.max()])
    for j in range(ncols):
        if j != 3:
            axes[-1, j].set_xlabel('hour', size=14)
        axes[0, 3].set_xlabel('hour', size=14)

    plt.suptitle('%s: %s' % (title, ID), size=18)
    plt.savefig(out_fname % ID)
    plt.close()


"""
End compare_2sets_NR_sfc_stations.py
"""
