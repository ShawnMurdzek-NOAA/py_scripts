"""
Perform Comparisons Between Real and Simulated Surface Station Observations

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
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
import metpy.constants as const
import metpy.calc as mc
from metpy.units import units
import scipy.stats as ss

import pyDA_utils.bufr as bufr


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Parameters for real obs
real_obs_dir = '/work2/noaa/wrfruc/murdzek/real_obs/sfc_stations/winter'
station_ids = ['ABR', 'ALB', 'BHM', 'CRP', 'DTW', 'GJT', 'MIA', 'OAK', 'SLE', 'TUS']
years = np.arange(1993, 2023)
#startdate = '04290000'
#enddate = '05070000'
startdate = '02010000'
enddate = '02080000'

# Parameters for fake obs
# analysis_times are dt.timedelta objects relative to 0000
fake_obs_dir = '/work2/noaa/wrfruc/murdzek/nature_run_winter/sfc_stat_obs_csv/perfect'
#analysis_days = [dt.datetime(2022, 4, 29) + dt.timedelta(days=i) for i in range(8)]
analysis_days = [dt.datetime(2022, 2, 1) + dt.timedelta(days=i) for i in range(7)]
analysis_times = [dt.timedelta(hours=i) for i in range(24)]

# Maximum time allowed between analysis_times and either the real or fake ob (sec)
max_time_allowed = 450.

# Option to plot ceiling obs as a binary histogram ('binary') or CDF of ceiling height ('hgt')
ceil_plot = 'hgt'

# Type of z-score to use for z-score plots ('regular' or 'modified')
# Also select cutoff value that denotes outliers
zscore = 'modified'
zcutoff = 3

# Output file name (include %s placeholder for station ID)
out_fname = '%s_sfc_station_compare_perfect.png'


#---------------------------------------------------------------------------------------------------
# Compare NR to Surface Station Obs
#---------------------------------------------------------------------------------------------------

# Variables to extract
real_varnames = ['lon', 'lat', 'tmpf', 'dwpf', 'drct', 'sknt', 'alti', 'vsby', 'elevation']
fake_varnames = ['TOB', 'QOB', 'POB', 'UOB', 'VOB', 'WSPD', 'WDIR']
fake_varnames_noplot = ['ELV', 'ceil']

# Extract some values that will be used a lot
ndays = len(analysis_days)
ntimes = len(analysis_times)
analysis_year = analysis_days[0].year
min_hr = -max_time_allowed / 3600.
max_hr = max_time_allowed / 3600.

# Extract real surface station obs first
real_stations = {}
for ID in station_ids:
    print('Extracting real obs at %s' % ID)
    real_stations[ID] = {}
    for v in real_varnames:
        real_stations[ID][v] = np.ones([len(years) * ndays, ntimes]) * np.nan
    real_stations[ID]['ceil'] = np.ones([len(years), ntimes * ndays]) * np.nan
    real_stations[ID]['n_skyl'] = np.zeros(len(years))
    real_stations[ID]['frac_ceil'] = np.zeros(len(years))
    for j, yr in enumerate(years):
        tmp_df = pd.read_csv('%s/%s_%d%s_%d%s.txt' % (real_obs_dir, ID, yr, startdate, yr, enddate), 
                             skiprows=5)

        # Only retain times with thermodynamic and kinematic obs (this eliminates special obs for 
        # gusts, etc.)
        ss_df = tmp_df.loc[(tmp_df['tmpf'] != 'M') & (tmp_df['drct'] != 'M')].copy()
        ss_df.reset_index(inplace=True, drop=True)

        station_times = [dt.datetime.strptime(s, '%Y-%m-%d %H:%M').replace(year=analysis_year) 
                         for s in ss_df['valid'].values]
        station_times_tot = np.array([(s - dt.datetime(analysis_year, 1, 1)).total_seconds() 
                                      for s in station_times])

        # Temporary variables for determining fraction of time with a ceiling
        ceil_obs = 0
        tot_obs = 0

        for k, d in enumerate(analysis_days):
            for l, time in enumerate(analysis_times):
                diff = np.abs(station_times_tot - ((d + time) - dt.datetime(analysis_year, 1, 1)).total_seconds())
                idx = np.argmin(diff)
                for v in real_varnames:
                    if (diff.min() > max_time_allowed) or (ss_df.loc[idx, v] == 'M'):
                        real_stations[ID][v][(j*ndays) + k, l] = np.nan
                    else:
                        real_stations[ID][v][(j*ndays) + k, l] = ss_df.loc[idx, v]

                # Determine if there is a ceiling
                skyc = ss_df.loc[idx, ['skyc%d' % n for n in range(1, 5)]].values
                skyl = ss_df.loc[idx, ['skyl%d' % n for n in range(1, 5)]].values
                skyl[skyl == 'M'] = np.nan
                skyl = np.float64(skyl)
                ceil_idx = np.where((skyc == 'OVC') | (skyc == 'BKN'))[0]
                if len(ceil_idx) > 0:
                    ceil_obs = ceil_obs + 1
                    if not np.all(np.isnan(skyl[ceil_idx])):
                        real_stations[ID]['n_skyl'][j] = real_stations[ID]['n_skyl'][j] + 1
                        real_stations[ID]['ceil'][j, (k*ntimes) + l] = np.nanmin(skyl[ceil_idx])
                if not np.all(skyc == 'M'):
                    tot_obs = tot_obs + 1
                    if len(ceil_idx) == 0:
                        real_stations[ID]['n_skyl'][j] = real_stations[ID]['n_skyl'][j] + 1

        real_stations[ID]['frac_ceil'][j] = ceil_obs / tot_obs

# Convert real obs to same units/variables as fake obs
# Note that the surface stations report the altimeter setting rather than station-level pressure.
# The difference between these two is discussed here: https://www.weather.gov/bou/pressure_definitions
for ID in station_ids:
    real_stations[ID]['POB'] = mc.altimeter_to_station_pressure(real_stations[ID]['alti'] * units.inHg, 
                                                                real_stations[ID]['elevation'] * units.m).to('mbar').magnitude
    real_stations[ID]['TOB'] = (real_stations[ID]['tmpf'] * units.degF).to('degC').magnitude       
    real_stations[ID]['QOB'] = mc.specific_humidity_from_dewpoint(real_stations[ID]['POB'] * units.mbar,
                                                                  real_stations[ID]['dwpf'] * units.degF).magnitude * 1e6
    wnd_tmp = mc.wind_components(real_stations[ID]['sknt'] * units.kt, real_stations[ID]['drct'] * units.deg)
    real_stations[ID]['UOB'] = wnd_tmp[0].to(units.m / units.s).magnitude
    real_stations[ID]['VOB'] = wnd_tmp[1].to(units.m / units.s).magnitude
    real_stations[ID]['WSPD'] = (real_stations[ID]['sknt'] * units.kt).to(units.m / units.s).magnitude
    real_stations[ID]['WDIR'] = real_stations[ID]['drct']
    real_stations[ID]['lon'] = real_stations[ID]['lon'] + 360.
    real_stations[ID]['ceil'] = (real_stations[ID]['ceil'] * units.ft).to(units.km).magnitude    

# Extract fake observations
print()
fake_stations = {}
for ID in station_ids:
    fake_stations[ID] = {}
    for v in fake_varnames + fake_varnames_noplot:
        fake_stations[ID][v] = np.ones([ndays, ntimes]) * np.nan
for i, d in enumerate(analysis_days):
    for j, t in enumerate(analysis_times):
        time = d + t
        print(time.strftime('%Y%m%d %H:%M'))
        try:
            full_bufr_csv = bufr.bufrCSV('%s/%s.sfc.fake.prepbufr.csv' % 
                                         (fake_obs_dir, time.strftime('%Y%m%d%H%M')), use_all_col=True)
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
                    fake_stations[ID][v][i, j] = red_csv.loc[idx, v]

# Convert ceiling from gpm to km and compute binary ceiling
for ID in station_ids:
    fake_stations[ID]['ceil'] = (mc.geopotential_to_height(fake_stations[ID]['ceil'] * units.m * const.g).to('km').magnitude -
                                 mc.geopotential_to_height(fake_stations[ID]['ELV'] * units.m * const.g).to('km').magnitude)
    fake_stations[ID]['bin_ceil'] = np.float64(~np.isnan(fake_stations[ID]['ceil']))
    fake_stations[ID]['bin_ceil'][np.isnan(fake_stations[ID]['TOB'])] = np.nan
    fake_stations[ID]['frac_ceil'] = (np.nansum(fake_stations[ID]['bin_ceil']) / 
                                      np.count_nonzero(~np.isnan(fake_stations[ID]['bin_ceil'])))

# Compute z-scores for fake stations
# For modified z-scores explanation, see here: 
# https://medium.com/analytics-vidhya/anomaly-detection-by-modified-z-score-f8ad6be62bac
fake_stations_z = {}
all_zscores = {}
for v in fake_varnames:
    all_zscores[v] = np.zeros([ndays*len(station_ids), ntimes])
for i, ID in enumerate(station_ids):
    fake_stations_z[ID] = {}
    for v in fake_varnames:
        fake_stations_z[ID][v] = np.zeros([ndays, ntimes])
        for k in range(ndays):
            if zscore == 'regular':
                ctr = np.nanmean(real_stations[ID][v][k::ndays, :], axis=0)
                spd = np.nanstd(real_stations[ID][v][k::ndays, :], axis=0)
            elif zscore == 'modified':
                ctr = np.nanmedian(real_stations[ID][v][k::ndays, :], axis=0)
                spd = 1.4826 * ss.median_abs_deviation(real_stations[ID][v][k::ndays, :], axis=0, 
                                                       nan_policy='omit')
            spd[np.isclose(spd, 0)] = np.nan
            fake_stations_z[ID][v][k, :] = (fake_stations[ID][v][k, :] - ctr) / spd
            all_zscores[v][i*ndays+k, :] = fake_stations_z[ID][v][k, :]

# Plot results
ncols = 4
plot_hr = np.array([t.total_seconds() / 3600. for t in analysis_times])
fig_ztot, axes_ztot = plt.subplots(nrows=2, ncols=ncols, figsize=(12, 8))
plt.subplots_adjust(left=0.06, bottom=0.08, right=0.99, top=0.9, wspace=0.4)
for ID in station_ids:
    fig, axes = plt.subplots(nrows=2, ncols=ncols, figsize=(12, 8))
    plt.subplots_adjust(left=0.06, bottom=0.08, right=0.99, top=0.9, wspace=0.4)
    fig_z, axes_z = plt.subplots(nrows=2, ncols=ncols, figsize=(12, 8))
    plt.subplots_adjust(left=0.06, bottom=0.08, right=0.99, top=0.9, wspace=0.4)
    for j, v in enumerate(fake_varnames):
        ax = axes[int(j/ncols), j%ncols]
        ax_z = axes_z[int(j/ncols), j%ncols]

        # Plot fake station data
        for k in range(ndays):
            ax.plot(plot_hr, fake_stations[ID][v][k, :], 'k-', lw=0.75)
            ax_z.plot(plot_hr, fake_stations_z[ID][v][k, :], 'k-', lw=0.75)
    
        ax_z.axhline(zcutoff, c='r', ls='-', lw=1.5)
        ax_z.axhline(-zcutoff, c='r', ls='-', lw=1.5)

        # Plot real station data percentiles
        pcts = [0, 10, 25, 50, 75, 90, 100]
        var_percentiles = {}
        for p in pcts:
            var_percentiles[p] = np.zeros(ntimes)
            for l in range(ntimes):
                var_percentiles[p][l] = np.nanpercentile(real_stations[ID][v][:, l], p) 
        
        ax.plot(plot_hr, var_percentiles[50], 'r-', lw=2)
        ax.fill_between(plot_hr, var_percentiles[25], var_percentiles[75], color='r', alpha=0.4)
        ax.fill_between(plot_hr, var_percentiles[10], var_percentiles[90], color='r', alpha=0.2)
        ax.plot(plot_hr, var_percentiles[0], 'r-', lw=0.5)
        ax.plot(plot_hr, var_percentiles[100], 'r-', lw=0.5)

        ax.set_ylabel('%s (%s)' % (v, full_bufr_csv.meta[v]['units']), size=14)
        ax_z.set_ylabel('%s z-score' % v, size=14)
        for a in [ax, ax_z]:
            a.grid()
            a.set_xlim([plot_hr.min(), plot_hr.max()])
    for j in range(ncols):
        if j != 3:
            axes[-1, j].set_xlabel('hour', size=14)
            axes_z[-1, j].set_xlabel('hour', size=14)
        axes[0, 3].set_xlabel('hour', size=14)
        axes_z[0, 3].set_xlabel('hour', size=14)

    # Add fraction of time with a ceiling in the final subplot
    ax = axes[1, 3]
    if ceil_plot == 'binary':
        ax.hist(real_stations[ID]['frac_ceil'], color='red')
        ax.axvline(fake_stations[ID]['frac_ceil'], c='k', lw=2)
        ax.set_xlabel('fraction of time with ceiling', size=10)
        ax.set_ylabel('number of years', size=14)
        ax.set_xlim([0, 1])
    elif ceil_plot == 'hgt':
        bins = np.arange(0, 12.5, 0.5)
        ctrs = 0.5 * (bins[1:] + bins[:-1]) 
        
        real_cdf = np.zeros([len(years), len(ctrs)])
        for j in range(len(years)):
            tmp, _ = np.histogram(real_stations[ID]['ceil'][j, :], bins=bins)
            real_cdf[j, :] = (np.array([np.sum(tmp[:(n+1)]) for n in range(len(tmp))]) / 
                              real_stations[ID]['n_skyl'][j])
        pcts = [0, 10, 25, 50, 75, 90, 100]
        var_percentiles = {}
        for p in pcts:
            var_percentiles[p] = np.zeros(len(ctrs))
            for l in range(len(ctrs)):
                var_percentiles[p][l] = np.nanpercentile(real_cdf[:, l], p) 
        
        ax.plot(ctrs, var_percentiles[50], 'r-', lw=2)
        ax.fill_between(ctrs, var_percentiles[25], var_percentiles[75], color='r', alpha=0.4)
        ax.fill_between(ctrs, var_percentiles[10], var_percentiles[90], color='r', alpha=0.2)
        ax.plot(ctrs, var_percentiles[0], 'r-', lw=0.5)
        ax.plot(ctrs, var_percentiles[100], 'r-', lw=0.5)
        
        tmp, _ = np.histogram(fake_stations[ID]['ceil'], bins=bins)
        fake_cdf = (np.array([np.sum(tmp[:(n+1)]) for n in range(len(tmp))]) / 
                    np.count_nonzero(~np.isnan(fake_stations[ID]['bin_ceil'])))
        ax.plot(ctrs, fake_cdf, 'k-')
        ax.set_xlabel('cloud ceiling (km)', size=14)
        ax.set_ylabel('cumulative fraction of days', size=14)
        ax.set_xlim([0, 12])
        ax.set_ylim([0, 1])

    ax.grid()

    for f, tag, ttl in zip([fig, fig_z], [ID, '%s_z%s' % (ID, zscore)], 
                           [ID, '%s: %s z-score' % (ID, zscore)]):
        plt.figure(f.number)
        plt.suptitle(ttl, size=18)
        plt.savefig(out_fname % tag)
        plt.close()

# Make ztot figure
fig, axes = plt.subplots(nrows=2, ncols=ncols, figsize=(12, 8))
plt.subplots_adjust(left=0.06, bottom=0.08, right=0.99, top=0.9, wspace=0.4)
for i, v in enumerate(fake_varnames):
    ax = axes[int(i/ncols), i%ncols]
        
    pcts = [0, 10, 25, 50, 75, 90, 100]
    var_percentiles = {}
    for p in pcts:
        var_percentiles[p] = np.zeros(ntimes)
        for l in range(ntimes):
            var_percentiles[p][l] = np.nanpercentile(all_zscores[v][:, l], p) 
        
    ax.plot(plot_hr, var_percentiles[50], 'k-', lw=2)
    ax.fill_between(plot_hr, var_percentiles[25], var_percentiles[75], color='k', alpha=0.4)
    ax.fill_between(plot_hr, var_percentiles[10], var_percentiles[90], color='k', alpha=0.2)
    ax.plot(plot_hr, var_percentiles[0], 'k-', lw=0.5)
    ax.plot(plot_hr, var_percentiles[100], 'k-', lw=0.5)

    ax.grid()
    ax.axhline(zcutoff, c='r', ls='-', lw=1.5)
    ax.axhline(-zcutoff, c='r', ls='-', lw=1.5)
    ax.set_ylabel('%s z-score' % v, size=14)
    ax.set_xlim([plot_hr.min(), plot_hr.max()])
for i in range(ncols):
    if i != 3:
        axes[-1, i].set_xlabel('hour', size=14)
    axes[0, 3].set_xlabel('hour', size=14)
plt.suptitle('all stations: %s z-score' % zscore, size=18)
plt.savefig(out_fname % ('all_z_%s' % zscore))
plt.close()


"""
End compare_NR_real_sfc_stations.py
"""
