"""
Perform Comparisons Between Real and Simulated Surface Station Observations

Real surface station obs come from the Iowa Environmental Mesonet database:
https://mesonet.agron.iastate.edu/request/download.phtml

shawn.s.murdzek@noaa.gov
Date Created: 7 April 2023
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
import geopy.distance as gd


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Parameters for real obs
real_obs_dir = '/work2/noaa/wrfruc/murdzek/real_obs/sfc_stations'
station_ids = ['ABR', 'ALB', 'BHM', 'CRP', 'DTW', 'GJT', 'MIA', 'OAK', 'SLE', 'TUS']
station_ids = ['TUS']
years = np.arange(1993, 2023) 
startdate = '04290000'
enddate = '05070000'

# Parameters for fake obs
# analysis_times are dt.timedelta objects relative to 0000
fake_obs_dir = '/work2/noaa/wrfruc/murdzek/nature_run_spring/synthetic_obs_csv/perfect'
analysis_days = [dt.datetime(2022, 4, 29) + dt.timedelta(days=i) for i in range(8)]
analysis_times = [dt.timedelta(hours=i) for i in range(24)]

# Maximum time allowed between analysis_times and either the real or fake ob (sec)
max_time_allowed = 480.

# Maximum distance allowed between real surface station and simulated surface station (km)
max_dist_allowed = 10.

# Output file name (include %s placeholder for station ID)
out_fname = '%s_sfc_station_compare.png'


#---------------------------------------------------------------------------------------------------
# Compare NR to Surface Station Obs
#---------------------------------------------------------------------------------------------------

# Variables to extract
real_varnames = ['lon', 'lat', 'tmpf', 'dwpf', 'drct', 'sknt', 'alti', 'vsby']
fake_var_thermo = ['TOB', 'QOB', 'POB']
fake_var_wind = ['UOB', 'VOB']
fake_varnames = fake_var_thermo + fake_var_wind

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
    for j, yr in enumerate(years):
        tmp_df = pd.read_csv('%s/%s_%d%s_%d%s.txt' % (real_obs_dir, ID, yr, startdate, yr, enddate), 
                             skiprows=5)
        station_times = [dt.datetime.strptime(s, '%Y-%m-%d %H:%M').replace(year=analysis_year) 
                         for s in tmp_df['valid'].values]
        station_times_tot = np.array([(s - dt.datetime(analysis_year, 1, 1)).total_seconds() 
                                      for s in station_times])
        for k, d in enumerate(analysis_days):
            for l, time in enumerate(analysis_times):
                diff = np.abs(station_times_tot - ((d + time) - dt.datetime(analysis_year, 1, 1)).total_seconds())
                idx = np.argmin(diff)
                for v in real_varnames:
                    if (diff.min() > max_time_allowed) or (tmp_df.loc[idx, v] == 'M'):
                        real_stations[ID][v][(j*ndays) + k, l] = np.nan
                    else:
                        real_stations[ID][v][(j*ndays) + k, l] = tmp_df.loc[idx, v]

# Convert real obs to same units/variables as fake obs
for ID in station_ids:
    real_stations[ID]['POB'] = (real_stations[ID]['alti'] * units.inHg).to('mbar').magnitude
    real_stations[ID]['TOB'] = (real_stations[ID]['tmpf'] * units.degF).to('degC').magnitude       
    real_stations[ID]['QOB'] = mc.specific_humidity_from_dewpoint(real_stations[ID]['POB'] * units.mbar,
                                                                  real_stations[ID]['dwpf'] * units.degF).magnitude * 1e6
    wnd_tmp = mc.wind_components(real_stations[ID]['sknt'] * units.kt, real_stations[ID]['drct'] * units.deg)
    real_stations[ID]['UOB'] = wnd_tmp[0].to(units.m / units.s).magnitude
    real_stations[ID]['VOB'] = wnd_tmp[1].to(units.m / units.s).magnitude
    real_stations[ID]['lon'] = real_stations[ID]['lon'] + 360.

# Extract fake observations
print()
fake_stations = {}
for ID in station_ids:
    fake_stations[ID] = {}
    for v in fake_varnames:
        fake_stations[ID][v] = np.ones([ndays, ntimes]) * np.nan
for i, d in enumerate(analysis_days):
    for j, t in enumerate(analysis_times):
        time = d + t
        print(time.strftime('%Y%m%d %H:%M'))
        try:
            full_bufr_csv = bufr.bufrCSV('%s/%s.rap.fake.prepbufr.csv' % 
                                         (fake_obs_dir, time.strftime('%Y%m%d%H%M')))
        except FileNotFoundError:
            print('file not found, continuing to next time')
            continue
        red_csv = full_bufr_csv.df.loc[((full_bufr_csv.df['subset'] == 'ADPSFC') | 
                                        (full_bufr_csv.df['subset'] == 'MSONET')) & 
                                       ((full_bufr_csv.df['DHR'] >= min_hr) &
                                        (full_bufr_csv.df['DHR'] <= max_hr))]
        red_csv.reset_index(inplace=True, drop=True)
        csv_thermo = red_csv.loc[red_csv['TYP'] < 200]
        csv_wind = red_csv.loc[red_csv['TYP'] > 200]
        csv_thermo.reset_index(inplace=True, drop=True)
        csv_wind.reset_index(inplace=True, drop=True)
        for csv, names in zip([csv_thermo, csv_wind], [fake_var_thermo, fake_var_wind]):
            YOB = csv['YOB'].values
            XOB = csv['XOB'].values
            for ID in station_ids:
                dist2 = (real_stations[ID]['lon'][0, 0] - XOB)**2 + (real_stations[ID]['lat'][0, 0] - YOB)**2
                idx = np.argmin(dist2)
                true_distance = gd.distance((real_stations[ID]['lat'][0, 0], real_stations[ID]['lon'][0, 0]), 
                                            (YOB[idx], XOB[idx])).km 
                if true_distance < max_dist_allowed:
                    for v in names:
                        fake_stations[ID][v][i, j] = red_csv.loc[idx, v]

# Plot results
plot_hr = np.array([t.total_seconds() / 3600. for t in analysis_times])
for ID in station_ids:
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(10, 8), sharex=True)
    plt.subplots_adjust(left=0.08, bottom=0.08, right=0.97, top=0.9, wspace=0.35)
    for j, v in enumerate(fake_varnames):
        ax = axes[int(j/3), j%3]
        for k in range(ndays):
            ax.plot(plot_hr, fake_stations[ID][v][k, :], 'k-')
        ax.grid()
        ax.set_ylabel('%s (%s)' % (v, full_bufr_csv.meta[v]['units']), size=14)
    for j in range(3):
        axes[-1, j].set_xlabel('hour', size=14)
        axes[-1, j].set_xlim([plot_hr.min(), plot_hr.max()])
    plt.suptitle(ID, size=18)
    plt.savefig(out_fname % ID)
    plt.close()


"""
End compare_NR_real_sfc_stations.py
"""
