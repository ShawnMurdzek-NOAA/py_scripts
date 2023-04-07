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

# Parameter for fake obs
fake_obs_dir = '/work2/noaa/wrfruc/murdzek/nature_run_spring/synthetic_obs_csv/perfect'
analysis_times = [dt.datetime(2022, 4, 29, 12) + dt.timedelta(hours=i) for i in range(168)]

# Maximum time allowed between matched real and fake obs (sec)
max_time_allowed = 600.

# Output file name (include %s placeholder for station ID)
out_fname = '%s_sfc_station_compare.png'


#---------------------------------------------------------------------------------------------------
# Compare NR to Surface Station Obs
#---------------------------------------------------------------------------------------------------

# Variables to extract
real_varnames = ['lon', 'lat', 'tmpf', 'dwpf', 'drct', 'sknt', 'alti', 'vsby']
fake_varnames = ['TOB', 'QOB', 'UOB', 'VOB', 'POB']

# Convert analysis times to seconds since 1 January
analysis_year = analysis_times[0].year
analysis_days = np.arange(analysis_times[0].day, analysis_times[-1].day+1)
analysis_times_tot_sec = np.array([(t - dt.datetime(analysis_year, 1, 1)).total_seconds() 
                                   for t in analysis_times])

# Extract real surface station obs first
real_stations = {}
for ID in station_ids:
    print('Extracting real obs at %s' % ID)
    real_stations[ID] = {}
    for v in real_varnames:
        real_stations[ID][v] = np.ones([len(years), len(analysis_times)]) * np.nan
    for j, yr in enumerate(years):
        tmp_df = pd.read_csv('%s/%s_%d%s_%d%s.txt' % (real_obs_dir, ID, yr, startdate, yr, enddate), 
                             skiprows=5)
        station_times = [dt.datetime.strptime(s, '%Y-%m-%d %H:%M').replace(year=analysis_year) 
                         for s in tmp_df['valid'].values]
        station_times_tot_sec = np.array([(t - dt.datetime(analysis_year, 1, 1)).total_seconds()
                                          for t in station_times]) 
        for k, time in enumerate(analysis_times_tot_sec):
            diff = np.abs(station_times_tot_sec - time)
            idx = np.argmin(diff)
            for v in real_varnames:
                if (diff.min() > max_time_allowed) or (tmp_df.loc[idx, v] == 'M'):
                    real_stations[ID][v][j, k] = np.nan
                else:
                    real_stations[ID][v][j, k] = tmp_df.loc[idx, v]

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
        fake_stations[ID][v] = np.ones(len(analysis_times)) * np.nan
for i, t in enumerate(analysis_times):
    print(t.strftime('%Y%m%d %H:%M'))
    full_bufr_csv = bufr.bufrCSV('%s/%s.rap.fake.prepbufr.csv' % (fake_obs_dir, t.strftime('%Y%m%d%H%M')))
    red_csv = full_bufr_csv.df.loc[np.logical_or(full_bufr_csv.df['subset'] == 'ADPSFC', 
                                                 full_bufr_csv.df['subset'] == 'MSONET')]
    red_csv.reset_index(inplace=True, drop=True)
    YOB = red_csv['YOB'].values
    XOB = red_csv['XOB'].values
    for ID in station_ids:
        dist2 = (real_stations[ID]['lon'][0, 0] - XOB)**2 + (real_stations[ID]['lat'][0, 0] - YOB)**2
        idx = np.argmin(dist2)
        true_distance = gd.distance((real_stations[ID]['lat'][0, 0], real_stations[ID]['lon'][0, 0]), 
                                    (YOB[idx], XOB[idx])).km 
        if true_distance < 10:
            for v in fake_varnames:
                fake_stations[ID][v][i] = red_csv.loc[idx, v]

# Create diurnal cycles
"""
Current problem:

analysis_days is an empty array b/c the starting day (29) is > the ending day (6). A Julian day 
would solve this, but then we would need to contend with leap years... Maybe I can create my own
"pseudo" Julian day using the month and day? This is starting to become a pain... Maybe I can do
this during data extraction to make things simpler?
"""
diurnal_cycle = {}
for typ in ['fake', 'real']:
    diurnal_cycle[typ] = {}
    for ID in station_ids:
        diurnal_cycle[typ][ID] = {}
        for v in fake_varnames:
            if typ == 'fake':
               yrs = years
            elif typ == 'real':
               yrs = np.array([analysis_year]) 
            diurnal_cycle[typ][ID][v] = np.ones([24, len(yrs) * len(analysis_days)]) * np.nan
            for iyr, yr in enumerate(yrs):
                for itime, t in enumerate(analysis_times):
                    hr = t.hour
                    iday = np.where(t.day == analysis_days)[0]
                    icycle = iyr * len(analysis_days) + iday
                    if typ == 'fake':
                        diurnal_cycle[typ][ID][v][hr, icycle] = fake_stations[ID][v][itime]
                    elif typ == 'real':
                        diurnal_cycle[typ][ID][v][hr, icycle] = real_stations[ID][v][iyr, itime]
         
# Plot results
for ID in station_ids:
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(10, 8), sharey=True)
    for j, v in enumerate(fake_varnames):
        ax = axes[int(j/3), j%3]
        for k in range(diurnal_cycle['fake'][ID][v].shape[1]):
            ax.plot(np.arange(24), diurnal_cycle['fake'][ID][v][:, k], 'k-')
        ax.grid()
        ax.set_ylim('%s (%s)' % (v, full_bufr_csv.meta[v]['units']), size=14)
    for j in range(3):
        axes[-1, j].set_xlim('hour', size=14)
    plt.savefig(out_fname % ID)
    plt.close()


"""
End compare_NR_real_sfc_stations.py
"""
