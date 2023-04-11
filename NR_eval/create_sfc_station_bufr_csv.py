"""
Create PrepBUFR CSVs for Iowa Environmental Mesonet Surface Station Data

Unlike actual prepBUFR files, wind and thermodynamic obs are not separated. PrepBUFR files are
created for every hour and include obs from the previous hour (so -1 <= DHR <= 0).

shawn.s.murdzek@noaa.gov
Date Created: 11 April 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import glob
import bufr
import pandas as pd
import numpy as np
import datetime as dt

import metpy.calc as mc
from metpy.units import units


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Surface station files (from IEM)
sfc_station_files = glob.glob('/work2/noaa/wrfruc/murdzek/real_obs/sfc_stations/*2022*.txt')

# Output directory and file suffix
out_dir = '/work2/noaa/wrfruc/murdzek/real_obs/sfc_stations/bufr_csv'
out_suffix = '.sfc.prepbufr.csv'


#---------------------------------------------------------------------------------------------------
# Create BUFR Files
#---------------------------------------------------------------------------------------------------

# Concatenate all the DataFrames together
ss_df_list = []
for ss_file in sfc_station_files:
    ss_df_list.append(pd.read_csv(ss_file, skiprows=5))
ss_df = pd.concat(ss_df_list)
ss_cols = ss_df.columns

# Switch missing values ("M") to NaN
for c in ss_cols:
    ss_df.loc[ss_df[c] == 'M', c] = np.nan

# Convert surface station DataFrame to correct format for prepBUFR CSVs
ss_df['POB'] = mc.altimeter_to_station_pressure(np.float64(ss_df['alti'].values) * units.inHg, 
                                                ss_df['elevation'].values * units.m).to('mbar').magnitude
ss_df['TOB'] = (np.float64(ss_df['tmpf'].values) * units.degF).to('degC').magnitude
ss_df['QOB'] = mc.specific_humidity_from_dewpoint(ss_df['POB'].values * units.mbar,
                                                  np.float64(ss_df['dwpf'].values) * units.degF).magnitude * 1e6
wnd_tmp = mc.wind_components(np.float64(ss_df['sknt'].values) * units.kt, 
                             np.float64(ss_df['drct'].values) * units.deg)
ss_df['UOB'] = wnd_tmp[0].to(units.m / units.s).magnitude
ss_df['VOB'] = wnd_tmp[1].to(units.m / units.s).magnitude
ss_df['XOB'] = ss_df['lon'] + 360.
ss_df['YOB'] = ss_df['lat']
ss_df['ELV'] = ss_df['elevation']

bufr_fixed_cols = {'nmsg':1, 'subset':'ADPSFC', 'ntb':2, 'SID':ss_file.split('/')[-1][:3], 
                   'TYP':300, 'T29':512, 'CAT':0, 'PQM':2, 'QQM':2, 'TQM':2, 'WQM':2}
for c in bufr_fixed_cols.keys():
    ss_df[c] = [bufr_fixed_cols[c]] * len(ss_df)
bufr_nan_cols = ['SAID', 'ZOB', 'PWO', 'PRSS', 'ZQM', 'NUL', 'PWQ', 'POE', 'QOE', 'TOE', 'NUL.1', 
                 'WOE', 'NUL.2', 'PWE']
for c in bufr_nan_cols:
    ss_df[c] = [np.nan] * len(ss_df)
                 
# Split up surface station DataFrame into hourly prepBUFR files
bufr_col_order = ['nmsg', 'subset', 'cycletime', 'ntb', 'SID', 'XOB', 'YOB', 'DHR', 'TYP', 'ELV', 
                  'SAID', 'T29',    'POB',       'QOB', 'TOB', 'ZOB', 'UOB', 'VOB', 'PWO', 'CAT',
                  'PRSS', 'PQM',    'QQM',       'TQM', 'ZQM', 'WQM', 'NUL', 'PWQ', 'POE', 'QOE',
                  'TOE',  'NUL',    'WOE',       'NUL', 'PWE']
ss_df['time'] = pd.to_datetime(ss_df['valid'])
bufr_times = np.unique(ss_df['time'].dt.strftime('%Y%m%d%H').values)
for t in bufr_times:
    t_datetime = dt.datetime.strptime(t, '%Y%m%d%H')
    full_df = ss_df.loc[(ss_df['time'] >= (t_datetime - dt.timedelta(hours=1))) & 
                        (ss_df['time'] <= t_datetime)].copy()
    full_df['cycletime'] = [t] * len(full_df)
    full_df['DHR'] = (t_datetime - full_df['time']).dt.total_seconds() / 3600.
    bufr_df = full_df.drop(ss_cols, axis=1)
    bufr_df = bufr_df[bufr_col_order]
    bufr.df_to_csv(bufr_df, '%s/%s00%s' % (out_dir, t, out_suffix)) 
        

"""
End create_sfc_station_bufr_csv.py 
"""
