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
import pandas as pd
import numpy as np
import datetime as dt
import metpy.calc as mc
from metpy.units import units

import pyDA_utils.bufr as bufr


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Surface station files (from IEM)
sfc_station_files = glob.glob('/work2/noaa/wrfruc/murdzek/real_obs/sfc_stations/*/*2022*.txt')

# Output directory and file suffix
out_dir = '/work2/noaa/wrfruc/murdzek/real_obs/sfc_stations/bufr_csv'
out_suffix = '.sfc.prepbufr.csv'

# Add fake obs at DHR = 0? This ensures that simulated surface station data is created at the top
# of every hour
fake_dhr0 = True


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
ss_df['PMO'] = ss_df['mslp']
ss_df['TOB'] = (np.float64(ss_df['tmpf'].values) * units.degF).to('degC').magnitude
ss_df['QOB'] = mc.specific_humidity_from_dewpoint(ss_df['POB'].values * units.mbar,
                                                  np.float64(ss_df['dwpf'].values) * units.degF).magnitude * 1e6
ss_df['TDO'] = (np.float64(ss_df['dwpf'].values) * units.degF).to('degC').magnitude
wnd_tmp = mc.wind_components(np.float64(ss_df['sknt'].values) * units.kt, 
                             np.float64(ss_df['drct'].values) * units.deg)
ss_df['UOB'] = wnd_tmp[0].to(units.m / units.s).magnitude
ss_df['VOB'] = wnd_tmp[1].to(units.m / units.s).magnitude
ss_df['HOVI'] = (np.float64(ss_df['vsby'].values) * units.miles).to('m').magnitude
ss_df['XOB'] = ss_df['lon'] + 360.
ss_df['YOB'] = ss_df['lat']
ss_df['ELV'] = ss_df['elevation']
ss_df['SID'] = ss_df['station']

bufr_fixed_cols = {'nmsg':1, 'subset':'ADPSFC', 'ntb':2, 'TYP':300, 'T29':512, 'CAT':0, 'PQM':2, 
                   'QQM':2, 'TQM':2, 'WQM':2, 'PMQ':2, 'tvflg':1, 'vtcd':0}
for c in bufr_fixed_cols.keys():
    ss_df[c] = [bufr_fixed_cols[c]] * len(ss_df)
bufr_nan_cols = ['SAID','ZOB','PWO', 'MXGS', 'PRSS','PMO','ZQM', 
                 'NUL', 'PWQ', 'POE',   'QOE',    'TOE', 'NUL.1',    'WOE',    'NUL.2', 'PWE',
                 'XDR', 'YDR',   'HRDR',     'MSST',  'DBSS',   'SST1','SSTQM',  'SSTOE',  'TFC', 'UFC',
                 'VFC', 'VSSO',  'CLAM',     'HOCB',  'CDTP',   'TOCC','GCDTT',  'CDTP_QM','MXTM','MITM', 
                 'POAF','IALR',  'PRWE',     'PRVSTG','SPRVSTG','HOWV','CEILING','QIFN',   'HBLCS','TSB', 
                 'ACID']
for c in bufr_nan_cols:
    ss_df[c] = [np.nan] * len(ss_df)

# Create fake obs for DHR = 0 (one for each station)
if fake_dhr0:
    tmp_series = []
    for ID in np.unique(ss_df['SID'].values):
        tmp_series.append(ss_df.loc[ss_df['SID'] == ID].iloc[0])
    fake_df = pd.concat(tmp_series, axis=1).T
    fake_df['DHR'] = 0.
    for v in ['TOB', 'QOB', 'UOB', 'VOB', 'POB']:
        fake_df[v] = 300.    
                 
# Split up surface station DataFrame into hourly prepBUFR files
bufr_col_order = ['nmsg','subset','cycletime','ntb',   'SID',    'XOB', 'YOB',    'DHR',    'TYP', 'ELV', 
                  'SAID','T29',   'POB',      'QOB',   'TOB',    'ZOB', 'UOB',    'VOB',    'PWO', 'MXGS',
                  'HOVI','CAT',   'PRSS',     'TDO',   'PMO',    'PQM', 'QQM',    'TQM',    'ZQM', 'WQM', 
                  'NUL', 'PWQ',   'PMQ',      'POE',   'QOE',    'TOE', 'NUL',    'WOE',    'NUL', 'PWE',
                  'XDR', 'YDR',   'HRDR',     'MSST',  'DBSS',   'SST1','SSTQM',  'SSTOE',  'TFC', 'UFC',
                  'VFC', 'VSSO',  'CLAM',     'HOCB',  'CDTP',   'TOCC','GCDTT',  'CDTP_QM','MXTM','MITM', 
                  'POAF','IALR',  'PRWE',     'PRVSTG','SPRVSTG','HOWV','CEILING','QIFN',   'HBLCS','TSB', 
                  'ACID','tvflg', 'vtcd']
ss_df['time'] = pd.to_datetime(ss_df['valid'])
bufr_times = np.unique(ss_df['time'].dt.strftime('%Y%m%d%H').values)
for t in bufr_times:
    t_datetime = dt.datetime.strptime(t, '%Y%m%d%H')
    full_df = ss_df.loc[(ss_df['time'] >= (t_datetime - dt.timedelta(hours=1))) & 
                        (ss_df['time'] <= t_datetime)].copy()
    full_df['cycletime'] = t
    full_df['DHR'] = (full_df['time'] - t_datetime).dt.total_seconds() / 3600.
    if fake_dhr0:
        fake_df['cycletime'] = t
        full_df = pd.concat([full_df, fake_df])
    bufr_df = full_df.drop(ss_cols, axis=1)
    bufr_df = bufr_df[bufr_col_order]
    bufr.df_to_csv(bufr_df, '%s/%s00%s' % (out_dir, t, out_suffix)) 
        

"""
End create_sfc_station_bufr_csv.py 
"""
