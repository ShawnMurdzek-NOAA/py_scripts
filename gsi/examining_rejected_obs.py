"""
Various Plots and Analyses to Try to Explain Why Certain Obs are Not Assimilated

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import pickle
import copy

import pyDA_utils.gsi_fcts as gsi


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# O-Bs are found in the "ges" files and O-As are found in the "anl" files
# Can have any number of datasets. Key is the name of the dataset
tmpl_real_winter = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/real_red_data/winter_updated/NCO_dirs/ptmp/prod/rrfs.%s/%s/'
tmpl_real_spring = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/real_red_data/spring/NCO_dirs/ptmp/prod/rrfs.%s/%s/'
tmpl_osse_winter = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data/winter_updated/NCO_dirs/ptmp/prod/rrfs.%s/%s/'
tmpl_osse_spring = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data/spring/NCO_dirs/ptmp/prod/rrfs.%s/%s/'
dates_winter = [dt.datetime(2022, 2, 1, 9) + dt.timedelta(hours=i) for i in range(159)]
dates_spring = [dt.datetime(2022, 4, 29, 21) + dt.timedelta(hours=i) for i in range(159)]
dates_spring = [dates_spring[0]]

path_tmpl = {}
path_tmpl['real'] = [tmpl_real_spring % (d.strftime('%Y%m%d'), d.strftime('%H')) for d in dates_spring]
path_tmpl['OSSE'] = [tmpl_osse_spring % (d.strftime('%Y%m%d'), d.strftime('%H')) for d in dates_spring]
dates = dates_spring

# Variable to plot
omf_var = 't'

# Observation type
ob_typ = 188

# Output directory and string to add to output file names
out_dir = './'
out_str = 'spring'


#---------------------------------------------------------------------------------------------------
# Extract Data
#---------------------------------------------------------------------------------------------------

data_names = list(path_tmpl.keys())
if (omf_var == 'u' or omf_var == 'v'):
    vname = '%s_Obs_Minus_Forecast_adjusted' % omf_var
else:
    vname = 'Obs_Minus_Forecast_adjusted'

# Extract data
omf_df = {}
omf_df_prep0 = {}
for key in data_names:
    print('Dataset = %s' % key)
    omf_df[key] = {}
    if (var == 'u' or var == 'v'):
        tmp_fnames = ['%s/diag_conv_uv_ges.%s.nc4' % (path, label, d.strftime('%Y%m%d%H')) 
                      for path, d in zip(path_tmpl[key], dates)]
    else:
        tmp_fnames = ['%s/diag_conv_ges_%s.%s.nc4' % (path, omf_var, label, d.strftime('%Y%m%d%H')) 
                      for path, d in zip(path_tmpl[key], dates)]
    tmp_df = gsi.read_diag(tmp_fnames, mesonet_uselist=sfcobs_uselist)
    omf_df[key] = tmp_df.loc[tmp_df['Observation_Type'] == ob_typ].copy()
    omf_df_prep0[key] = omf_df[key].loc[np.isclose(omf_df[key]['Prep_Use_Flag'].values, 0)].copy()
omf_dates = np.unique(omf_df[data_names[0]]['date_time'])
    

#---------------------------------------------------------------------------------------------------
# Plot histogram of O-Bs for a single ob type
#---------------------------------------------------------------------------------------------------
        
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
omb_str = 'Obs_Minus_Forecast_adjusted'        

# Determine bins
for i, key in enumerate(data_names):
    if i == 0:
        min_bin = np.percentile(omf_df_prep0[key][omb_str], 0.1)     
        max_bin = np.percentile(omf_df_prep0[key][omb_str], 99.9)
    else:
        min_bin = min(min_bin, np.percentile(omf_df_prep0[key][omb_str], 0.1))       
        max_bin = max(max_bin, np.percentile(omf_df_prep0[key][omb_str], 99.9))       
    bins = np.linspace(min_bin, max_bin, 51)

for key, c in zip(data_names, colors):
    for flag, label, ls in zip([1, -1], ['assim', 'rej/mon'], ['-', '--']):
        hvals = np.histogram(omf_df_prep0[key].loc[omf_df_prep0[key]['Analysis_Use_Flag'] == flag][omb_str],
                             bins=bins, density=False)[0]
        ax.plot(0.5*(bins[1:] + bins[:-1]), hvals, c=c, ls=ls, label='%s %s' % (key, label))
        
ax.grid()
ax.legend()
ax.set_xlabel('O$-$B', size=16)
ax.set_ylabel('counts', size=16)
plt.suptitle('Type = %d, Var = %s, Prep_Use_Flag = 0' % (hist_ob_typ, var), size=20)
plt.savefig('%s/omb_hist_%d_%s_%s_%s_%d_%d.png' % 
            (out_dir, hist_ob_typ, out_str, var, data_subset, np.amin(omf_dates), np.amax(omf_dates)))
plt.close() 

# Print Prep_Use_Flags
for key in data_names:
    print()
    print('Prep_Use_Flag counts for %s:' % key)
    print(gsi.gsi_flags_table(omf_df[key]['oma'], field='Prep_Use_Flag'))


"""
End examining_rejected_obs.py
"""
