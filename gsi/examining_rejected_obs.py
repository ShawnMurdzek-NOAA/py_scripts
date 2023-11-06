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
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import pyDA_utils.gsi_fcts as gsi


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# O-Bs are found in the "ges" files and O-As are found in the "anl" files
# Can have any number of datasets. Key is the name of the dataset
tmpl_real_winter = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/real_red_data/winter_updated/NCO_dirs/ptmp/prod/rrfs.%s/%s_spinup/'
tmpl_real_spring = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/real_red_data/spring/NCO_dirs/ptmp/prod/rrfs.%s/%s/'
tmpl_osse_winter = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data/winter_updated/NCO_dirs/ptmp/prod/rrfs.%s/%s_spinup/'
tmpl_osse_spring = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data/spring/NCO_dirs/ptmp/prod/rrfs.%s/%s/'
dates_winter = [dt.datetime(2022, 2, 1, 9) + dt.timedelta(hours=i) for i in range(159)]
dates_spring = [dt.datetime(2022, 4, 29, 21) + dt.timedelta(hours=i) for i in range(159)]
dates_winter = [dt.datetime(2022, 2, 1, 3)]

path_tmpl = {}
path_tmpl['real'] = [tmpl_real_winter % (d.strftime('%Y%m%d'), d.strftime('%H')) for d in dates_winter]
path_tmpl['OSSE'] = [tmpl_osse_winter % (d.strftime('%Y%m%d'), d.strftime('%H')) for d in dates_winter]
dates = dates_winter

# Variable to plot
omf_var = 't'

# Observation type
ob_typ = 188

# Output directory and string to add to output file names
out_dir = './'
out_str = 'winter'


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
    if (omf_var == 'u' or omf_var == 'v'):
        tmp_fnames = ['%s/diag_conv_uv_ges.%s.nc4' % (path, d.strftime('%Y%m%d%H')) 
                      for path, d in zip(path_tmpl[key], dates)]
    else:
        tmp_fnames = ['%s/diag_conv_%s_ges.%s.nc4' % (path, omf_var, d.strftime('%Y%m%d%H')) 
                      for path, d in zip(path_tmpl[key], dates)]
    tmp_df = gsi.read_diag(tmp_fnames)
    omf_df[key] = tmp_df.loc[tmp_df['Observation_Type'] == ob_typ].copy()
    omf_df_prep0[key] = omf_df[key].loc[np.isclose(omf_df[key]['Prep_Use_Flag'].values, 0)].copy()
omf_dates = np.unique(omf_df[data_names[0]]['date_time'])
    

#---------------------------------------------------------------------------------------------------
# Plot histogram of O-Bs
#---------------------------------------------------------------------------------------------------
        
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
omb_str = 'Obs_Minus_Forecast_adjusted'        
colors = ['b', 'r']

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
plt.suptitle('Type = %d, Var = %s, Prep_Use_Flag = 0' % (ob_typ, omf_var), size=20)
plt.savefig('%s/omb_hist_%d_%s_%s_%d_%d.png' % 
            (out_dir, ob_typ, out_str, omf_var, np.amin(omf_dates), np.amax(omf_dates)))
plt.close() 


#---------------------------------------------------------------------------------------------------
# Plot histogram of elevation differences
#---------------------------------------------------------------------------------------------------

elev_diff = (omf_df_prep0[data_names[0]]['Station_Elevation'].values - 
             omf_df_prep0[data_names[1]]['Station_Elevation'].values)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
colors = ['b', 'r', 'g', 'c']

# Determine bins
min_bin = np.percentile(elev_diff, 0.5)
max_bin = np.percentile(elev_diff, 99.5)
bins = np.linspace(min_bin, max_bin, 51)

i = 0
for flag1, label1 in zip([1, -1], ['assim', 'rej/mon']):
    for flag2, label2 in zip([1, -1], ['assim', 'rej/mon']):
        idx = np.where(np.logical_and(omf_df_prep0[data_names[0]]['Analysis_Use_Flag'] == flag1,
                                      omf_df_prep0[data_names[1]]['Analysis_Use_Flag'] == flag2))[0]
        hvals = np.histogram(elev_diff[idx], bins=bins, density=True)[0]
        ax.plot(0.5*(bins[1:] + bins[:-1]), hvals, c=colors[i], ls='-', 
                label='%s %s, %s %s' % (data_names[0], label1, data_names[1], label2))
        i = i + 1     
   
ax.grid()
ax.legend()
ax.set_xlabel('elevation difference (m) (%s $-$ %s)' % (data_names[0], data_names[1]), size=16)
ax.set_ylabel('density', size=16)
plt.suptitle('Type = %d, Var = %s, Prep_Use_Flag = 0' % (ob_typ, omf_var), size=20)
plt.savefig('%s/elev_hist_%d_%s_%s_%d_%d.png' % 
            (out_dir, ob_typ, out_str, omf_var, np.amin(omf_dates), np.amax(omf_dates)))
plt.close() 


#---------------------------------------------------------------------------------------------------
# Plot histogram of pressure differences
#---------------------------------------------------------------------------------------------------

prs_diff = {}
for key in data_names:
    prs = omf_df[key]['Pressure'].values
    tmp = []
    ct = 0
    for sid in omf_df[key]['Station_ID'].unique():
        s_idx = np.where(omf_df[key]['Station_ID'].values == sid)[0]
        if len(s_idx) > 1:
            ct = ct + 1
            for j, s_j in enumerate(s_idx[:-1]):
                for s_k in s_idx[j+1:]:
                    tmp.append(prs[s_j] - prs[s_k])
    prs_diff[key] = np.array(tmp)
    print(ct)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
colors = ['b', 'r', 'g', 'c']

# Determine bins
bins = np.arange(-1, 1.01, 0.1)

for key, c in zip(data_names, colors):
    hvals = np.histogram(prs_diff[key], bins=bins, density=False)[0]
    ax.plot(0.5*(bins[1:] + bins[:-1]), hvals, c=c, ls='-', label=key)
   
ax.grid()
ax.legend()
ax.set_xlabel('pressure differences between adjacent\nmeasurements at the same station (hPa)', size=16)
ax.set_ylabel('count', size=16)
plt.suptitle('Type = %d, Var = %s' % (ob_typ, omf_var), size=20)
plt.savefig('%s/prs_hist_%d_%s_%s_%d_%d.png' % 
            (out_dir, ob_typ, out_str, omf_var, np.amin(omf_dates), np.amax(omf_dates)))
plt.close() 


#---------------------------------------------------------------------------------------------------
# Plot Spatial Distribution of Rejected Obs
#---------------------------------------------------------------------------------------------------

scale='50m'
fig = plt.figure(figsize=(6, 4+(2*len(data_names))))
for i, key in enumerate(data_names):
    ax = fig.add_subplot(len(data_names), 1, i+1, projection=ccrs.LambertConformal())
    rej_data = omf_df_prep0[key].loc[omf_df_prep0[key]['Analysis_Use_Flag'] == flag].copy()
    ax.plot(rej_data['Longitude'], rej_data['Latitude'], 'b.', transform=ccrs.PlateCarree())
    ax.coastlines(scale)
    borders = cfeature.NaturalEarthFeature(category='cultural',
                                           scale=scale,
                                           facecolor='none',
                                           name='admin_1_states_provinces')
    ax.add_feature(borders)
    ax.set_title('%s rej/mon' % key, size=18)

plt.savefig('%s/omb_map_%d_%s_%s_%d_%d.png' % 
            (out_dir, ob_typ, out_str, omf_var, np.amin(omf_dates), np.amax(omf_dates)))
plt.close() 

# Print Prep_Use_Flags
for key in data_names:
    print()
    print('Prep_Use_Flag counts for %s:' % key)
    print(gsi.gsi_flags_table(omf_df[key], field='Prep_Use_Flag'))


"""
End examining_rejected_obs.py
"""
