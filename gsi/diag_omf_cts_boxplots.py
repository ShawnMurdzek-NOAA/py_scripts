"""
Plot Distributions of Ob Counts and O-A and O-B Boxplots for all Types for a Single Variable

shawn.s.murdzek@noaa.gov
Date Created: 18 May 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import gsi_fcts as gsi


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# O-Bs are found in the "ges" files and O-As are found in the "anl" files
# Can 1 or 2 datasets. Key is the name of the dataset
omb_tmpl_real = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/real_data/winter/NCO_dirs/ptmp/prod/rrfs.%s/%s/diag_conv_pw_ges.%s.nc4'
oma_tmpl_real = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/real_data/winter/NCO_dirs/ptmp/prod/rrfs.%s/%s/diag_conv_pw_anl.%s.nc4'
omb_tmpl_osse = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/tune_conv_ob_err/winter1/NCO_dirs/ptmp/prod/rrfs.%s/%s/diag_conv_pw_ges.%s.nc4'
oma_tmpl_osse = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/tune_conv_ob_err/winter1/NCO_dirs/ptmp/prod/rrfs.%s/%s/diag_conv_pw_anl.%s.nc4'
dates = [dt.datetime(2022, 2, 1, 9) + dt.timedelta(hours=i) for i in range(18)]

omb_fnames = {}
oma_fnames = {}
omb_fnames['real'] = [omb_tmpl_real % (d.strftime('%Y%m%d'), d.strftime('%H'), d.strftime('%Y%m%d%H')) for d in dates]
oma_fnames['real'] = [oma_tmpl_real % (d.strftime('%Y%m%d'), d.strftime('%H'), d.strftime('%Y%m%d%H')) for d in dates]
omb_fnames['OSSE'] = [omb_tmpl_osse % (d.strftime('%Y%m%d'), d.strftime('%H'), d.strftime('%Y%m%d%H')) for d in dates]
oma_fnames['OSSE'] = [oma_tmpl_osse % (d.strftime('%Y%m%d'), d.strftime('%H'), d.strftime('%Y%m%d%H')) for d in dates]

# Variable to plot (only necessary for the uv diag files)
omf_var = 'v'

# Subset of each observation type to plot ('all' - all obs, 'assim' - only obs that are assimilated)
data_subset = 'assim'

# Output directory and string to add to output file names
out_dir = './'
out_str = ''


#---------------------------------------------------------------------------------------------------
# Compute Statistics
#---------------------------------------------------------------------------------------------------

data_names = list(omb_fnames.keys())
if omb_fnames[data_names[0]][0].split('_')[-2] == 'uv':
    vname = '%s_Obs_Minus_Forecast_adjusted' % omf_var
else:
    vname = 'Obs_Minus_Forecast_adjusted'
    omf_var = omb_fnames[data_names[0]][0].split('_')[-2]

# Extract data
omf_df = {}
for key in data_names:
    print('Dataset = %s' % key)
    omf_df[key] = {}
    omf_df[key]['omb'] = gsi.read_diag(omb_fnames[key])
    omf_df[key]['oma'] = gsi.read_diag(oma_fnames[key])
omf_dates = np.unique(omf_df[data_names[0]]['omb']['date_time'])

# Create list of ob types
ob_typ = np.unique(omf_df[data_names[0]]['omb']['Observation_Type'].values)
if len(data_names) > 1:
    for key in data_names[1:]:
        ob_typ = np.unique(np.concatenate([ob_typ, omf_df[key]['omb']['Observation_Type'].values]))

# Determine ob counts and O-B and O-A distributions
plot_data = {}
for key in data_names:
    plot_data[key] = {}
    plot_data[key]['n_obs'] = np.zeros(len(ob_typ))
    plot_data[key]['n_assim'] = np.zeros(len(ob_typ))
    plot_data[key]['omb'] = []
    plot_data[key]['oma'] = []
    for j, t in enumerate(ob_typ):
        omb_subset = omf_df[key]['omb'].loc[omf_df[key]['omb']['Observation_Type'] == t]
        oma_subset = omf_df[key]['oma'].loc[omf_df[key]['oma']['Observation_Type'] == t]
        plot_data[key]['n_obs'][j] = len(oma_subset)
        plot_data[key]['n_assim'][j] = np.sum(oma_subset['Analysis_Use_Flag'] == 1)
        if data_subset == 'all':
            plot_data[key]['omb'].append(omb_subset[vname].values)
            plot_data[key]['oma'].append(oma_subset[vname].values)
        elif data_subset == 'assim':
            plot_data[key]['omb'].append(omb_subset[vname].loc[omb_subset['Analysis_Use_Flag'] == 1].values)
            plot_data[key]['oma'].append(oma_subset[vname].loc[oma_subset['Analysis_Use_Flag'] == 1].values)

# Plot data
fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(12, 9), sharey=True)
ylocs = np.arange(len(ob_typ))
bar_hgt = 0.38
if len(data_names) == 1:
    offsets = [0]
    colors = ['b']
elif len(data_names) == 2:
    offsets = [bar_hgt / 2, -bar_hgt / 2]
    colors = ['b', 'r']
for key, off, c in zip(data_names, offsets, colors):
    for j, (n_name, xlabel) in enumerate(zip(['n_obs', 'n_assim'], ['# obs', '# assimilated'])):
        ax = axes[j]
        ax.barh(ylocs+off, plot_data[key][n_name], height=bar_hgt, label=key, color=c)
        ax.set_xlabel(xlabel, size=14)
        ax.set_xscale('log')
        ax.set_xlim(left=1)
    for j, (omf, xlabel) in enumerate(zip(['omb', 'oma'], ['O$-$B', 'O$-$A'])):
        ax = axes[j+2]
        ax.boxplot(plot_data[key][omf], positions=(ylocs+off), widths=bar_hgt, vert=False, patch_artist=True,
                   boxprops={'facecolor':c}, showfliers=False)
        ax.set_xlabel(xlabel, size=14)

axes[0].set_ylabel('observation type', size=14)
axes[0].legend(fontsize=12)
plt.sca(axes[0])
plt.yticks(ticks=ylocs, labels=ob_typ)

for i in range(4):
    axes[i].grid(axis='x')

plt.subplots_adjust(wspace=0.1, left=0.07, bottom=0.07, right=0.97, top=0.92)

# Add additional metrics
ttl = ('start = %d, end = %d, var = %s, subset = %s' % 
       (np.amin(omf_dates), np.amax(omf_dates), omf_var, data_subset)) 

plt.suptitle(ttl, size=18) 

plt.savefig('%s/omf_diag_%s_%s_%s_%d_%d.png' % 
            (out_dir, out_str, omf_var, data_subset, np.amin(omf_dates), np.amax(omf_dates)))
plt.close() 


"""
End diag_omf_distributions.py
"""
