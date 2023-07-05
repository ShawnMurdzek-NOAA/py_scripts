"""
Plot O-B, O-A, and Number of Obs Assimilated for A Single Platform (Or Collection of Platforms) 
and A Single Variable

shawn.s.murdzek@noaa.gov
Date Created: 5 July 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import pandas as pd

import pyDA_utils.gsi_fcts as gsi


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# O-Bs are found in the "ges" files and O-As are found in the "anl" files
# Can 1 or 2 datasets. Key is the name of the dataset
omb_tmpl_real = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/real_data/winter/NCO_dirs/ptmp/prod/rrfs.%s/%s/diag_conv_uv_ges.%s.nc4'
oma_tmpl_real = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/real_data/winter/NCO_dirs/ptmp/prod/rrfs.%s/%s/diag_conv_uv_anl.%s.nc4'
omb_tmpl_osse = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/tune_conv_ob_err/winter_DEBUG/NCO_dirs/ptmp/prod/rrfs.%s/%s/diag_conv_uv_ges.%s.nc4'
oma_tmpl_osse = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/tune_conv_ob_err/winter_DEBUG/NCO_dirs/ptmp/prod/rrfs.%s/%s/diag_conv_uv_anl.%s.nc4'
dates = [dt.datetime(2022, 2, 1, 9) + dt.timedelta(hours=i) for i in range(4)]

omb_fnames = {}
oma_fnames = {}
omb_fnames['real'] = [omb_tmpl_real % (d.strftime('%Y%m%d'), d.strftime('%H'), d.strftime('%Y%m%d%H')) for d in dates]
oma_fnames['real'] = [oma_tmpl_real % (d.strftime('%Y%m%d'), d.strftime('%H'), d.strftime('%Y%m%d%H')) for d in dates]
omb_fnames['OSSE'] = [omb_tmpl_osse % (d.strftime('%Y%m%d'), d.strftime('%H'), d.strftime('%Y%m%d%H')) for d in dates]
oma_fnames['OSSE'] = [oma_tmpl_osse % (d.strftime('%Y%m%d'), d.strftime('%H'), d.strftime('%Y%m%d%H')) for d in dates]

# Variable to plot (only necessary for the uv diag files)
omf_var = 'v'

# Subset of each observation type to plot ('all' - all obs, 'assim' - only obs that are assimilated)
data_subset = 'all'

# Vertical coordinate for vertical profile ('Pressure' or 'Height')
vcoord = 'Height'

# Observation types
ob_types = [227]

# Output directory and string to add to output file names
out_dir = './'
out_str = 'proflr'


#---------------------------------------------------------------------------------------------------
# Compute Statistics
#---------------------------------------------------------------------------------------------------

# Define variable name in GSI diag file
data_names = list(omb_fnames.keys())
if omb_fnames[data_names[0]][0].split('_')[-2] == 'uv':
    vname = '%s_Obs_Minus_Forecast_adjusted' % omf_var
else:
    vname = 'Obs_Minus_Forecast_adjusted'
    omf_var = omb_fnames[data_names[0]][0].split('_')[-2]

# Settings based on desired vertical coordinate
if vcoord == 'Pressure':
    vbins = np.arange(1050, 95, -10)
    bin_ctr = np.log10(vbins[:-1] + (0.5 * (vbins[1:] - vbins[:-1])))
    vunits = 'hPa'
elif vcoord == 'Height':
    vbins = np.arange(0, 10000, 100)
    bin_ctr = vbins[:-1] + (0.5 * (vbins[1:] - vbins[:-1]))
    vunits = 'm'
nbins = len(vbins) - 1

# Extract data
omf_df = {}
for key in data_names:
    print('Dataset = %s' % key)
    omf_df[key] = {}
    omf_df[key]['omb'] = gsi.read_diag(omb_fnames[key])
    omf_df[key]['oma'] = gsi.read_diag(oma_fnames[key])
omf_dates = np.unique(omf_df[data_names[0]]['omb']['date_time'])

# Only retain desired ob types
all_ob_typ = np.unique(omf_df[data_names[0]]['omb']['Observation_Type'].values)
for o in all_ob_typ:
    if o not in ob_types:
        for key in data_names:
            for omf in ['oma', 'omb']:
                omf_df[key][omf] = omf_df[key][omf].loc[omf_df[key][omf]['Observation_Type'] != o].copy()

# Determine ob counts, O-B, and O-A distributions
print('Computing distributions...')
plot_data = {}
pcts = [0, 10, 25, 50, 75, 90, 100]
for key in data_names:
    plot_data[key] = {}
    plot_data[key]['n_obs'] = np.zeros(nbins)
    plot_data[key]['n_assim'] = np.zeros(nbins)
    plot_data[key]['omb'] = {}
    plot_data[key]['oma'] = {}
    for p in pcts:
        plot_data[key]['omb'][p] = np.zeros([nbins]) * np.nan 
        plot_data[key]['oma'][p] = np.zeros([nbins]) * np.nan
    for j in range(nbins):
        omb_subset = omf_df[key]['omb'].loc[np.logical_and(omf_df[key]['omb'][vcoord] >= vbins[j],
                                                           omf_df[key]['omb'][vcoord] < vbins[j+1])]
        oma_subset = omf_df[key]['oma'].loc[np.logical_and(omf_df[key]['oma'][vcoord] >= vbins[j],
                                                           omf_df[key]['oma'][vcoord] < vbins[j+1])]
        plot_data[key]['n_obs'][j] = len(oma_subset)
        plot_data[key]['n_assim'][j] = np.sum(oma_subset['Analysis_Use_Flag'] == 1)
        if data_subset == 'assim':
            omb_subset = omb_subset.loc[omb_subset['Analysis_Use_Flag'] == 1]
            oma_subset = oma_subset.loc[oma_subset['Analysis_Use_Flag'] == 1]
        if len(omb_subset) > 0:
            for p in pcts:
                plot_data[key]['omb'][p][j] = np.percentile(omb_subset[vname].values, p)
                plot_data[key]['oma'][p][j] = np.percentile(oma_subset[vname].values, p)

# Create Plot
print('Making plot...')
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 8), sharey=True, sharex='col')
plt.subplots_adjust(left=0.12, bottom=0.08, right=0.98, top=0.88, hspace=0.25, wspace=0.1)
    
for key, c in zip(data_names, ['b', 'r', 'g', 'c', 'purple']):
 
    # Ob counts in first column
    for j, n in enumerate(['n_obs', 'n_assim']):
        axes[j, 0].plot(plot_data[key][n], bin_ctr, c=c, lw=2, label=key)

    # O-B and O-A in second column
    for j, omf in enumerate(['oma', 'omb']):
        ax = axes[j, 1]
        ax.plot(plot_data[key][omf][50], bin_ctr, c=c, lw=2)
        ax.fill_betweenx(bin_ctr, plot_data[key][omf][25], plot_data[key][omf][75], color=c, alpha=0.4)
        ax.fill_betweenx(bin_ctr, plot_data[key][omf][10], plot_data[key][omf][90], color=c, alpha=0.2)
        ax.plot(plot_data[key][omf][0], bin_ctr, c=c, lw=0.75)
        ax.plot(plot_data[key][omf][100], bin_ctr, c=c, lw=0.75)

axes[0, 0].set_xlabel('# obs', size=14)
axes[1, 0].set_xlabel('# assimilated', size=14)
axes[0, 1].set_xlabel('O$-$B', size=14)
axes[1, 1].set_xlabel('O$-$A', size=14)
axes[0, 0].legend(fontsize=12)
for i in range(2):
    axes[i, 1].axvline(0, c='k', lw=1)
    axes[i, 0].set_xlim(left=0)
    axes[i, 0].set_ylabel('%s (%s)' % (vcoord, vunits), size=14)
    for j in range(2):
        axes[i, j].grid()    

ttl = ('start = %d, end = %d, var = %s\nsubset = %s, types = ' %
       (np.amin(omf_dates), np.amax(omf_dates), omf_var, data_subset))
for t in ob_types:
    ttl = ttl + ('%s, ' % t)
plt.suptitle(ttl, size=18)

plt.savefig('%s/omf_diag_%s_%s_%s_%d_%d.png' %
            (out_dir, out_str, omf_var, data_subset, np.amin(omf_dates), np.amax(omf_dates)))
plt.close()


"""
End diag_omf_vprof.py  
"""
