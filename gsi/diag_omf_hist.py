"""
Plot Distributions of O-A and O-B For a Single Variable

shawn.s.murdzek@noaa.gov
Date Created: 18 April 2023
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
omb_template = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/tune_conv_ob_err/winter1/NCO_dirs/ptmp/prod/rrfs.%s/%s/diag_conv_uv_ges.%s.nc4'
oma_template = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/tune_conv_ob_err/winter1/NCO_dirs/ptmp/prod/rrfs.%s/%s/diag_conv_uv_anl.%s.nc4'
dates = [dt.datetime(2022, 2, 1, 9) + dt.timedelta(hours=i) for i in range(18)]
omb_fnames = [omb_template % (d.strftime('%Y%m%d'), d.strftime('%H'), d.strftime('%Y%m%d%H')) for d in dates]
oma_fnames = [oma_template % (d.strftime('%Y%m%d'), d.strftime('%H'), d.strftime('%Y%m%d%H')) for d in dates]

# Observation types
ob_types = [280, 281, 282, 284, 287]

# Output directory and string to add to output file names
out_dir = './'
out_str = 'adpsfc_sfcshp_fake'


#---------------------------------------------------------------------------------------------------
# Compute Statistics
#---------------------------------------------------------------------------------------------------

omf_var = omb_fnames[0].split('/')[-1].split('_')[2]
print('Var = %s' % omf_var)

# Extract data
omf_df = {}
omf_df['omb'] = gsi.read_diag(omb_fnames)
omf_df['oma'] = gsi.read_diag(oma_fnames)
omf_dates = np.unique(omf_df['omb']['date_time'])

for omf in ['omb', 'oma']:
    partial_df = []
    for typ in ob_types:
        partial_df.append(omf_df[omf].loc[omf_df[omf]['Observation_Type'] == typ].copy())
    omf_df[omf] = pd.concat(partial_df)

maxvals = []
minvals = []
for omf in ['omb', 'oma']:
    if omf_var == 'uv':
        var_strs = ['u_', 'v_']
    else:
        var_strs = ['']
    for v in var_strs:
        maxvals.append(np.percentile(omf_df[omf]['%sObs_Minus_Forecast_adjusted' % v], 99.5))
        minvals.append(np.percentile(omf_df[omf]['%sObs_Minus_Forecast_adjusted' % v], 0.5))
plot_lims = [np.amin(np.array(minvals)), np.amax(np.array(maxvals))]

# Create figure
nrows = 2
if omf_var == 'uv':
    ncols = 4
    fig = plt.figure(figsize=(14, 8))
else:
    ncols = 2
    fig = plt.figure(figsize=(8, 8))

# Plot data
print('plotting data...')
max_hist_val = [0, 0]
axes = np.empty([nrows, ncols], dtype=object)
for i, (omf, xlabel) in enumerate(zip(['omb', 'oma'], ['O$-$B', 'O$-$A'])):
    for j, v in enumerate(var_strs):
        nx = i + 2*j
        for k, subset in enumerate(['all', 'assim']):
            ax = fig.add_subplot(nrows, ncols, (ncols*k)+nx+1)
            vname = '%sObs_Minus_Forecast_adjusted' % v
            if subset == 'all':
                data = omf_df[omf][vname].values
            elif subset == 'assim':
                data = omf_df[omf][vname].loc[omf_df[omf]['Analysis_Use_Flag'] == 1].values
            hist_out = ax.hist(data, bins=30, range=plot_lims)
            max_hist_val[k] = max(hist_out[0].max(), max_hist_val[k])
            ax.set_xlabel(xlabel, size=12)
            ax.set_title('%s%s ($\mu$=%.3f, $\sigma$=%.3f)' % 
                         (v, subset, np.mean(data), np.std(data)), size=14)
            ax.grid() 
            axes[k, nx] = ax

for i in range(nrows):
    for ax in axes[i, :]:
        ax.set_ylim([0, 1.05*max_hist_val[i]])
    axes[i, 0].set_ylabel('counts', size=12)

plt.subplots_adjust(hspace=0.3, wspace=0.3, left=0.07, bottom=0.07, right=0.97, top=0.85)

# Add additional metrics
ttl = 'TYP ='
for typ in ob_types:
    ttl = '%s %d,' % (ttl, typ)
ttl = ('%s\nstart=%d, end=%d,   n_obs=%d, n_assimilated=%d' % 
       (ttl, np.amin(omf_dates), np.amax(omf_dates), len(omf_df['oma']), 
        np.sum(omf_df['oma']['Analysis_Use_Flag'] == 1))) 

plt.suptitle(ttl, size=14) 

plt.savefig('%s/omf_diag_%s_%s_%d_%d.png' % (out_dir, out_str, omf_var, np.amin(omf_dates), np.amax(omf_dates)))
plt.close() 


"""
End diag_omf_distributions.py
"""
