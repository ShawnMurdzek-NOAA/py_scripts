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

import pyDA_utils.gsi_fcts as gsi


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# O-Bs are found in the "ges" files and O-As are found in the "anl" files
# Can have any number of datasets. Key is the name of the dataset
tmpl_real = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/real_red_data/winter/NCO_dirs/ptmp/prod/rrfs.%s/%s/'
tmpl_osse = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data/winter_perfect/NCO_dirs/ptmp/prod/rrfs.%s/%s/'
tmpl_osse2 = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data/winter_1st_iter_tuning/NCO_dirs/ptmp/prod/rrfs.%s/%s/'
tmpl_osse3 = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data/winter_1st_iter_whitenoise/NCO_dirs/ptmp/prod/rrfs.%s/%s/'
dates = [dt.datetime(2022, 2, 1, 9) + dt.timedelta(hours=i) for i in range(4*24)]

path_tmpl = {}
path_tmpl['real'] = [tmpl_real % (d.strftime('%Y%m%d'), d.strftime('%H')) for d in dates]
path_tmpl['perfect'] = [tmpl_osse % (d.strftime('%Y%m%d'), d.strftime('%H')) for d in dates]
path_tmpl['corr_errors'] = [tmpl_osse2 % (d.strftime('%Y%m%d'), d.strftime('%H')) for d in dates]
#path_tmpl['uncorr_errors'] = [tmpl_osse3 % (d.strftime('%Y%m%d'), d.strftime('%H')) for d in dates]

# Variables to plot
omf_vars = ['t', 'q', 'u', 'v', 'pw', 'ps']

# Subset of each observation type to plot ('all' - all obs, 'assim' - only obs that are assimilated)
data_subset = 'assim'

# GSD sfcobs uselist file. Used to add the provider string for mesonet obs. Set to None to not 
# use this feature (Note: Using this feature causes the program to run much slower!)
#sfcobs_uselist = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/tune_conv_ob_err/winter_DEBUG/NCO_dirs/stmp/2022020103/anal_conv_gsi_spinup/gsd_sfcobs_uselist.txt'
sfcobs_uselist = None

# Option to group observation classes into larger subsets
use_subsets = False
ob_subsets = {'raob':[120, 122, 132, 220, 221, 222],
              'aircft':[130, 131, 133, 134, 135, 230, 231, 232, 233, 234, 235],
              'sfc':[180, 181, 182, 183, 187, 188, 192, 193, 194, 195, 280, 281, 282, 284, 287, 
                     288, 292, 293, 294, 295],
              'gps':[153]}

# Output directory and string to add to output file names
out_dir = './'
out_str = ''


#---------------------------------------------------------------------------------------------------
# Compute Statistics
#---------------------------------------------------------------------------------------------------

for var in omf_vars:

    print()
    print('-------------------')
    print('Making plots for %s' % var)
    print()

    data_names = list(path_tmpl.keys())
    if (var == 'u' or var == 'v'):
        vname = '%s_Obs_Minus_Forecast_adjusted' % var
    else:
        vname = 'Obs_Minus_Forecast_adjusted'

    # Extract data
    omf_df = {}
    for key in data_names:
        print('Dataset = %s' % key)
        omf_df[key] = {}
        for key2, label in zip(['omb', 'oma'], ['ges', 'anl']):
            if (var == 'u' or var == 'v'):
                tmp_fnames = ['%s/diag_conv_uv_%s.%s.nc4' % (path, label, d.strftime('%Y%m%d%H')) 
                              for path, d in zip(path_tmpl[key], dates)]
            else:
                tmp_fnames = ['%s/diag_conv_%s_%s.%s.nc4' % (path, var, label, d.strftime('%Y%m%d%H')) 
                              for path, d in zip(path_tmpl[key], dates)]
            omf_df[key][key2] = gsi.read_diag(tmp_fnames, mesonet_uselist=sfcobs_uselist)
    omf_dates = np.unique(omf_df[data_names[0]]['omb']['date_time'])
    
    # Create list of ob types
    ob_typ = np.unique(omf_df[data_names[0]]['omb']['Observation_Type'].values)
    if len(data_names) > 1:
        for key in data_names[1:]:
            ob_typ = np.unique(np.concatenate([ob_typ, omf_df[key]['omb']['Observation_Type'].values]))
    ob_groups = ob_typ    

    # Create column for ob subsets
    if use_subsets:
        for key in data_names:
            for omf in ['oma', 'omb']:
                subset = np.array(['a'*50] * len(omf_df[key][omf]))
                for s in ob_subsets.keys():
                    for typ in ob_subsets[s]:
                        subset[omf_df[key][omf]['Observation_Type'] == typ] = s
                omf_df[key][omf]['subset'] = subset
        ob_groups = list(ob_subsets.keys())

    # Determine ob counts and O-B and O-A distributions
    plot_data = {}
    for key in data_names:
        plot_data[key] = {}
        plot_data[key]['n_obs'] = np.zeros(len(ob_groups))
        plot_data[key]['n_assim'] = np.zeros(len(ob_groups))
        plot_data[key]['omb'] = []
        plot_data[key]['oma'] = []
        for j, g in enumerate(ob_groups):
            if use_subsets:
                omb_subset = omf_df[key]['omb'].loc[omf_df[key]['omb']['subset'] == g]
                oma_subset = omf_df[key]['oma'].loc[omf_df[key]['oma']['subset'] == g]
            else:
                omb_subset = omf_df[key]['omb'].loc[omf_df[key]['omb']['Observation_Type'] == g]
                oma_subset = omf_df[key]['oma'].loc[omf_df[key]['oma']['Observation_Type'] == g]
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
    ylocs = np.arange(len(ob_groups))
    nsims = len(data_names)
    bar_hgt = 0.8 / nsims
    colors = ['b', 'r', 'c', 'g', 'goldenrod', 'purple']
    if (nsims % 2) == 0:
        offsets = bar_hgt * np.arange(-(nsims-1), nsims, 2)  / 2.
    else:
        offsets = bar_hgt * np.arange(-int(nsims / 2), int(nsims / 2) + 0.5)
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
    plt.yticks(ticks=ylocs, labels=ob_groups)

    for i in range(4):
        axes[i].grid(axis='x')

    plt.subplots_adjust(wspace=0.1, left=0.07, bottom=0.07, right=0.97, top=0.92)

    # Add additional metrics
    ttl = ('start = %d, end = %d, var = %s, subset = %s' % 
           (np.amin(omf_dates), np.amax(omf_dates), var, data_subset)) 

    plt.suptitle(ttl, size=18) 

    plt.savefig('%s/omf_diag_%s_%s_%s_%d_%d.png' % 
                (out_dir, out_str, var, data_subset, np.amin(omf_dates), np.amax(omf_dates)))
    plt.close() 

    # Print Prep_Use_Flags
    for key in data_names:
        print()
        print('Prep_Use_Flag counts for %s:' % key)
        print(gsi.gsi_flags_table(omf_df[key]['oma'], field='Prep_Use_Flag'))


"""
End diag_omf_distributions.py
"""
