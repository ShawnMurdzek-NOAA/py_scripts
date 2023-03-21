"""
Plot Observation Differences Between Two Prepbufr CSV

This script should ideally be used to compare synthetic obs vs. real obs

shawn.s.murdzek@noaa.gov
Date Created: 14 February 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import bufr
import matplotlib.pyplot as plt
import numpy as np


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Input BUFR CSV files
#fname1 = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_spring/synthetic_obs/202204291200.fake.prepbufr.csv'
#fname2 = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_spring/synthetic_obs/202204291200.real_red.prepbufr.csv'
#fname1 = '/scratch1/BMC/wrfruc/murdzek/nature_run_tests/nature_run_spring_v2/synthetic_obs/202204291200.fake.prepbufr.csv'
#fname2 = '/scratch1/BMC/wrfruc/murdzek/nature_run_tests/nature_run_spring_v2/synthetic_obs/202204291200.real_red.prepbufr.csv'
fname1 = '/work2/noaa/wrfruc/murdzek/nature_run_spring/synthetic_obs/202204291200.fake.prepbufr.csv'
fname2 = '/work2/noaa/wrfruc/murdzek/nature_run_spring/synthetic_obs/202204291200.real_red.prepbufr.csv'

# Output file name
save_fname = './ob_diffs_vprof.png'

# Observation subsets
subsets = ['ADPUPA', 'AIRCFT', 'AIRCAR']

# Variables to plot
obs_vars = ['ELV', 'POB', 'TOB', 'QOB', 'UOB', 'VOB', 'ZOB']

# Bins for binning each variable in the vertical. First entry is left-most edge whereas all 
# subsequent entries are the right edge of the bin. Must be in descending order.
pbins = np.arange(1050, 95, -10)

# Title
title = 'synthetic obs $-$ real obs'


#---------------------------------------------------------------------------------------------------
# Plot BUFR Observation Differences
#---------------------------------------------------------------------------------------------------

# Dictionary giving quality marker fields for each variable (only plot quality markers 0-2)
qm = {'POB':'PQM',
      'QOB':'QQM',
      'TOB':'TQM',
      'ZOB':'ZQM',
      'UOB':'WQM',
      'VOB':'WQM',
      'PWQ':'PWO'}

# Open files
bufr_df1 = bufr.bufrCSV(fname1)
bufr_df2 = bufr.bufrCSV(fname2)

# Apply rounding so precision in simulated obs matches real obs
#bufr_df1.df = bufr.match_bufr_prec(bufr_df1.df)
#bufr_df2.df = bufr.match_bufr_prec(bufr_df2.df)

# Only retain obs from desired subset
boo = np.zeros(len(bufr_df1.df))
for s in subsets:
    boo[bufr_df1.df['subset'] == s] = 1
ind = np.where(boo)
bufr_df1.df = bufr_df1.df.loc[ind]
bufr_df2.df = bufr_df2.df.loc[ind]

fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(8, 6), sharey=True)
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, hspace=0.3, wspace=0.3)
for i, v in enumerate(obs_vars):
    print('Plotting %s' % v)
    ax = axes[int(i/4), i%4]

    # Only plot if the quality marker is <= 2
    if v in qm.keys():
        cond = np.logical_and(bufr_df1.df[qm[v]] <= 2, bufr_df2.df[qm[v]])
        diff = bufr_df1.df.loc[cond, v] - bufr_df2.df.loc[cond, v]
        pres = bufr_df1.df.loc[cond, 'POB'].values
    else:
        diff = bufr_df1.df[v] - bufr_df2.df[v]
        pres = bufr_df1.df['POB'].values
    
    # Compute various percentiles for each bin along the pressure coordinate
    var_percentiles = {}
    pcts = [0, 10, 25, 50, 75, 90, 100]
    for p in pcts:
        var_percentiles[p] = np.ones(len(pbins)-1) * np.nan
    for j, (prs1, prs2) in enumerate(zip(pbins[:-1], pbins[1:])):
        subset = diff[np.logical_and(pres <= prs1, pres > prs2)]
        if len(subset) > 0:
            for p in pcts:
                var_percentiles[p][j] = np.percentile(subset, p)

    bin_ctr = np.log10(pbins[:-1] + (0.5 * (pbins[1:] - pbins[:-1])))
    ax.plot(var_percentiles[50], bin_ctr, 'b-', lw=2)
    ax.fill_betweenx(bin_ctr, var_percentiles[25], var_percentiles[75], color='b', alpha=0.4)
    ax.fill_betweenx(bin_ctr, var_percentiles[10], var_percentiles[90], color='b', alpha=0.2)
    ax.plot(var_percentiles[0], bin_ctr, 'b-', lw=0.75)
    ax.plot(var_percentiles[100], bin_ctr, 'b-', lw=0.75)
    ax.axvline(0, c='k', lw='1')
    ax.grid()

    ax.set_ylim([np.log10(pbins.max()), np.log10(pbins.min())])
    ticks = np.array([1000, 850, 700, 500, 400, 300, 200, 100])
    ax.set_yticks(np.log10(ticks))
    ax.set_yticklabels(ticks)
    ax.set_xlabel('%s (%s)' % (v, bufr_df1.meta[v]['units']), size=12)

for i in range(2):
    axes[i, 0].set_ylabel('pressure (%s)' % bufr_df1.meta['POB']['units'], size=12)

plt.suptitle(title, size=16)
plt.savefig(save_fname)
plt.close()


"""
End plot_ob_diffs_vprof.py  
"""
