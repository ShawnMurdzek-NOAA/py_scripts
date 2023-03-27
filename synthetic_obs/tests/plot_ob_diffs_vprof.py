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

# Input BUFR CSV directory
bufr_dr = '/work2/noaa/wrfruc/murdzek/nature_run_spring/synthetic_obs_csv/perfect'

# Range of datetimes to use for the comparison
date_range = [dt.datetime(2022, 4, 29, 12) + dt.timedelta(hours=i) for i in range(12)]

# Output file name
save_fname = './ob_diffs_aircraft_vprof.png'

# Observation subsets
subsets = ['ADPUPA', 'AIRCFT', 'AIRCAR']
#subsets = ['ADPUPA']

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
real_ob_dfs = []
sim_ob_dfs = []
for d in date_range:
    date_str = d.strftime('%Y%m%d%H%M')
    real_bufr_csv = bufr.bufrCSV('%s/%s.rap.real_red.prepbufr.csv' % (bufr_dir, date_str))
    real_ob_dfs.append(real_bufr_df.df)
    sim_bufr_csv = bufr.bufrCSV('%s/%s.rap.fake.prepbufr.csv' % (bufr_dir, date_str))
    sim_ob_dfs.append(real_bufr_df.df)
    meta = sim_bufr_csv.meta
bufr_df_real = pd.concat(real_ob_dfs)
bufr_df_sim = pd.concat(sim_ob_dfs)

# Only retain obs with DHR between 0 and -1 to prevent double-counting
bufr_df_real = bufr_df_real.loc[np.logical_and(bufr_df_real['DHR'] > -1, bufr_df_real['DHR'] <= 0)]
bufr_df_sim = bufr_df_sim.loc[np.logical_and(bufr_df_sim['DHR'] > -1, bufr_df_sim['DHR'] <= 0)]

# Only retain obs from desired subset
boo = np.zeros(len(bufr_df_sim))
for s in subsets:
    boo[bufr_df_sim['subset'] == s] = 1
ind = np.where(boo)
bufr_df_sim = bufr_df_sim.loc[ind]
bufr_df_real = bufr_df_real.loc[ind]

fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(8, 6), sharey=True)
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, hspace=0.3, wspace=0.3)
for i, v in enumerate(obs_vars):
    print('Plotting %s' % v)
    ax = axes[int(i/4), i%4]

    # Only plot if the quality marker is <= 2
    if v in qm.keys():
        cond = np.logical_and(bufr_df_sim[qm[v]] <= 2, bufr_df_real[qm[v]])
        diff = bufr_df_sim.loc[cond, v] - bufr_df_real.loc[cond, v]
        pres = bufr_df_sim.loc[cond, 'POB'].values
    else:
        diff = bufr_df_sim[v] - bufr_df_real[v]
        pres = bufr_df_sim['POB'].values
    
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
    ax.set_xlabel('%s (%s)' % (v, meta[v]['units']), size=12)

for i in range(2):
    axes[i, 0].set_ylabel('pressure (%s)' % meta['POB']['units'], size=12)

plt.suptitle(title, size=16)
plt.savefig(save_fname)
plt.close()


"""
End plot_ob_diffs_vprof.py  
"""
