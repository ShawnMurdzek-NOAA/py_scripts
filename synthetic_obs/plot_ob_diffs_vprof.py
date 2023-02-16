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
fname1 = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_spring_v2/synthetic_obs/202204291200.fake.prepbufr.csv'
fname2 = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_spring_v2/synthetic_obs/202204291200.real_red.prepbufr.csv'

# Output file name
save_fname = './ob_diffs_vprof.png'

# Observation subsets
subsets = ['ADPUPA', 'AIRCAR', 'AIRCFT']

# Variables to plot
obs_vars = ['ELV', 'POB', 'TOB', 'QOB', 'UOB', 'VOB', 'ZOB']

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
    
    ax.plot(diff, np.log10(pres), 'b.')
    ax.axvline(0, c='k', lw='1')
    ax.grid()

    ax.set_ylim([np.log10(1050), np.log10(100)])
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
