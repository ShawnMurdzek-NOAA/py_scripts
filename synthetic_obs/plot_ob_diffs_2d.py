"""
Plot Observation Differences Between Two Prepbufr CSV

This script should ideally be used to compare synthetic obs vs. real obs

shawn.s.murdzek@noaa.gov
Date Created: 7 February 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import bufr
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Input BUFR CSV files
fname1 = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_spring_v2/synthetic_obs/202204291200.fake.prepbufr.csv'
fname2 = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_spring_v2/synthetic_obs/202204291200.real_red.prepbufr.csv'

# Output file name (include %s placeholder for variable)
save_fname = './ob_diffs_%s.png'

# Observation subsets
subsets = ['SHPSFC', 'ADPSFC', 'MSONET']

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

for v in obs_vars:
    print('Plotting %s' % v)
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    # Only plot if the quality marker is <= 2
    if v in qm.keys():
        cond = np.logical_and(bufr_df1.df[qm[v]] <= 2, bufr_df2.df[qm[v]])
        diff = bufr_df1.df.loc[cond, v] - bufr_df2.df.loc[cond, v]
        lat = bufr_df1.df.loc[cond, 'YOB'].values
        lon = bufr_df1.df.loc[cond, 'XOB'].values
    else:
        diff = bufr_df1.df[v] - bufr_df2.df[v]
        lat = bufr_df1.df['YOB'].values
        lon = bufr_df1.df['XOB'].values

    maxdiff = np.amax(np.abs(diff))

    cax = ax.scatter(lon, lat, s=75, c=diff, cmap='bwr', vmin=-maxdiff, vmax=maxdiff,
                     transform=ccrs.PlateCarree(), edgecolors='k', linewidths=0.5)
    cbar = plt.colorbar(cax, ax=ax, orientation='horizontal')
    cbar.set_label('%s diffs (%s)' % (v, bufr_df1.meta[v]['units']), size=14)

    ax.coastlines('50m')
    borders = cfeature.NaturalEarthFeature(category='cultural',
                                           scale='50m',
                                           facecolor='none',
                                           name='admin_1_states_provinces')
    ax.add_feature(borders, linewidth=0.5, edgecolor='k')
    ax.set_xlim([-135, -50])
    ax.set_ylim([20, 55])
    ax.set_title(title, size=18)

    plt.savefig(save_fname % v)
    plt.close()


"""
End plot_ob_diffs_2d.py  
"""
