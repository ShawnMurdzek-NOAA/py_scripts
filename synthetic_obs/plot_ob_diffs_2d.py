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
#fname1 = '/scratch1/BMC/wrfruc/murdzek/nature_run_tests/nature_run_spring_v2/synthetic_obs/202204291200.fake.prepbufr.csv'
#fname2 = '/scratch1/BMC/wrfruc/murdzek/nature_run_tests/nature_run_spring_v2/synthetic_obs/202204291200.real_red.prepbufr.csv'

# Dataset names
name1 = 'Sim Obs'
name2 = 'Real Obs'

# Output file name (include %s placeholder for variable)
save_fname = './ob_diffs_%s.png'

# Observation subsets
subsets = ['SHPSFC', 'ADPSFC', 'MSONET']

# Variables to plot
obs_vars = ['ELV', 'POB', 'TOB', 'QOB', 'UOB', 'VOB', 'ZOB', 'XOB', 'YOB']


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

# Load borders
borders = cfeature.NaturalEarthFeature(category='cultural',
                                       scale='50m',
                                       facecolor='none',
                                       name='admin_1_states_provinces')

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
    fig = plt.figure(figsize=(12, 10))

    # Only plot if the quality marker is <= 2
    if v in qm.keys():
        cond = np.logical_and(bufr_df1.df[qm[v]] <= 2, bufr_df2.df[qm[v]])
        field1 = bufr_df1.loc[cond, v]
        field2 = bufr_df2.loc[cond, v]
        diff = field1 - field2
        lat = bufr_df1.df.loc[cond, 'YOB'].values
        lon = bufr_df1.df.loc[cond, 'XOB'].values
    else:
        field1 = bufr_df1.df[v]
        field2 = bufr_df2.df[v]
        diff = field1 - field2
        lat = bufr_df1.df['YOB'].values
        lon = bufr_df1.df['XOB'].values

    # Plot actual values
    maxval = max(np.amax(field1), np.amax(field2))
    minval = max(np.amin(field1), np.amin(field2))
    for j, (f, ttl) in enumerate(zip([field1, field2], [name1, name2])):
        ax1 = fig.add_subplot(2, 3, j+1, projection=ccrs.PlateCarree())
        cax = ax1.scatter(lon, lat, s=25, c=f, cmap='plasma', vmin=minval, vmax=maxval,
                          transform=ccrs.PlateCarree())
        ax1.coastlines()
        ax1.add_feature(borders, linewidth=0.5, edgecolor='k')
        ax1.set_xlim([-135, -50])
        ax1.set_ylim([20, 55])
        ax1.set_title(ttl, size=18)

        ax2 = fig.add_subplot(2, 3, j+3)
        ax2.hist(f, bins=20, range=(minval, maxval))
        ax2.set_xlabel('%s (%s)' % (v, bufr_df1.meta[v]['units']), size=12)
        ax2.set_ylabel('counts', size=12)
        ax2.grid()

    # Plot differences
    maxdiff = np.amax(np.abs(diff))
    ax1 = fig.add_subplot(2, 3, 3, projection=ccrs.PlateCarree())
    cax = ax1.scatter(lon, lat, s=75, c=diff, cmap='bwr', vmin=-maxdiff, vmax=maxdiff,
                      transform=ccrs.PlateCarree(), edgecolors='k', linewidths=0.5)
    cbar = plt.colorbar(cax, ax=ax1, orientation='horizontal')
    cbar.set_label('%s diffs (%s)' % (v, bufr_df1.meta[v]['units']), size=14)

    ax1.coastlines('50m')
    ax1.add_feature(borders, linewidth=0.5, edgecolor='k')
    ax1.set_xlim([-135, -50])
    ax1.set_ylim([20, 55])
    ax1.set_title('%s $-$ %s' % (name1, name2), size=18)
        
    ax2 = fig.add_subplot(2, 3, 6)
    ax2.hist(diff, bins=20, range=(-maxdiff, maxdiff))
    ax2.set_xlabel('%s diff (%s)' % (v, bufr_df1.meta[v]['units']), size=12)
    ax2.set_ylabel('counts', size=12)
    ax2.grid()

    plt.savefig(save_fname % v)
    plt.close()


"""
End plot_ob_diffs_2d.py  
"""
