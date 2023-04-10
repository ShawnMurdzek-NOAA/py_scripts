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
import datetime as dt
import pandas as pd


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Input BUFR CSV directory
bufr_dir = '/work2/noaa/wrfruc/murdzek/nature_run_spring/synthetic_obs_csv/perfect'

# Prepbufr file tag (e.g., 'rap', 'rap_e', 'rap_p')
bufr_tag = 'rap_p'

# Range of datetimes to use for the comparison
date_range = [dt.datetime(2022, 4, 29, 12) + dt.timedelta(hours=i) for i in range(13)]
#date_range = [dt.datetime(2022, 4, 29, 12)]

# Dataset names
name1 = 'Sim Obs'
name2 = 'Real Obs'

# Output file name (include %s placeholders for bufr_tag, variable name, and start and end of date 
# range)
save_fname = './ob_diffs_%s_%s_%s_%s.png'

# Observation subsets
subsets = ['SFCSHP', 'ADPSFC', 'MSONET', 'GPSIPW']

# Variables to plot
obs_vars = ['PWO', 'ELV', 'POB', 'TOB', 'QOB', 'UOB', 'VOB', 'ZOB', 'PWO']


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
      'PWO':'PWQ'}

# Load borders
borders = cfeature.NaturalEarthFeature(category='cultural',
                                       scale='50m',
                                       facecolor='none',
                                       name='admin_1_states_provinces')

# Open files
real_ob_dfs = []
sim_ob_dfs = []
for d in date_range:
    date_str = d.strftime('%Y%m%d%H%M')
    try:
        real_bufr_csv = bufr.bufrCSV('%s/%s.%s.real_red.prepbufr.csv' % (bufr_dir, date_str, bufr_tag))
    except FileNotFoundError:
        # Skip to next file
        continue
    real_ob_dfs.append(real_bufr_csv.df)
    sim_bufr_csv = bufr.bufrCSV('%s/%s.%s.fake.prepbufr.csv' % (bufr_dir, date_str, bufr_tag))
    sim_ob_dfs.append(sim_bufr_csv.df)
    meta = sim_bufr_csv.meta
bufr_df_real = pd.concat(real_ob_dfs, ignore_index=True)
bufr_df_sim = pd.concat(sim_ob_dfs, ignore_index=True)

# Only retain obs with DHR between 0 and -1 to prevent double-counting
# UNLESS ob is from GPSIPW, in which case keep all obs b/c DHR is always -1 for these obs
bufr_df_real = bufr_df_real.loc[np.logical_or(np.logical_and(bufr_df_real['DHR'] > -1, bufr_df_real['DHR'] <= 0),
                                              bufr_df_real['subset'] == 'GPSIPW')]
bufr_df_sim = bufr_df_sim.loc[np.logical_or(np.logical_and(bufr_df_sim['DHR'] > -1, bufr_df_sim['DHR'] <= 0),
                                            bufr_df_sim['subset'] == 'GPSIPW')]

bufr_df_real.reset_index(inplace=True)
bufr_df_sim.reset_index(inplace=True)

# Apply rounding so precision in simulated obs matches real obs
bufr_df_sim = bufr.match_bufr_prec(bufr_df_sim)
bufr_df_real = bufr.match_bufr_prec(bufr_df_real)

# Only retain obs from desired subset
boo = np.zeros(len(bufr_df_sim))
for s in subsets:
    boo[bufr_df_sim['subset'] == s] = 1
ind = np.where(boo)
bufr_df_sim = bufr_df_sim.loc[ind]
bufr_df_real = bufr_df_real.loc[ind]

for v in obs_vars:
    print('Plotting %s' % v)
    fig = plt.figure(figsize=(12, 8))
    plt.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.95, wspace=0.3)

    # Only plot if the quality marker is <= 2
    if v in qm.keys():
        cond = np.logical_and(bufr_df_sim[qm[v]] <= 2, bufr_df_real[qm[v]] <= 2)
        field1 = bufr_df_sim.loc[cond, v]
        field2 = bufr_df_real.loc[cond, v]
        diff = (field1 - field2).values
        lat = bufr_df_sim.loc[cond, 'YOB'].values
        lon = bufr_df_sim.loc[cond, 'XOB'].values
    else:
        field1 = bufr_df_sim[v]
        field2 = bufr_df_real[v]
        diff = (field1 - field2).values
        lat = bufr_df_sim['YOB'].values
        lon = bufr_df_sim['XOB'].values

    # Plot actual values
    maxval = max(np.amax(field1), np.amax(field2))
    minval = min(np.amin(field1), np.amin(field2))
    axlist = []
    for j, (f, ttl) in enumerate(zip([field1, field2], [name1, name2])):
        ax1 = fig.add_subplot(2, 3, j+1, projection=ccrs.PlateCarree())
        cax = ax1.scatter(lon, lat, s=2, c=f, cmap='plasma', vmin=minval, vmax=maxval,
                          transform=ccrs.PlateCarree())
        ax1.coastlines('50m')
        ax1.add_feature(borders, linewidth=0.5, edgecolor='k')
        ax1.set_xlim([-135, -50])
        ax1.set_ylim([20, 55])
        ax1.set_title(ttl, size=18)
        axlist.append(ax1)

        ax2 = fig.add_subplot(2, 3, j+4)
        ax2.hist(f, bins=40, range=(minval, maxval))
        ax2.set_xlabel('%s (%s)' % (v, meta[v]['units']), size=12)
        ax2.set_ylabel('counts', size=12)
        ax2.grid()

    cbar = plt.colorbar(cax, ax=axlist, orientation='horizontal')
    cbar.set_label('%s (%s)' % (v, meta[v]['units']), size=12)

    # Plot differences
    dlim = np.percentile(np.abs(diff), 99)
    ax1 = fig.add_subplot(2, 3, 3, projection=ccrs.PlateCarree())
    sort_idx = np.argsort(np.abs(diff))
    cax = ax1.scatter(lon[sort_idx], lat[sort_idx], s=2, c=diff[sort_idx], cmap='bwr', 
                      vmin=-dlim, vmax=dlim, transform=ccrs.PlateCarree())
    cbar = plt.colorbar(cax, ax=ax1, orientation='horizontal')
    cbar.set_label('%s diffs (%s)' % (v, meta[v]['units']), size=14)

    ax1.coastlines('50m')
    ax1.add_feature(borders, linewidth=0.5, edgecolor='k')
    ax1.set_xlim([-135, -50])
    ax1.set_ylim([20, 55])
    ax1.set_title('%s $-$ %s' % (name1, name2), size=18)
        
    ax2 = fig.add_subplot(2, 3, 6)
    ax2.hist(diff, bins=40, range=(-dlim, dlim))
    ax2.set_xlabel('%s diff (%s)' % (v, meta[v]['units']), size=12)
    ax2.set_ylabel('counts', size=12)
    ax2.grid()

    # Include some statistics
    plt.suptitle('Mean Diff = %.3e, Std Dev = %.3e, RMSD = %.3e, n = %d' % 
                 (np.mean(diff), np.std(diff), np.sqrt(np.mean(diff**2)), 
                  len(diff) - np.isnan(diff).sum()), size=16)

    plt.savefig(save_fname % (bufr_tag, v, date_range[0].strftime('%Y%m%d%H'), 
                              date_range[-1].strftime('%Y%m%d%H')))
    plt.close()


"""
End plot_ob_diffs_2d.py  
"""
