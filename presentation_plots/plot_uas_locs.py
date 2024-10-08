"""
Plot UAS Observation Locations for Two Configurations

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

file1 = ['/work2/noaa/wrfruc/murdzek/src/osse_ob_creator/fix_data/uas_site_locs_300km.txt']
file2 = ['/work2/noaa/wrfruc/murdzek/src/osse_ob_creator/fix_data/uas_site_locs_35km.txt1',
         '/work2/noaa/wrfruc/murdzek/src/osse_ob_creator/fix_data/uas_site_locs_35km.txt2',
         '/work2/noaa/wrfruc/murdzek/src/osse_ob_creator/fix_data/uas_site_locs_35km.txt3']

title1 = '300-km Spacing (84 sites)'
title2 = '35-km Spacing (6335 sites)'

ms1 = 8
ms2 = 1

save_fname = 'uas_locs.png'


#---------------------------------------------------------------------------------------------------
# Create Plots
#---------------------------------------------------------------------------------------------------

fig = plt.figure(figsize=(10, 2.75))

for i, (all_f, ttl, ms) in enumerate(zip([file1, file2], [title1, title2], [ms1, ms2])):

    # Read in UAS locations
    all_df = []
    for f in all_f:
        all_df.append(pd.read_csv(f))
    locs_df = pd.concat(all_df, ignore_index=True)
   
    print(f"number of UAS sites = {len(locs_df)}")

    # Plot UAS locations
    ax = fig.add_subplot(1, 2, i+1, projection=ccrs.LambertConformal())
    ax.plot(locs_df['lon (deg E)'].values, locs_df['lat (deg N)'].values, 
            ms=ms, c='r', lw=0, marker='.', transform=ccrs.PlateCarree())

    # Add coastlines and state boundaries
    scale = '10m'
    ax.coastlines(scale)
    borders = cfeature.NaturalEarthFeature(category='cultural',
                                           scale=scale,
                                           facecolor='none',
                                           name='admin_1_states_provinces')
    ax.add_feature(borders)
    lakes = cfeature.NaturalEarthFeature(category='physical',
                                           scale=scale,
                                           facecolor='none',
                                           name='lakes')
    ax.add_feature(lakes)
    ax.set_extent([-127, -65, 22, 49])

    # Add title
    ax.set_title(ttl, size=20)

plt.subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.92, wspace=0.02)
if save_fname[-3:] == 'png':
    plt.savefig(save_fname, dpi=750)
else:
    plt.savefig(save_fname)


"""
End plot_uas_locs.py
"""
