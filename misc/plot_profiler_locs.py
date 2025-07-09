"""
Plot Profiler Locations from prepBUFR CSV Files

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import pandas as pd
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from pyDA_utils import bufr


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# BUFR files
bufr_path = '/work/noaa/wrfruc/murdzek/nature_run_spring/obs/perfect_conv/real_csv'
bufr_start = dt.datetime(2022, 4, 29, 12)
bufr_files = [f"{bufr_path}/{(bufr_start + dt.timedelta(hours=i)).strftime('%Y%m%d%H%M')}.rap.prepbufr.csv"
              for i in range(169)]
#bufr_files = [f"{bufr_path}/{(bufr_start + dt.timedelta(hours=i)).strftime('%Y%m%d%H%M')}.rap.prepbufr.csv"
#              for i in range(3)]

# BUFR observation types
bufr_types = [126, 223, 224, 227, 228, 229]

# Output file (include {t} placeholder for type)
out_fname = "{t}_spatial_dist.png"


#---------------------------------------------------------------------------------------------------
# Main Program
#---------------------------------------------------------------------------------------------------

# Open prepBUFR CSV files and combine into a single DataFrame
# It takes ~20 s to read each CSV file
all_df = []
for f in bufr_files:
    print(f"Reading {f}")
    bufr_obj = bufr.bufrCSV(f)
    bufr_obj.select_obtypes(bufr_types)
    all_df.append(bufr_obj.df)
ob_df = pd.concat(all_df)

# Create CartoPy feature with state and lake borders
borders = cfeature.NaturalEarthFeature(category='cultural',
                                       scale='50m',
                                       facecolor='none',
                                       name='admin_1_states_provinces')
lakes = cfeature.NaturalEarthFeature(category='physical',
                                     scale='50m',
                                     facecolor='none',
                                     name='lakes')

print()
all_cts = {}
for t in bufr_types:
    # Create a summary dictionary with ob counts
    all_sid = ob_df.loc[ob_df['TYP'] == t, 'SID'].unique()
    all_cts[t] = {'SID':all_sid, 'XOB':[], 'YOB':[], 'N':[]}
    print(f"Ob type {t} has {len(all_sid)} unique SIDs")
    if len(all_sid) == 0:
        continue
    for sid in all_sid:
        subset = ob_df.loc[(ob_df['TYP'] == t) & (ob_df['SID'] == sid), :]
        all_cts[t]['XOB'].append(subset['XOB'].values[0] - 360.)
        all_cts[t]['YOB'].append(subset['YOB'].values[0])
        all_cts[t]['N'].append(len(subset['cycletime'].unique()))

    # Create plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
    plt.subplots_adjust(left=0.03, bottom=0.1, right=0.97, top=0.92)

    cax = ax.scatter(all_cts[t]['XOB'], all_cts[t]['YOB'], c=all_cts[t]['N'], cmap='plasma_r', s=16, 
                     vmin=0, vmax=len(bufr_files), transform=ccrs.PlateCarree())
    cbar = plt.colorbar(cax, ax=ax, orientation='horizontal')
    cbar.set_label('number of cycles with ob available', size=14)
    ax.set_title(f"Ob type = {t} (n = {len(all_sid)})", size=18)
    ax.coastlines('50m', linewidth=0.5, edgecolor='k')
    ax.add_feature(borders, linewidth=0.3, edgecolor='gray')
    ax.add_feature(lakes, linewidth=0.3, edgecolor='gray')

    ax.set_extent([-127, -63, 20, 53])

    plt.savefig(out_fname.format(t=t))
    plt.close()


"""
End plot_profiler_locs.py
"""
