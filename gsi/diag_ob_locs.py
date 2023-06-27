"""
Plot Ob Locations for all Types for a Single Variable

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
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import pyDA_utils.gsi_fcts as gsi


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# O-Bs are found in the "ges" files and O-As are found in the "anl" files
# Can 1 or 2 datasets. Key is the name of the dataset
tmpl = '/work2/noaa/wrfruc/murdzek/tmp/diag_conv_q_ges.%s.nc4'
dates = [dt.datetime(2023, 6, 21, 17) ]

fnames = {}
fnames['real'] = [tmpl % d.strftime('%Y%m%d%H') for d in dates]

# Variable to plot (only necessary for the uv diag files)
diag_var = 'u'

# Subset of each observation type to plot ('all' - all obs, 'assim' - only obs that are assimilated)
data_subset = 'assim'

# Output directory and string to add to output file names
out_dir = './'
out_str = ''


#---------------------------------------------------------------------------------------------------
# Compute Statistics
#---------------------------------------------------------------------------------------------------

data_names = list(fnames.keys())
if fnames[data_names[0]][0].split('_')[-2] == 'uv':
    vname = '%s_Obs_Minus_Forecast_adjusted' % diag_var
else:
    vname = 'Obs_Minus_Forecast_adjusted'
    diag_var = fnames[data_names[0]][0].split('_')[-2]

# Extract data
diag_df = {}
for key in data_names:
    print('Dataset = %s' % key)
    diag_df[key] = gsi.read_diag(fnames[key])
    if data_subset == 'assim':
        diag_df[key] = diag_df[key].loc[diag_df[key]['Analysis_Use_Flag'] == 1]
diag_dates = np.unique(diag_df[data_names[0]]['date_time'])

# Create list of ob types
ob_typ = np.unique(diag_df[data_names[0]]['Observation_Type'].values)
if len(data_names) > 1:
    for key in data_names[1:]:
        ob_typ = np.unique(np.concatenate([ob_typ, diag_df[key]['Observation_Type'].values]))

# Plot data
for typ in ob_typ:
    fig = plt.figure(figsize=(10, 10))

    for i, key in enumerate(data_names):
        ax = fig.add_subplot(i+1, 1, len(data_names), projection=ccrs.PlateCarree()) 
        lat = diag_df[key].loc[diag_df[key]['Observation_Type'] == typ, 'Latitude'].values
        lon = diag_df[key].loc[diag_df[key]['Observation_Type'] == typ, 'Longitude'].values
        ax.plot(lon, lat, 'b.', transform=ccrs.PlateCarree())
        ax.coastlines('50m')
        borders = cfeature.NaturalEarthFeature(category='cultural',
                                               scale='50m',
                                               facecolor='none',
                                               name='admin_1_states_provinces')
        ax.add_feature(borders)
        ax.set_extent([230, 295, 23, 47])

    # Add additional metrics
    ttl = ('ob type = %d\nstart = %d, end = %d, var = %s, subset = %s' % 
           (typ, np.amin(diag_dates), np.amax(diag_dates), diag_var, data_subset)) 

    plt.suptitle(ttl, size=18) 

    plt.savefig('%s/loc_ob%s_%s_%s_%s_%d_%d.png' % 
                (out_dir, typ, out_str, diag_var, data_subset, np.amin(diag_dates), np.amax(diag_dates)))
    plt.close() 


"""
End diag_ob_locs.py
"""
