"""
Plot Surface Station Locations

shawn.s.murdzek@noaa.gov
Date Created: 1 March 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import matplotlib.pyplot as plt


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Surface station files
station_files = []
path = '/scratch1/BMC/wrfruc/murdzek/sample_real_obs/sfc_stations'
for s in ['OAK', 'TUS', 'GJT', 'DTW', 'ABR', 'ALB', 'MIA', 'BHM', 'CRP', 'SLE']:
    station_files.append('%s/%s_202204290000_202205070000.txt' % (path, s))

# Output file
out_fname = 'sfc_station_locs.png'


#---------------------------------------------------------------------------------------------------
# Plot Surface Station Locations
#---------------------------------------------------------------------------------------------------

fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
for f in station_files:
    df = pd.read_csv(f, skiprows=5)
    ax.plot(df['lon'].loc[0], df['lat'].loc[0], 'bo', markersize=10, transform=ccrs.PlateCarree())

ax.set_extent([-130, -65, 20, 50])
ax.set_title('Surface Station Locations', size=18)
ax.coastlines('50m')
borders = cfeature.NaturalEarthFeature(category='cultural',
                                       scale='50m',
                                       facecolor='none',
                                       name='admin_1_states_provinces')
ax.add_feature(borders, linewidth=0.5, edgecolor='k')

plt.savefig(out_fname)
plt.close()


"""
End plot_sfc_station_locs.py
"""
