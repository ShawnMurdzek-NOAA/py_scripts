"""
Plot Observation Locations and Model Gridpoint Locations

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import pyproj

from pyDA_utils import bufr


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

bufr_csv_fname = '/work2/noaa/wrfruc/murdzek/nature_run_spring/obs/perfect_conv/real_csv/202204302000.rap.prepbufr.csv'
obtypes = list(range(180, 190))

# Model grid spacing (m)
model_dx = 13000

# Plotting projection and scale for features (e.g., coastlines)
plot_proj = ccrs.LambertConformal()
scale = '10m'
figsize = (8, 8)
markersize = 20

# Plotting domain
lat_lim = [37.5, 39.75]
lon_lim = [-78, -75]

#---------------------------------------------------------------------------------------------------
# Create Plots
#---------------------------------------------------------------------------------------------------

# Plot observation locations
fig = plt.figure(figsize=figsize)
bufr_obj = bufr.bufrCSV(bufr_csv_fname)
bufr_obj.select_obtypes(obtypes)
ax = bufr.plot_obs(bufr_obj.df, fig=fig, proj=plot_proj, scale=scale, s=markersize, c='b')
ax.set_extent([lon_lim[0], lon_lim[1], lat_lim[0], lat_lim[1]])
ax.set_title('Surface Observations', size=24)
plt.savefig('sfc_obs.png')

# Create model gridpoints in (lat, lon) space
proj_str = '+proj=lcc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45'
proj = pyproj.Proj(proj_str)
npts_we = 23250000 / model_dx
npts_sn = 13650000 / model_dx
x_pt, y_pt = np.meshgrid(np.arange(-0.5*npts_we*model_dx, 0.5*npts_we*model_dx + (0.1*model_dx), model_dx),
                         np.arange(-0.5*npts_sn*model_dx, 0.5*npts_sn*model_dx + (0.1*model_dx), model_dx))
lon_pt, lat_pt = proj(x_pt.ravel(), y_pt.ravel(), inverse=True)

# Plot model gridpoints
fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(1, 1, 1, projection=plot_proj)
ax.coastlines(scale)
borders = cfeature.NaturalEarthFeature(category='cultural',
                                       scale=scale,
                                       facecolor='none',
                                       name='admin_1_states_provinces')
ax.add_feature(borders)
ax.scatter(lon_pt, lat_pt, transform=ccrs.PlateCarree(), s=markersize, c='g')
ax.set_extent([lon_lim[0], lon_lim[1], lat_lim[0], lat_lim[1]])
ax.set_title('Model Gridpoints ({dx}-km Grid)'.format(dx=int(1e-3*model_dx)), size=24)
plt.savefig('model_gridpoints_green.png')


"""
End plot_ob_locs_and_model_gridpts.py
"""
