"""
Determine Locations of UAS Observation Sites

List of CartoPy map projections: https://scitools.org.uk/cartopy/docs/latest/reference/projections.html
Primer on map projections: https://www.e-education.psu.edu/geog160/node/1918

shawn.s.murdzek@noaa.gov
Date Created: 28 March 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import matplotlib.pyplot as plt


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Spacing between UAS sites (km)
dx = 35.

# Lower-left and upper-right corners of domain (deg N, deg E)
ll_lat =
ll_lon =
ur_lat =
ur_lon =

# Map projection
proj = ccrs.LambertConformal()

# Output text file to dump UAS site (lat, lon) coordinates
out_file = 'uas_site_locs.txt'

# Options for plotting UAS sites
make_plot = True


#---------------------------------------------------------------------------------------------------
# Compute (lat, lon) coordinates of UAS Sites
#---------------------------------------------------------------------------------------------------

# Can use proj.transform_points() to convert (x, y) coordinates to (lat, lon) coordinates (in m, I 
# think)

# The trick will be how to convert the bounds (ll_lat, ll_lon, ur_lat, ur_lon) to Cartesian points.
# There doesn't appear to be a function that does that. This will be needed in order to create the
# initial arrays of (x, y) coordinates

# Also need a shapefile containing the land mask in order to remove sites that occur over lakes and
# oceans.

#---------------------------------------------------------------------------------------------------
# Make Plot
#---------------------------------------------------------------------------------------------------



"""
End uas_sites.py
"""
