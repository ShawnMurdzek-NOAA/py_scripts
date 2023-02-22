"""
Functions for Map Projections

This code is based on the map_utils module in WRF, which can be found in the following WPS file:

WPS/geogrid/src/module_map_utils.F

Only the Lambert Conformal map projection (which is used for the HRRR) is currently supported.

Even though this code is essentially copied from WPS, it still cannot perfectly reproduce the map 
projection from WRF. For example, if a certain (lat, lon) coordinate is associated with a certain
(i, j) index, plugging that (lat, lon) coordinate into ll_to_xy_lc does not yield the same (i, j)
index. This is true even after switching to single precision from double precision (WRF uses single
precision, I believe). The snippet of code after the ll_to_xy_lc function graphically shows the
difference between the actual gridpoint indices from WRF and those computed using ll_to_xy_lc.

shawn.s.murdzek@noaa.gov
Date Created: 22 February 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np


#---------------------------------------------------------------------------------------------------
# Functions
#---------------------------------------------------------------------------------------------------

def ll_to_xy_lc(lat, lon, ref_lat=38.497246, ref_lon=-97.505974, truelat1=38.5, truelat2=38.5, 
                stand_lon=-97.5, dx=1., e_we=5400, e_sn=3180, knowni=2699, knownj=1589):
    """
    Convert (lat, lon) coordinates to (x, y) coordinates usig a Lambert Conformal map projection

    Parameters
    ----------
    lat : 1-D array
        Array of latitude coordinates in the range -90 to 90 deg N
    lon : 1-D array
        Array of longitude coordinates in the range -180 to 180 deg E
    ref_lat : float, optional
        Latitude coordinate of grid center (in range -90 to 90 deg N). Aka lat1.
    ref_lon : float, optional
        Longitude coordinate of grid center (in range -180 to 180 deg E). Aka lon1.
    truelat1 : float, optional
        First true latitude (in range -90 to 90 deg N)
    truelat2 : float, optional
        Second true latitude (in range -90 to 90 deg N). Can be the same as truelat1
    stand_lon : float, optional
        The longitude parallel to the y-axis (in range -180 to 180 deg E)
    dx : float, optional
        Grid spacing (km)
    e_we : integer, optional
        Number of points in x direction
    e_sn : integer, optional
        Number of points in y direction
    knowni : integer, optional
        Index of (ref_lon, ref_lat). If NaN, determine knowni using e_we
    knownj : integer, optional
        Index of (ref_lon, ref_lat). If NaN, determine knownj using e_sn

    Returns
    -------
    xi : array
        Array of x coordinates
    yi : array
        Array of y coordinates

    Notes
    -----
    The (x, y) "coordinates" are the decimal indices for the given (lat, lon) coordinates for the 
    specified map projection. To get coordinates in km, set dx = 1. The origin is assumed to be the
    SW corner.

    If knowni and knownj are specified, ref_lat and ref_lon do not necessarily need to be the center 
    of the grid as long as this coordinate aligns with knowni and knownj

    Default values for optional parameters come from the values used in the OSSE Nature Run (which
    is identical to the HRRR, except the gridspacing is reduced from 3 to 1 km) 

    Code is adapted from the llij_lc() subroutine in module_map_utils.F

    """

    # Set radius of the Earth (value comes from constants_module.F in WPS)
    re = 6370. 
    rebydx = re / dx

    # Compute deltalon
    deltalon = lon - stand_lon
    ind1 = np.where(deltalon > 180.)[0]
    ind2 = np.where(deltalon < -180.)[0]
    if len(ind1) > 0:
        deltalon[ind1] = deltalon[ind1] - 360.
    if len(ind2) > 0:
        deltalon[ind2] = deltalon[ind2] + 360.

    # Compute cosine of truelat1
    tl1r = np.deg2rad(truelat1)
    ctl1r = np.cos(tl1r)

    # Set hemisphere factor
    if truelat1 < 0:
        hemi = -1.
    else:
        hemi = 1.

    # Determine center gridpoint indices (from metgrid/src/process_domain_module.F)
    if np.isnan(knowni):
        if e_we % 2 == 0:
            knowni = 0.5 * (e_we + 1) - 1
        else:
            knowni = 0.5 * e_we - 1
        if e_sn % 2 == 0:
            knownj = 0.5 * (e_sn + 1) - 1
        else:
            knownj = 0.5 * e_sn - 1

    # Compute cone factor (round to 7 decimal places to match WRF)
    if np.isclose(truelat1, truelat2):
        cone = np.sin(np.abs(np.deg2rad(truelat1)))
    else:
        cone = np.log10(np.cos(np.deg2rad(truelat1))) - np.log10(np.cos(np.deg2rad(truelat2)))
        cone = cone / (np.log10(np.tan(np.deg2rad(45.-np.abs(truelat1)/2.))) -
                       np.log10(np.tan(np.deg2rad(45.-np.abs(truelat2)/2.))))

    # Find pole point
    deltalon1 = ref_lon - stand_lon
    if deltalon1 > 180:
        deltalon1 = deltalon1 - 360.
    elif deltalon1 < -180:
        deltalon1 = deltalon1 + 360.
    rsw = rebydx * (ctl1r/cone) * (np.tan(np.deg2rad(90.*hemi-ref_lat) / 2.) /
                                   np.tan(np.deg2rad(90.*hemi-truelat1) / 2.)) ** cone
    arg = cone * np.deg2rad(deltalon1)
    polei = hemi * knowni - (hemi * rsw * np.sin(arg))
    polej = hemi * knownj + (rsw * np.cos(arg))

    # Determine radii to each of the desired points
    rm = rebydx * (ctl1r/cone) * (np.tan(np.deg2rad(90.*hemi-lat) / 2.) /
                                  np.tan(np.deg2rad(90.*hemi-truelat1) / 2.)) ** cone 

    # Compute (xi, yi) values for each point (subtract 1 b/c Python uses 0 to represent the first 
    # index)
    arg = cone * np.deg2rad(deltalon)   
    xi = polei + (hemi * rm * np.sin(arg))
    yi = polej - (rm * np.cos(arg))

    # Flip xi, yi values in Southern Hemisphere
    xi = xi * hemi
    yi = yi * hemi

    return xi, yi


#---------------------------------------------------------------------------------------------------
# Plot Showing Differences Between WRF Gridpoint Locations and Calculated Gridpoint Locations
#---------------------------------------------------------------------------------------------------
'''
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# Extract data
fname = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_spring_v2/wrfnat_202204291200.grib2'
ds = xr.open_dataset(fname, engine='pynio')

# Compute expected gridpoint locations
lat = ds['gridlat_0'].values
lon = ds['gridlon_0'].values
shape = lat.shape
xi, yi = ll_to_xy_lc(lat.ravel(), lon.ravel())
xi2d = np.reshape(xi, shape)
yi2d = np.reshape(yi, shape)

# Plot differences
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

xit, yit = np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
diff = np.sqrt((xi2d - xit)**2 + (yi2d - yit)**2)

cax = ax.contourf(lon, lat, diff, np.linspace(0, np.amax(diff), 20), transform=ccrs.PlateCarree())
cbar = plt.colorbar(cax, ax=ax, orientation='horizontal')
cbar.set_label('diff (km)', size=14)

ax.set_title('MAE = %.4f km' % np.mean(diff), size=16)
ax.coastlines('50m')

plt.savefig('gridpoint_diff.png')
plt.close()
'''

"""
End map_proj.py
"""
