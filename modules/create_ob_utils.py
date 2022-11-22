"""
Utilities for Creating Synthetic Observations

shawn.s.murdzek@noaa.gov
Date Created: 22 November 2022
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np


#---------------------------------------------------------------------------------------------------
# Functions
#---------------------------------------------------------------------------------------------------

def wrf_coords(lat, lon, ds):
    """
    Find the decimal (i.e., non-integer) x and y indices for a given (lat, lon) coordinate in WRF
 
    Parameters
    ----------
    lat, lon : float
        Latitude and longitude coordinate (deg E and deg W)
    ds : XArray.dataset
        Dataset containing the WRF output data

    Returns
    -------
    xi, yi : float
        Indices in the west-east and south-north directions for the input (lat, lon) coordinate

    Notes
    -----
    Timing tests indicate that one call to this function takes ~0.025 s.

    """

    # Find nearest neighbor in horizontal direction
    close = np.unravel_index(np.argmin((ds['XLONG'] - lon)**2 + (ds['XLAT'] - lat)**2), ds['XLONG'].shape)
    clat = ds['XLAT'][close[0], close[1]].values
    clon = ds['XLONG'][close[0], close[1]].values

    # Perform pseudo-bilinear interpolation
    if lat < clat:
        yi = close[0] + (lat - clat) / (clat - ds['XLAT'][close[0]-1, close[1]].values)
    else:
        yi = close[0] + (lat - clat) / (ds['XLAT'][close[0]+1, close[1]].values - clat)
    if lon < clon:
        xi = close[1] + (lon - clon) / (clon - ds['XLONG'][close[0], close[1]-1].values)
    else:
        xi = close[1] + (lon - clon) / (ds['XLONG'][close[0], close[1]+1].values - clon)

    return xi, yi


"""
End create_ob_utils.py
"""
