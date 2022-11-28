"""
Utilities for Creating Synthetic Observations

shawn.s.murdzek@noaa.gov
Date Created: 22 November 2022
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import pandas as pd


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


def read_ob_errors(fname):
    """
    Parse out obbservation errors from an errtable file in GSI

    Parameters
    ----------
    fname : string
        Name of errtable text file

    Returns
    -------
    errors: dictionary
        A dictionary of pd.DataFrame objects containing the observation errors

    Notes
    -----
    More information about the errtable format in GSI can be found here: 
    https://dtcenter.ucar.edu/com-GSI/users/docs/users_guide/html_v3.7/gsi_ch4.html#conventional-observation-errors

    """

    # Extract contents of file
    fptr = open(fname, 'r')
    contents = fptr.readlines()
    fptr.close()

    # Loop over each line
    errors = {}
    headers = ['prs', 'Terr', 'RHerr', 'UVerr', 'PSerr', 'PWerr']
    for l in contents:
        if l[5:21] == 'OBSERVATION TYPE':
            key = int(l[1:4])
            errors[key] = {}
            for h in headers:
                errors[key][h] = []
        else:
            vals = l.strip().split(' ')
            for k, h in enumerate(headers):
                errors[key][h].append(float(vals[k]))

    # Convert to DataFrame
    for key in errors.keys():
        errors[key] = pd.DataFrame(errors[key])
        for h in headers[1:]:
            errors[key][h].where(errors[key][h] < 5e8, inplace=True)

    return errors


"""
End create_ob_utils.py
"""
