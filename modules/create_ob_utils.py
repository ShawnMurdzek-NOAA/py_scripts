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

def determine_twgt(wrf_hr, dhr):
    """
    Compute the weight and index needed for time interpolation using 1-D linear interpolation

    Inputs
    ------
    wrf_hr : array
        Decimal hours for each WRF UPP file
    dhr : float
        Decimal hour for the given observation

    Returns
    -------
    ihr : integer
        Index for wrf_hr associated with twgt 
    twgt : float
        Interpolation weight associated with ihr

    """

    ihr = np.where((wrf_hr - dhr) <= 0)[0][-1]
    twgt = (wrf_hr[ihr+1] - dhr) / (wrf_hr[ihr+1] - wrf_hr[ihr])

    return ihr, twgt


def _bilinear_interp_horiz(field, iwgt, jwgt, i0, j0, threeD=False):
    """
    Interpolate the given field in the horizontal

    Inputs
    ------
    field : array 
        2D field to interpolate
    iwgt : float
        Interpolation weight for i0
    jwgt : float
        Interpolation weight fro j0
    i0 : integer
        Lower-left index for the first dimension of field
    j0 : integer
        Lower-left index for the second dimension of field
    threeD : boolean, optional
        Is this field actually 3D?

    Returns
    -------
    val : float
        Interpolated value

    """

    iind = np.array([i0,        i0,            i0+1,          i0+1])
    jind = np.array([j0,        j0+1,          j0,            j0+1])
    wgts = np.array([iwgt*jwgt, iwgt*(1-jwgt), (1-iwgt)*jwgt, (1-iwgt)*(1-jwgt)])

    if threeD:
        val = np.zeros(field.shape[0])
        for w, i, j in zip(wgts, iind, jind):
            val = val + w * field[:, i, j]
    else:
        val = 0.
        for w, i, j in zip(wgts, iind, jind):
            val = val + w * field[i, j]

    return val


def _log_interp(val1, val2, wgt1):

    return (val1**wgt1) * (val2**(1.-wgt1))


def _interp_x_y_t(field1, field2, iwgt, jwgt, twgt, i0, j0, threeD=False):
    """
    Interpolate linearly in the horizontal and time dimensions

    Inputs
    ------
    field1 : array
        2D field to interpolate with time weight = twgt
    field2 : array
        2D field to interpolate with time weight = 1 - twgt
    iwgt, jwgt, twgt : float
        Interpolation weights for the (i0, j0) gridpoint in field1
    i0, j0 : integers
        Lower-left index for field1 and field2
    threeD : boolean, optional
        Is this field actually 3D?

    Returns
    -------
    val : float
        Interpolated value

    """

    val = (twgt * _bilinear_interp_horiz(field1, iwgt, jwgt, i0, j0, threeD=threeD) +
           (1.-twgt) * _bilinear_interp_horiz(field2, iwgt, jwgt, i0, j0, threeD=threeD))

    return val


def _interp_x_y_z(field, iwgt, jwgt, pwgt, i0, j0, pi0):
    """
    Interpolate linearly in the horizontal dimensions and logarithmically in pressure

    Inputs
    ------
    field : array
        3D field to interpolate
    iwgt, jwgt, pwgt : float
        Interpolation weights for the (pi0, i0, j0) gridpoint in field
    i0, j0, pi0 : integers
        Lower-left index for field

    Returns
    -------
    val : float
        Interpolated value

    """

    val = _log_interp(_bilinear_interp_horiz(field[pi0, :, :], iwgt, jwgt, i0, j0),
                      _bilinear_interp_horiz(field[pi0+1, :, :], iwgt, jwgt, i0, j0), pwgt)

    return val


def interp_wrf_to_obs(wrf_data, wrf_hr, var, ob_subset, ihr, twgt, threeD=False):
    """
    Wrapper function for interpolation in x, y, and t

    Inputs
    ------
    wrf_data : dictionary
        Dictionary of Xarray datasets containing UPP output data
    wrf_hr : list
        List of valid times (in decimal hours) for wrf_data. These are the keys for wrf_data
    var : string
        Variable from the UPP dataset to interpolate
    ob_subset : series
        Pandas series containing a single observation
    ihr : integer
        Index for the UPP dataset immediately before the observation time
    twgt : float
        Interpolation weight associated with the ihr UPP dataset
    threeD : boolean, optional
        Is this UPP field 3D?

    Returns
    -------
    val : float
        UPP output linearly interpolated in x, y, and t to the observation location

    """

    val = _interp_x_y_t(wrf_data[wrf_hr[ihr]][var], wrf_data[wrf_hr[ihr+1]][var], twgt,
                        ob_subset['iwgt'], ob_subset['jwgt'], ob_subset['i0'], ob_subset['j0'],
                        threeD=threeD)

    return val


def interp_wrf_3d(wrf3d, ob_subset):
    """
    Wrapper function for interpolation in x, y, and p

    Inputs
    ------
    wrf3d : array
        WRF 3D output array to interpolate
    ob_subset : series
        Pandas series containing a single observation

    Returns
    -------
    val : float
        UPP output linearly interpolated in x, y, and p to the observation location

    """

    val = _interp_x_y_z(wrf3d, ob_subset['iwgt'], ob_subset['jwgt'], ob_subset['pwgt'],
                        ob_subset['i0'], ob_subset['j0'], ob_subset['pi0'])

    return val


def interp_wrf_p1d(p1d, ob_subset):
    """
    Wrapper function for interpolation of pressure in one dimension

    Inputs
    ------
    p1d : array
        1D pressure array to interpolate
    ob_subset : series
        Pandas series containing a single observation

    Returns
    -------
    val : float
        UPP output logarithmically interpolated in p to the observation location
    pwgt : float
        Weight used to logarithmic interpolation

    """

    pi0 = ob_subset['pi0']
    pwgt = (p1d[pi0+1] - ob_subset['POB']) / (p1d[pi0+1] - p1d[pi0])
    val = _log_interp(p1d[pi0], p1d[pi0+1], pwgt)

    return val, pwgt


"""
End create_ob_utils.py
"""
