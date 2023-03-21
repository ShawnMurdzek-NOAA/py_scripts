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


def _linear_interp(val1, val2, wgt):
    """
    Interpolate linearly between two values given a specific weight

    Inputs
    ------
    val1 : float or array
        Value(s) with weight 'wgt'
    val2 : float or array
        Value(s) with weight '1-wgt'
    wgt : float
        Interpolation weight

    Returns
    -------
    val : float
        Interpolated value

    """

    return (val1 * wgt) + (val2 * (1. - wgt))


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


def interp_x_y(field, ob_subset, threeD=False):
    """
    Wrapper function for interpolation in x and y

    Inputs
    ------
    field : array
        Array containing WRF data to be interpolated
    ob_subset : series
        Pandas series containing a single observation
    threeD : boolean, optional
        Is this UPP field 3D?

    Returns
    -------
    val : float
        UPP output linearly interpolated in x and y to the observation location

    """

    iwgt = ob_subset['iwgt']
    jwgt = ob_subset['jwgt']
    i0 = ob_subset['i0']
    j0 = ob_subset['j0']

    val = _bilinear_interp_horiz(field, iwgt, jwgt, i0, j0, threeD=threeD)

    return val


def interp_x_y_t(wrf_data, wrf_hr, var, ob_subset, ihr, twgt, threeD=False):
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

    field1 = wrf_data[wrf_hr[ihr]][var]
    field2 = wrf_data[wrf_hr[ihr+1]][var]
    iwgt = ob_subset['iwgt']
    jwgt = ob_subset['jwgt']
    i0 = ob_subset['i0']
    j0 = ob_subset['j0']

    val = _linear_interp(_bilinear_interp_horiz(field1, iwgt, jwgt, i0, j0, threeD=threeD),
                         _bilinear_interp_horiz(field2, iwgt, jwgt, i0, j0, threeD=threeD), twgt)

    return val


def interp_x_y_z(wrf3d, ob_subset):
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
    
    iwgt = ob_subset['iwgt']
    jwgt = ob_subset['jwgt']
    pwgt = ob_subset['pwgt']
    i0 = ob_subset['i0']
    j0 = ob_subset['j0']
    pi0 = ob_subset['pi0']

    val = _linear_interp(_bilinear_interp_horiz(wrf3d[pi0, :, :], iwgt, jwgt, i0, j0),
                         _bilinear_interp_horiz(wrf3d[pi0+1, :, :], iwgt, jwgt, i0, j0), pwgt)

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
        Weight used for interpolation

    """

    pi0 = ob_subset['pi0']
    pwgt = np.log10(p1d[pi0+1] / ob_subset['POB']) / np.log10(p1d[pi0+1] / p1d[pi0])
    val = _linear_interp(p1d[pi0], p1d[pi0+1], pwgt)

    return val, pwgt


"""
End create_ob_utils.py
"""
