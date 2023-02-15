"""
General Meteorology Utilities

shawn.s.murdzek@noaa.gov
Date Created: 15 February 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np


#---------------------------------------------------------------------------------------------------
# Define Functions
#---------------------------------------------------------------------------------------------------

def equil_vapor_prs(T):
    """
    Equilibrium vapor pressure with respect to liquid water

    Parameters
    ----------
    T : float
        Temperature (K)

    Returns
    -------
    es : float
        Equilibrium vapor pressure (Pa)

    Notes
    -----
    Equation is based on that from Bolton (1980) and reproduced in Markowski and Richardson (2010),
    eqn (2.16)

    """

    T = T - 273.15
    es = 611.2 * np.exp((17.67 * T) / (T + 243.5))

    return es


def equil_mix(T, p):
    """
    Equilibrium mixing ratio with respect to liquid water

    Parameters
    ----------
    T : float
        Temperature (K)
    p : float
        Pressure (Pa)

    Returns
    -------
    qvs : float
        Equilibrium mixing ratio (kg / kg)

    Notes
    -----
    Equation is based on Markowski and Richardson (2010), eqn (2.11)

    """

    epn = 0.622
    es = equil_vapor_prs(T)
    qvs = (epn * es) / (p - es)

    return qvs


"""
End meteo_util.py
"""
