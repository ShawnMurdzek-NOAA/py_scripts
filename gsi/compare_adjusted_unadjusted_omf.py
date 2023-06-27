"""
Quick Script to Compare Unadjusted and Adjusted O-B and O-A Values from GSI Diag Files

shawn.s.murdzek@noaa.gov
Date Created: 19 May 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import pandas as pd
import numpy as np
import glob

import pyDA_utils.gsi_fcts as gsi


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

rrfs_dir = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/real_data/winter/NCO_dirs/ptmp/prod/rrfs.20220201/12'


#---------------------------------------------------------------------------------------------------
# Compare Adjusted and Unadjusted O-A and O-B for Each Ob Type
#---------------------------------------------------------------------------------------------------

fnames = glob.glob('%s/*diag*' % rrfs_dir)

for f in fnames:
    print()
    print(f)
    diag = gsi.read_diag([f])
    ob_types = np.unique(diag['Observation_Type'].values)
    for t in ob_types:
        cond = (diag['Observation_Type'] == t)
        if diag['var'].iloc[0] == 'uv':
            rmsd = np.sqrt(np.mean((diag['u_Obs_Minus_Forecast_adjusted'].loc[cond] - 
                                    diag['u_Obs_Minus_Forecast_unadjusted'].loc[cond])**2))
        else:
            rmsd = np.sqrt(np.mean((diag['Obs_Minus_Forecast_adjusted'].loc[cond] - 
                                    diag['Obs_Minus_Forecast_unadjusted'].loc[cond])**2))
        print('%d   %.2e' % (t, rmsd))


"""
End compare_adjusted_unadjusted_omf.py
"""
