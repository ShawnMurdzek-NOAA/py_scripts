"""
Create a Sample Prepbufr CSV

shawn.s.murdzek@noaa.gov
Date Created: 18 October 2022
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import bufr
import numpy as np


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

fname = '/mnt/lfs4/BMC/wrfruc/murdzek/sample_real_obs/obs_rap/prepbufr.csv'
save_fname = '/mnt/lfs4/BMC/wrfruc/murdzek/sample_real_obs/obs_rap/sample_prepbufr.csv'


#---------------------------------------------------------------------------------------------------
# Create Sample PrepBUFR CSV
#---------------------------------------------------------------------------------------------------

bufr_df = bufr.bufrCSV(fname)
bufr_df.sample(save_fname)


"""
End create_sample_csv.py  
"""
