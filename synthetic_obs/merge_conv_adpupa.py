"""
Combine Conventional Obs with ADPUPA Obs That Account for Radiosonde Drift

Passed Arguments:
    sys[1] = Timestamp on prpebufr CSV (YYYYMMDDHHMM)

shawn.s.murdzek@noaa.gov
Date Created: 21 March 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import pandas as pd
import bufr
import sys


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

work = '/work2/noaa/'

# Conventional ob prepbufr
conv_fname = work + 'wrfruc/murdzek/nature_run_spring/synthetic_obs_csv/conv/%s.fake.prepbufr.csv' % sys.argv[1]

# ADPUPA ob prepbufr
adpupa_fname = work + 'wrfruc/murdzek/nature_run_spring/synthetic_obs_csv/adpupa/%s.fake.adpupa.csv' % sys.argv[1]

# Output prepbufr
out_fname = work + 'wrfruc/murdzek/nature_run_spring/synthetic_obs_csv/perfect/%s.fake.prepbufr.csv' % sys.argv[1]


#---------------------------------------------------------------------------------------------------
# Combine CSV Files
#---------------------------------------------------------------------------------------------------

conv_bufr = bufr.bufrCSV(conv_fname)
adpupa_bufr = bufr.bufrCSV(adpupa_fname)

out_df = conv_bufr.df.loc[conv_bufr.df['subset'] != 'ADPUPA'].copy()
out_df = pd.concat([adpupa_bufr.df, out_df])

bufr.df_to_csv(out_df, out_fname)


"""
End merge_conv_adpupa.py
"""
