"""
Combine Conventional Obs with ADPUPA Obs That Account for Radiosonde Drift

Passed Arguments:
    sys[1] = Output simulated observation directory
    sys[2] = Timestamp on prpebufr CSV (YYYYMMDDHHMM)

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

out_dir = sys.argv[1]
time = sys.argv[2]

# Conventional ob prepbufr
conv_fname = '%s/conv/%s.fake.prepbufr.csv' % (out_dir, time)

# ADPUPA ob prepbufr
adpupa_fname = '%s/adpupa/%s.fake.adpupa.csv' % (out_dir, time)

# Output prepbufr
out_fname = '%s/perfect/%s.fake.prepbufr.csv' % (out_dir, time)


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
