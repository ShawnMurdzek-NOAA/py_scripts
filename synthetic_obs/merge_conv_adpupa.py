"""
Combine Conventional Obs with ADPUPA Obs That Account for Radiosonde Drift

Passed Arguments:
    argv[1] = Conventional ob prepbufr csv name
    argv[2] = ADPUPA ob prepbufr csv name
    argv[3] = Output prepbufr csv name

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

# Conventional ob prepbufr
conv_fname = ''

# ADPUPA ob prepbufr
adpupa_fname = ''

# Output prepbufr
out_fname = ''

# Used passed arguments, if they exist
if len(sys.argv) > 1:
    conv_fname = sys.argv[1]
    adpupa_fname = sys.argv[2]
    out_fname = sys.argv[3]


#---------------------------------------------------------------------------------------------------
# Combine CSV Files
#---------------------------------------------------------------------------------------------------

conv_bufr = bufr.bufrCSV(conv_fname)
adpupa_bufr = bufr.bufrCSV(adpupa_fname)

out_df = conv_bufr.df.copy()

# Drop ob types found in adpupa dataframe
conv_typ = conv_bufr.df['TYP'].unique()
adpupa_typ = adpupa_bufr.df['TYP'].unique()
for typ in conv_typ:
    if typ in adpupa_typ:
        out_df.drop(index=np.where(out_df['TYP'] == typ)[0], inplace=True)
        out_df.reset_index(drop=True, inplace=True)

out_df = pd.concat([adpupa_bufr.df, out_df])

bufr.df_to_csv(out_df, out_fname)


"""
End merge_conv_adpupa.py
"""
