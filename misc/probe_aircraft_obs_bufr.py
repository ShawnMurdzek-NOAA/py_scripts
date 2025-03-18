"""
Probe Aircraft Observations from PrepBUFR files

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy

from pyDA_utils import bufr


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Observation types
all_typ = [130, 131, 133, 134, 135]

# BUFR CSV file names
bufr_fnames = [f"/work/noaa/wrfruc/murdzek/nature_run_spring/obs/perfect_conv/real_csv/20220430{i:02d}00.rap.prepbufr.csv" for i in range(23)]
#bufr_fnames = [f"/work/noaa/wrfruc/murdzek/nature_run_winter/obs/perfect_conv/real_csv/20220201{i:02d}00.rap.prepbufr.csv" for i in range(23)]


#---------------------------------------------------------------------------------------------------
# Main Program
#---------------------------------------------------------------------------------------------------

# General information
bufr_obj = []
for f in bufr_fnames:
    print()
    print(f.split('/')[-1])
    bufr_obj.append(bufr.bufrCSV(f))
    for t in all_typ:
        subset = bufr_obj[-1].df.loc[bufr_obj[-1].df['TYP'] == t, :]
        print(f"Type {t} (n = {len(subset)}), DHR min = {subset['DHR'].min()}, DHR max = {subset['DHR'].max()}")

# Determine if any observations are repeated across prepBUFR files
print()
print(50*'-')
print('Checking for duplicates...')
for i in range(1, len(bufr_obj)):
    print()
    print(bufr_fnames[i].split('/')[-1])
    copy_df = copy.deepcopy(bufr_obj[i-1].df)
    copy_df.loc[:, 'DHR'] = copy_df.loc[:, 'DHR'] - 1
    concat_df = pd.concat([copy_df, bufr_obj[i].df])
    for t in all_typ:
        subset = concat_df.loc[concat_df['TYP'] == t, ['TYP', 'SID', 'DHR', 'XOB', 'YOB', 'ELV', 'TOB', 'QOB', 'UOB', 'VOB']]
        subset.reset_index(drop=True)
        duplicates = subset.duplicated()
        print(f"Type {t} has {duplicates.sum()} duplicates")
        
        # Focus on duplicates w/ DHR > -0.75. These have the potential of being assimilated twice
        # B/c aircraft obs use a +/- 45 min window for DA (other than type 131)
        obs_DA2 = subset.loc[np.logical_and(subset['DHR'] > -0.75, duplicates), :]
        print(f"Type {t} has {len(obs_DA2)} duplicates with DHR > -0.75")
        print('Examples =', np.unique(obs_DA2['SID'].values)[:3])


"""
End probe_aircraft_obs_bufr.py
"""
