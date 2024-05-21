"""
Count Number of Restricted PrepBUFR Obs

Count the number of observations for each 3-digit BUFR observation type, then determine which mesonet
obs specifically are restricted using the RSRD field.

NOTE: RSRD must be dumped from the prepBUFR file!

Inputs
------
argv[1] : Prepbufr CSV file

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
from pyDA_utils import bufr
import matplotlib.pyplot as plt
import sys
import pandas as pd
import copy


#---------------------------------------------------------------------------------------------------
# Main Program
#---------------------------------------------------------------------------------------------------

def count_ob_typ(bufr_csv, typ=[120]):
    """
    Count the number of observations that match the given types
    """

    n = 0
    for t in typ:
        n = n + np.sum(bufr_csv.df['TYP'] == t)

    return n
    

def count_ob_rsrd(bufr_csv):
    """
    Count the number of obs that meet each restriction criteria
    """

    rsrd_ct_dict = {'TYP':[], 'RSRD':[], 'COUNT':[]}
    rsrd = np.unique(bufr_csv.df['RSRD'])
    for r in rsrd:
        if not np.isnan(r):
            rcond = bufr_csv.df['RSRD'] == r
        else:
            rcond = np.isnan(bufr_csv.df['RSRD'])
        typ = np.unique(bufr_csv.df.loc[rcond, 'TYP'])
        for t in typ:
            rsrd_ct_dict['TYP'].append(int(t))
            rsrd_ct_dict['RSRD'].append(r)
            rsrd_ct_dict['COUNT'].append(len(bufr_csv.df.loc[rcond & (bufr_csv.df['TYP'] == t)]))

    return pd.DataFrame(rsrd_ct_dict)                
    

def count_all_obs(bufr_csv, typ=list(range(100, 200))):
    """
    Count each ob type and different ob restrictions for each ob type in the given window
    """

    # Count number of restricted obs first
    rsrd_df = count_ob_rsrd(bufr_csv)
    rsrd_typ = np.unique(rsrd_df['TYP'])
    for t in rsrd_typ:
        if t not in typ:
            rsrd_df = rsrd_df.loc[rsrd_df['TYP'] != t]
        # Don't care about RSRD = 1, b/c these obs are unrestricted
        elif ((rsrd_df.loc[rsrd_df['TYP'] == t, 'RSRD'].values[0] == 1) and 
              (len(rsrd_df.loc[rsrd_df['TYP'] == t]) == 1)):
            rsrd_df = rsrd_df.loc[rsrd_df['TYP'] != t]
    ob_cts = {'TYP':[], 'COUNT':[]}
    for i in range(len(rsrd_df)):
        ob_cts['TYP'].append(f"{rsrd_df['TYP'].values[i]} (RSRD={rsrd_df['RSRD'].values[i]})")
        ob_cts['COUNT'].append(rsrd_df['COUNT'].values[i])

    # Count number of obs of each type, excluding those ob types that include restricted obs
    all_typ = np.unique(bufr_csv.df['TYP'])
    for t in all_typ:
        if (t in typ) and (t not in rsrd_df['TYP'].values):
            ob_cts['TYP'].append(f"{int(t)}")
            ob_cts['COUNT'].append(count_ob_typ(bufr_csv, typ=[t]))

    return pd.DataFrame(ob_cts)


def compute_ob_pct(ob_cts):
    """
    Compute the percent of obs in each category
    """

    tot = np.sum(ob_cts['COUNT'].values)
    ob_cts['PCT'] = 100 * ob_cts['COUNT'] / tot

    return ob_cts


def aggregate_ob_cts(ob_cts):
    """
    Combine counts for similar ob types
    """

    ob_groups = {'raob':['120', '122', '132', '220', '221', '222', '232'],
                 'profiler':['126', '227', '229'],
                 'satellite':['151', '153', '156', '158', '164', '174', 
                              '242', '243', '245', '246', '247', '250', '251', '252', '253', '254', '257', '258', '259', '289', '290'],
                 'other sfc':['180', '181', '182', '183', '187', '192', '193', '194', 
                              '280', '281', '282', '284', '287', '292', '293', '294'],
                 'VAD':['224']}

    ob_groups = {'other':['120', '122', '132', '220', '221', '222', '232',
                          '130', '230',
                          '126', '227', '229',
                          '151', '153', '156', '158', '164', '174', 
                          '242', '243', '245', '246', '247', '250', '251', '252', '253', '254', '257', '258', '259', '289', '290',
                          '224'],
                 'restricted aircraft':['131', '133', '134', '135', '231', '233', '234', '235'],
                 'other sfc':['180', '181', '182', '183', '187', '192', '193', '194', 
                              '280', '281', '282', '284', '287', '292', '293', '294']}

    # Count number of obs in each group
    ob_cts_new = copy.deepcopy(ob_cts)
    for g in ob_groups.keys():
        n = 0
        for t in ob_groups[g]:
            if t in ob_cts['TYP'].values:
                n = n + ob_cts.loc[ob_cts['TYP'] == t, 'COUNT'].values
                ob_cts_new = ob_cts_new.loc[ob_cts_new['TYP'].values != t, :]
        if n > 0:
            ob_cts_new = pd.concat([ob_cts_new, pd.DataFrame({'TYP':[g], 'COUNT':[n]})])

    return ob_cts_new


def make_pie_chart(ob_cts, ax=None, pie_kwargs={'autopct':'%1.1f%%'}):
    """
    Make pie chart of observation counts
    """

    if ax == None:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))

    ax.pie(ob_cts['COUNT'], labels=ob_cts['TYP'], **pie_kwargs)

    return ax


if __name__ == "__main__":
    fname = sys.argv[1]
    bufr_csv = bufr.bufrCSV(fname)

    # Thermodynamic obs
    print('counting thermodynamic obs...')
    bufr_copy = copy.deepcopy(bufr_csv)
    bufr_copy.df = bufr_copy.df.loc[~np.isnan(bufr_copy.df['TOB']), :]
    print(f'Length of BUFR CSV copy = {len(bufr_copy.df)}')
    ob_cts_thermo = count_all_obs(bufr_copy, typ=list(range(100, 190)))
    ob_cts_thermo = compute_ob_pct(ob_cts_thermo)
    print(ob_cts_thermo)
    ob_cts_thermo_agg = aggregate_ob_cts(ob_cts_thermo)
    ax = make_pie_chart(ob_cts_thermo_agg)
    ax.set_title(f'Thermodynamic Obs (n = {len(bufr_copy.df)})', size=14)
    plt.suptitle(bufr_csv.df['cycletime'].values[0], size=18)
    plt.savefig('thermo_obs.png')
    plt.show()
    ob_cts_thermo.to_csv('thermo_obs.csv')

    # Kinematic obs
    print()
    print('counting kinematic obs...') 
    bufr_copy = copy.deepcopy(bufr_csv)
    bufr_copy.df = bufr_copy.df.loc[~np.isnan(bufr_copy.df['UOB']), :]
    print(f'Length of BUFR CSV copy = {len(bufr_copy.df)}')
    ob_cts_kine = count_all_obs(bufr_copy, typ=list(range(200, 290)))
    ob_cts_kine = compute_ob_pct(ob_cts_kine)
    print(ob_cts_kine)
    ob_cts_kine_agg = aggregate_ob_cts(ob_cts_kine)
    ax = make_pie_chart(ob_cts_kine_agg)
    ax.set_title(f'Kinematic Obs (n = {len(bufr_copy.df)})', size=14)
    plt.suptitle(bufr_csv.df['cycletime'].values[0], size=18)
    plt.savefig('kine_obs.png')
    plt.show()
    ob_cts_kine.to_csv('kine_obs.csv')


"""
End count_restricted.py 
"""
