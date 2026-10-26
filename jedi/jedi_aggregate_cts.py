"""
Aggregate ob counts from several pyDAmonitor obs_count.txt files

Shawn Murdzek
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import sys
import argparse
import datetime as dt
import pandas as pd


#---------------------------------------------------------------------------------------------------
# Main Program
#---------------------------------------------------------------------------------------------------

def parse_in_args(argv):
    """
    Parse input arguments
    """

    parser = argparse.ArgumentParser(description='Script that parses obs_count.txt files from \
                                                  pyDAmonitor to get aggregate ob counts')
    
    # Positional arguments
    parser.add_argument('rrfs_exp', 
                        help='Path to rrfs.YYYYMMDD directories in RRFS run',
                        type=str)

    # Optional arguments
    parser.add_argument('--start',
                        dest='start_time',
                        default='2024050613',
                        help='Start time (YYYYMMDDHH)',
                        type=str)

    parser.add_argument('--end',
                        dest='end_time',
                        default='2024051223',
                        help='End time (YYYYMMDDHH)',
                        type=str)

    parser.add_argument('--exclude_00_12',
                        dest='exclude_00_12',
                        default=1,
                        help='Option to exclude 00 and 12z initializations',
                        type=int)

    return parser.parse_args(argv)


def total_ob_counts(path, start, end, skip=[0, 12]):
    """
    Create a main ob count DataFrame that concatenates output from all obs_count.txt files
    """

    # Iterate over all output times
    df_list = []
    current = start - dt.timedelta(hours=1)
    while current <= end:

        # Increment by 1 hr
        current = current + dt.timedelta(hours=1)

        # Skip if hours in skip list
        if current.hour in skip:
            continue

        # Open DataFrame and save to list
        fname = f"{path}/rrfs.{current.strftime('%Y%m%d')}/{current.strftime('%H')}/pyDAmonitor/det/web/obs_count.txt"
        try:
            tmp_df = pd.read_csv(fname, sep='\s+')
            tmp_df['timestamp'] = current.strftime('%Y%m%d%H')
            df_list.append(tmp_df)
        except FileNotFoundError:
            print(f"Missing file: {fname}")
            continue       

    return pd.concat(df_list, ignore_index=True)


if __name__ == '__main__':

    # Read input parameters  
    in_param = parse_in_args(sys.argv[1:])
    start_dt = dt.datetime.strptime(in_param.start_time, '%Y%m%d%H')
    end_dt = dt.datetime.strptime(in_param.end_time, '%Y%m%d%H')

    # Create DataFrame with total ob counts
    if in_param.exclude_00_12 == 1:
        tot_ob_cts = total_ob_counts(in_param.rrfs_exp, start_dt, end_dt, skip=[0, 12])
    elif in_param.exclude_00_12 == 0:
        tot_ob_cts = total_ob_counts(in_param.rrfs_exp, start_dt, end_dt, skip=[])
    else:
        raise ValueError(f"Invalid --exclude_00_12 option : {in_param.exclude_00_12}")

    # Print results
    print()
    print(tot_ob_cts.groupby('observer')[['n_ioda', 'nobs', 'nobs_r', 'n_loop1', 'n_loop2']].sum())


"""
End jedi_aggregate_cts.py
"""
