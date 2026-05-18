"""
Plot Percent Diffs in Added Errors Between the 1st and 2nd Iterations

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import sys
import argparse
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np

import pyDA_utils.gsi_fcts as gsi


#---------------------------------------------------------------------------------------------------
# Main Program
#---------------------------------------------------------------------------------------------------

def parse_in_args(argv):
    """
    Parse input arguments

    Parameters
    ----------
    argv : list
        Command-line arguments from sys.argv[1:]

    Returns
    -------
    Parsed input arguments

    """

    parser = argparse.ArgumentParser(description='Script that plots percent differences between \
                                                  two GSI errtable files')

    # Positional arguments
    parser.add_argument('file1',
                        help='First GSI errtable file',
                        type=str)

    parser.add_argument('file2',
                        help='Second GSI errtable file',
                        type=str)

    # Optional arguments
    parser.add_argument('--tag',
                        dest='tag',
                        default='',
                        help='Tag to add to output file name',
                        type=str)

    return parser.parse_args(argv)


def plot_errtable_pct_diff(errtable1, errtable2, tag=''):
    """
    Create a bar chart showing percent differences between two errtables
    """

    all_colors = ['k', '#66CCEE', '#AA3377', '#228833', '#CCBB44', '#4477AA']

    # Compute percent diffs
    labels = []
    pct_diffs = []
    colors = []
    for typ in errtable1.keys():
        for v, c in zip(errtable1[typ].columns, all_colors):
            mean_e1 = np.nanmean(errtable1[typ][v])
            mean_e2 = np.nanmean(errtable2[typ][v])
            if v != 'prs' and (mean_e1 > 0): 
                labels.append(f"{v}_{typ}")
                pct_diffs.append(100*(mean_e2 - mean_e1) / mean_e1)
                colors.append(c)

    # Sort arrays for cosmetic reasons
    labels = np.array(labels)
    sort_idx = np.argsort(labels)
    labels = labels[sort_idx]
    pct_diffs = np.array(pct_diffs)[sort_idx]
    colors = np.array(colors)[sort_idx]

    # Create plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    ax.barh(labels, pct_diffs, color=colors)
    ax.axvline(0, color='k', linewidth=0.75)
    ax.invert_yaxis()
    ax.set_xlabel('GSI errtable % differences\n(errors are averaged vertically prior to differencing)', size=12)
    ax.grid()
    ax.set_title(tag, size=16)

    # Save plot
    if len(tag) > 0:
        tag = f"_{tag}"
    plt.savefig(f"errtable_pct_diff{tag}.png", dpi=400, bbox_inches="tight")

    return None


if __name__ == '__main__':

    start = dt.datetime.now()
    print('Starting plot_gsi_errtable_pct_diff.py')
    print(f"Time = {start.strftime('%Y%m%d %H:%M:%S')}")

    # Read in parameters and open GSI errtable files
    param = parse_in_args(sys.argv[1:])
    err1 = gsi.read_errtable(param.file1)
    err2 = gsi.read_errtable(param.file2)

    # Make plot
    plot_errtable_pct_diff(err1, err2, tag=param.tag)

    print('Program finished!')
    print(f"Elapsed time = {(dt.datetime.now() - start).total_seconds()} s")


"""
End plot_gsi_errtable_pct_diff.py
"""
