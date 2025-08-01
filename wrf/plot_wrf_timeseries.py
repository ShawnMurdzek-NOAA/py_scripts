"""
Plot Timeseries of Output from WRF History Dumps

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import datetime as dt
import argparse
import sys
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import glob


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

    parser = argparse.ArgumentParser(description='Script that creates timeseries plots from wrfout \
                                                  files')
    
    # Positional arguments
    parser.add_argument('file_list', 
                        help='CSV file containing the name of each simulation and the \
                              corresponding to the path containing the wrfout netCDF files',
                        type=str)
    
    parser.add_argument('field', 
                        help='Field used to to compute the timeseries',
                        type=str)

    parser.add_argument('reduction', 
                        help="Method for reducing the desired field to a single value for each \
                              time. Options: 'min', 'max', 'mean', 'median'",
                        type=str)

    parser.add_argument('outfile', 
                        help='Output file',
                        type=str)

    return parser.parse_args(argv)


def parse_sim_paths(file_list):
    """
    Parse file_list CSV into a dictionary for easier use
    """

    file_df = pd.read_csv(file_list, names=['name', 'path'])
    file_dict = {}
    for i in range(len(file_df)):
        file_dict[file_df['name'].values[i]] = file_df['path'].values[i]

    return file_dict


def create_timeseries(file_dict, field, red):
    """
    Compute timeseries using the desired reduction method
    """

    timeseries = {}
    for key in file_dict.keys():
        print(f"Opening wrfout files for {key}")
        files = glob.glob(f"{file_dict[key]}/wrfout*")
        files.sort()
        timeseries[key] = {'time': np.empty(len(files), dtype=object),
                           'values': np.zeros(len(files))}
        for j, f in enumerate(files):
            ds = xr.open_dataset(f)
            timeseries[key]['time'][j] = dt.datetime.strptime(ds['Times'].values[0].decode('utf-8'),
                                                              '%Y-%m-%d_%H:%M:%S')
            timeseries[key]['values'][j] = reduce_field(ds, field, red)

    # Save metadata
    metadata = ds[field].attrs
    metadata['reduction'] = red

    return timeseries, metadata


def reduce_field(ds, field, red):
    """
    Reduce a multidimensional field to a single value
    """

    array = ds[field].values
    if red == 'max':
        value = np.amax(array)
    elif red == 'min':
        value = np.amin(array)
    elif red == 'mean':
        value = np.mean(array)
    elif red == 'median':
        value = np.median(array)
    else:
        return ValueError(f"Variable reduction {red} is not supported")

    return value


def plot_timeseries(timeseries, metadata, outfile):
    """
    Plot time series and save to output file
    """

    # Define colors
    colors = ['k', '#66CCEE', '#AA3377', '#228833', '#CCBB44', '#4477AA']

    # Configure figure
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.97, top=0.97)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=15))

    # Plot data
    for key, c in zip(timeseries.keys(), colors):
        ax.plot(timeseries[key]['time'], timeseries[key]['values'], c=c, ls='-', label=key)

    # Make plot look nice
    ax.set_ylabel(f"{metadata['reduction']} {metadata['description']} {metadata['units']}", size=14)
    ax.legend(fontsize=14)
    ax.tick_params(axis='x', rotation=45)
    ax.tick_params(axis='both', labelsize=12)
    ax.grid()

    plt.savefig(outfile)

    return None


if __name__ == '__main__':

    start = dt.datetime.now()
    print('Starting plot_wrf_timeseries.py')
    print(f"Time = {start.strftime('%Y%m%d %H:%M:%S')}")

    # Read in data
    param = parse_in_args(sys.argv[1:])
    file_dict = parse_sim_paths(param.file_list)

    # Compute timeseries
    timeseries, metadata = create_timeseries(file_dict, param.field, param.reduction)

    # Plot timeseries
    plot_timeseries(timeseries, metadata, param.outfile)

    print('Program finished!')
    print(f"Elapsed time = {(dt.datetime.now() - start).total_seconds()} s")


"""
End plot_wrf_timeseries.py
"""
