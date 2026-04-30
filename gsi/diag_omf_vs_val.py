"""
Plot OmF vs Obs Value

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import datetime as dt
import sys
import argparse

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

    parser = argparse.ArgumentParser(description='Plot OmF vs obs value using a 2D histogram')
    
    # Positional arguments
    parser.add_argument('diag_tmpl', 
                        help='GSI diag file template. Include strftime formatters for date.',
                        type=str)
    
    parser.add_argument('diag_type', 
                        help="GSI diag file type ('text' or 'netcdf')",
                        type=str)

    parser.add_argument('start_time', 
                        help='Start time in YYYYMMDDHH format',
                        type=int)

    parser.add_argument('end_time', 
                        help='End time in YYYYMMDDHH format',
                        type=str)

    parser.add_argument('tag', 
                        help='Output file tag',
                        type=str)

    # Optional arguments
    parser.add_argument('--ob_type',
                        dest='ob_type',
                        default='all',
                        help="Observation types to use. Set to 'all' to retain all obs types",
                        type=str)

    parser.add_argument('--plot_log',
                        dest='plot_log',
                        default=0,
                        help="Option to make plot logarithmic",
                        type=int)

    parser.add_argument('--pseudo_rh',
                        dest='pseudo_rh',
                        default=0,
                        help="Option toconvert observations to pseudo relative humidities by \
                              dividing by the forecast saturation specific humidity",
                        type=int)

    return parser.parse_args(argv)


def read_diag(diag_tmpl, diag_type, start, end, ob_type='all'):
    """
    Read GSI diag file
    """

    # Create list of datetime objects
    start_dt = dt.datetime.strptime(str(start), '%Y%m%d%H')
    end_dt = dt.datetime.strptime(str(end), '%Y%m%d%H')
    dates = [start_dt]
    while dates[-1] < end_dt:
        dates.append(dates[-1] + dt.timedelta(hours=1))

    # Create list of file names
    fnames = [d.strftime(diag_tmpl) for d in dates]

    # Read in diag file and only keep desired ob types
    df = gsi.read_diag(fnames, ftype=diag_type)
    if ob_type != 'all':
        df = df.loc[df['Observation_Type'] == int(ob_type), :]

    return df


def plot_omf_vs_ob(df, ttl='', log=0, pseudo_rh=0):
    """
    Plot OmF vs observation value
    """

    # Compute pseudo RH
    if pseudo_rh == 1:
        df['pseudo_rh'] = df['Observation'] / df['Forecast_Saturation_Spec_Hum']
        xname = 'pseudo_rh'
    else:
        xname = 'Observation'

    # Determine min and max plotting value
    yname = 'Obs_Minus_Forecast_adjusted'
    max_xval = np.percentile(df[xname], 99.5)
    min_xval = np.percentile(df[xname], 0.5)
    max_yval = np.percentile(df[yname], 99.5)
    min_yval = np.percentile(df[yname], 0.5)

    # Make plot
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 6), sharex=True, sharey=True)
    plt.subplots_adjust(wspace=0.05, left=0.08, bottom=0.06, right=0.98, top=0.87)
    for i, subset in enumerate(['all', 'assim']):
        ax = axes[i]
        if subset == 'all':
            xdata = df[xname].values
            ydata = df[yname].values
        elif subset == 'assim':
            cond = df['Analysis_Use_Flag'] == 1
            xdata = df[xname].loc[cond].values
            ydata = df[yname].loc[cond].values
        hist = np.histogram2d(xdata, ydata, bins=40, range=[[min_xval, max_xval], [min_yval, max_yval]], density=True)
        hist[0][hist[0] == 0] = np.nan
        if subset == 'all': vmax = np.nanmax(hist[0])
        if log == 1:
            cax = ax.pcolormesh(hist[1], hist[2], hist[0].T, shading='flat', cmap='plasma', 
                                norm=colors.LogNorm(vmin=1, vmax=vmax))
        else:
            cax = ax.pcolormesh(hist[1], hist[2], hist[0].T, shading='flat', cmap='plasma', vmin=1, vmax=vmax)

        ax.ticklabel_format(style="sci", axis="both", scilimits=(0, 0))
        ax.set_title(f"{subset} (n = {len(xdata)})", size=14)
        ax.grid()
        ax.set_xlabel(xname, size=12)
        if i == 0: ax.set_ylabel(yname, size=12)

    # Add colorbar and other annotations
    cbar = plt.colorbar(cax, ax=axes, orientation='horizontal', aspect=35)
    cbar.set_label('density', size=12)
    var = df['Observation_Class'].values[0].decode('utf-8').strip()
    plt.suptitle(f"{ttl} {var}", size=16)

    return fig


if __name__ == '__main__':

    start = dt.datetime.now()
    print('Starting diag_omf_vs_val.py')
    print(f"Time = {start.strftime('%Y%m%d %H:%M:%S')}")

    # Read in GSI diag info
    param = parse_in_args(sys.argv[1:])
    omf_df = read_diag(param.diag_tmpl, param.diag_type, param.start_time, param.end_time, ob_type=param.ob_type)

    # Make plot
    ttl = f"{param.start_time}$-${param.end_time} type={param.ob_type}"
    fig = plot_omf_vs_ob(omf_df, ttl=ttl, log=param.plot_log, pseudo_rh=param.pseudo_rh)
    plt.savefig(f"omf_vs_ob_{param.start_time}_{param.end_time}_{param.ob_type}_{param.tag}.png")

    print('Program finished!')
    print(f"Elapsed time = {(dt.datetime.now() - start).total_seconds()} s")


"""
End diag_omf_vs_val.py
"""
