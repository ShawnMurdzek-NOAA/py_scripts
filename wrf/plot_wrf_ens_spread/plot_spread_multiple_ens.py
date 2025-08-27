"""
Create Line Plots of Ensemble Spread from Several Ensembles

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import datetime as dt
import sys
import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import plot_wrf_ens_spread as pwes


#---------------------------------------------------------------------------------------------------
# Main Program
#---------------------------------------------------------------------------------------------------

def modify_paths_multi(config):
    """
    Combine parent and ensemble paths in paths_multi

    Parameters
    ----------
    config : dictionary
        Program configuration (from YAML)

    Returns
    -------
    Updated config where parent and ensemble paths are combined

    """
    
    # Add parent directory to ensemble member directories
    for ens in config['paths_multi'].keys():
        config['paths_multi'][ens]['_mem'] = []
        if 'parent' in config['paths_multi'][ens]:
            root = config['paths_multi'][ens]['parent']
            for i in range(len(config['paths_multi'][ens]['ens_mem'])):
                config['paths_multi'][ens]['_mem'].append(f"{root}/{config['paths_multi'][ens]['ens_mem'][i]}")
        else:
            config['paths_multi'][ens]['_mem'] = config['paths_multi'][ens]['ens_mem']
    
    return config


def compute_wrf_field_stat_multi_ens(config):
    """
    Compute WRF model statistics for several ensembles

    Parameters
    ----------
    config : dictionary
        Program configuration (from YAML)

    Returns
    -------
    Dictionary of fields to plot and metadata

    """
    
    # Loop over each ensemble and compute stats
    out_dict = {}
    for ens in config['paths_multi'].keys():
        print(f'\nReading ensemble {ens}')
        config_1_ens = copy.deepcopy(config)
        config_1_ens['paths'] = config['paths_multi'][ens]
        out_dict[ens], times, metadata = pwes.main_time_loop(config_1_ens)
    
    return out_dict, times, metadata


def avg_dist_from_ctr(x, y):
    """
    Compute the average distance from the centroid of a collection of points

    Parameters
    ----------
    x : array
        X coordinates
    y : array
        Y coordinates

    Returns
    -------
    scalar

    """
    
    # Compute center
    xc = np.mean(x)
    yc = np.mean(y)
    
    # Compute distances from the center
    dist = np.sqrt((x - xc)**2 + (y - yc)**2)
    
    # Compute average distance
    return np.mean(dist)


def plot_ens_stat_timeseries(stat_dict, times, opt, meta, plot_dict):
    """
    Make timeseries plots of ensemble statistics

    Parameters
    ----------
    stat_dict : dictionary
        Statistics from each WRF ensemble member from each ensemble
    times : list of dt.datetimes
        Times corresponding to the statistics in stat_dict
    opt : dictionary
        Specific plotting options for each plot
    meta : dictionary
        Metadata for variables being plotted
    plot_dict : dictionary
        General plotting options

    Returns
    -------
    None.

    """
    
    # Determine spacing between x axis ticks
    delta_min = (times[1] - times[0]).total_seconds() / 60
    int_min = int(np.ceil(len(times) / 10) * delta_min)
    ntime = len(times)
    
    for f in mem_stat_dict[list(mem_stat_dict.keys())[0]]:
        spread = {}
        if opt[f]['ptype'] == 'timeseries':
            stat = 'std dev'
            for ens in stat_dict.keys():
                spread[ens] = np.std(stat_dict[ens][f], axis=1)
        elif opt[f]['ptype'] == 'spatial_2d':
            stat = 'avg dist'
            for ens in stat_dict.keys():
                spread[ens] = np.zeros(ntime)
                for k in range(ntime):
                    spread[ens][k] = avg_dist_from_ctr(stat_dict[ens][f][k, :, 0],
                                                       stat_dict[ens][f][k, :, 1])
        print(f"Making timeseries plot for {f}")
            
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
            
        colors = ['k', '#66CCEE', '#AA3377', '#228833', '#CCBB44', '#4477AA', 'gray']
        for i, ens in enumerate(spread.keys()):
            ax.plot(times, spread[ens], c=colors[i], ls='-', lw=1,
                    label=f"{ens} (avg {stat} = {np.mean(spread[ens]):.3e})")
            
        ax.grid(lw=0.25)
        ax.legend()
        ax.set_title(f, size=16)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
        ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=int_min))
        ax.tick_params(axis='x', rotation=45)
        ax.set_ylabel(f"{stat} {opt[f]['stat']}\n{meta[f]['description']} ({meta[f]['units']})", 
                      size=12)
            
        plt.savefig(f"{f}_ens_stat_timeseries.png")
        plt.close()


if __name__ == '__main__':

    start = dt.datetime.now()
    print('Starting plot_spread_multiple_ens.py')
    print(f"Time = {start.strftime('%Y%m%d %H:%M:%S')}\n")

    # Read inputs
    param = pwes.parse_in_args(sys.argv[1:])
    config = pwes.read_yaml(param.in_yaml)
    config = modify_paths_multi(config)

    # Compute WRF model stats for each ensemble
    mem_stat_dict, times, meta = compute_wrf_field_stat_multi_ens(config)
    print()

    # Create plots
    plot_ens_stat_timeseries(mem_stat_dict, times, config['opt'], meta, config['plot'])

    print('\nProgram finished!')
    print(f"Elapsed time = {(dt.datetime.now() - start).total_seconds()} s")


"""
End plot_spread_multiple_ens.py
"""
