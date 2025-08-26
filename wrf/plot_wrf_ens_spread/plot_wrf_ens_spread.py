"""
Plot WRF Ensemble Spread

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import datetime as dt
import yaml
import sys
import argparse
import copy
import numpy as np
import netCDF4 as nc
import wrf
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


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

    parser = argparse.ArgumentParser(description='Script that plots various measures of ensemble \
                                                  spread for and idealized WRF ensemble of a \
                                                  supercell.')
    
    # Positional arguments
    parser.add_argument('in_yaml', 
                        help='Input YAML file. Will overwrite contents of default.yml',
                        type=str)
    
    return parser.parse_args(argv)


def recursive_dict_update(defaults, updates):
    """
    Recursively update a default dictionary with updates from another dictionary

    Parameters
    ----------
    defaults : dictionary
        Default dictionary to be updated
    updates : dictionary
        Values used to update the default dictionary

    Returns
    -------
    Updated dictionary

    """
    
    out_dict = copy.deepcopy(defaults)
    for k, v in updates.items():
        if isinstance(v, dict) and (k in out_dict) and isinstance(out_dict[k], dict):
            out_dict[k] = recursive_dict_update(out_dict[k], v)
        else:
            out_dict[k] = v
            
    return out_dict


def read_yaml(in_yaml):
    """
    Read input YAML file, merge with default.yml, and reformat entries

    Parameters
    ----------
    in_yaml : string
        Input YAML file

    Returns
    -------
    Dictionary with YAML settings

    """
    
    # Read in both YAML files
    with open('default.yml', 'r') as fptr:
        defaults = yaml.safe_load(fptr)
        
    with open(in_yaml, 'r') as fptr:
        updates = yaml.safe_load(fptr)
    
    # Combine dictionaries
    out_dict = recursive_dict_update(defaults, updates)
    
    # Check out_dict for necessary parameters
    check_yaml_contents(out_dict)
    
    # Add parent directory to ensemble member directories
    if 'parent' in out_dict['paths']:
        root = out_dict['paths']['parent']
        for i in range(len(out_dict['paths']['ens_mem'])):
            out_dict['paths']['ens_mem'][i] = f"{root}/{out_dict['paths']['ens_mem'][i]}"
    
    # Change times to datetimes and time differences to timedeltas
    for key in ['start', 'stop']:
        out_dict['time'][key] = dt.datetime.strptime(out_dict['time'][key],
                                                     '%Y-%m-%d-%H-%M-%S')
    out_dict['time']['step'] = dt.timedelta(seconds=out_dict['time']['step'])
    
    return out_dict


def check_yaml_contents(yaml_dict):
    """
    Check that all necessary fields are in the input YAML file

    Parameters
    ----------
    yaml_dict : dictionary
        YAML input parameters

    Returns
    -------
    None

    """
    
    # Check out_dict for necessary parameters
    necessary_param = {'paths': ['ens_mem'],
                       'time': ['start', 'stop', 'step'],
                       'fields': []}
    for key in necessary_param.keys():
        if key not in yaml_dict:
            raise ValueError(f"Input YAML is missing a '{key}' section")
        for key2 in necessary_param[key]:
            if key2 not in yaml_dict[key]:
                raise ValueError(f"Input YAML is missing '{key2}' in the '{key}' section")
    
    # Check that each field has a corresponding opt section
    opt_keys = ['var', 'stat', 'ptype']
    for f in yaml_dict['fields']:
        if f in yaml_dict['opt']:
            for key in opt_keys:
                if key not in yaml_dict['opt'][f]:
                    raise ValueError(f"Input YAML is missing '{key}' for {f}")
        else:
            ValueError(f"Input YAML is missing an 'opt' section for {f}")
    
    return None


def compute_wrf_field_stat(fname, fields, opt):
    """
    Open WRF file, extract desired variables, and compute the desired statistics

    Parameters
    ----------
    fname : string
        WRF file name
    fields : list
        Fields to extract. Additional details are given in opt
    opt : dictionary
        Options related to extracting the desired WRF field and computing the statistic

    Returns
    -------
    Statistic and associated metadata

    """
    
    # Open WRF file
    fptr = nc.Dataset(fname)
    
    stats = {}
    meta = {}
    for f in fields:
        
        # Get desired field from WRF
        try:
            if 'get_kw' in opt[f]:
                v = wrf.getvar(fptr, opt[f]['var'], **opt[f]['get_kw'])
            else:
                v = wrf.getvar(fptr, opt[f]['var'])
        except ValueError:
            print(f"{opt[f]['var']} is not a valid WRF variable. Returning NaN.")
            stats[f] = np.nan
            continue
        
        # Compute statistic
        if opt[f]['stat'] == 'max':
            stats[f] = np.amax(v)
        elif opt[f]['stat'] == 'min':
            stats[f] = np.amin(v)
        elif opt[f]['stat'] == 'mean':
            stats[f] = np.mean(v)
        elif opt[f]['stat'] == 'sum':
            stats[f] = np.sum(v)
        elif opt[f]['stat'] == 'max_loc_2d':
            idx = np.unravel_index(np.argmax(v.values), v.shape)
        elif opt[f]['stat'] == 'min_loc_2d':
            idx = np.unravel_index(np.argmin(v.values), v.shape)
        else:
            print(f"Statistic option {opt[f]['stat']} is not supported. Returning NaN.")
            stats[f] = np.nan
        
        # Additional processing for 2D locations
        if 'loc_2d' in opt[f]['stat']:
            ny, nx = v.shape
            dx = wrf.extract_global_attrs(fptr, 'DX')['DX']
            dy = wrf.extract_global_attrs(fptr, 'DY')['DY']
            stats[f] = np.array([idx[1]*dx - 0.5*(nx-1)*dx, 
                                 idx[0]*dy - 0.5*(ny-1)*dy])
            v.attrs['x_grid'] = np.arange(nx)*dx - 0.5*(nx-1)*dx
            v.attrs['y_grid'] = np.arange(ny)*dy - 0.5*(ny-1)*dy
        
        # Extract metadata
        meta[f] = v.attrs
    
    fptr.close()
    
    return stats, meta  


def main_time_loop(config):
    """
    Main time loop over all WRF ensemble members

    Parameters
    ----------
    config : dictionary
        Program configuration (from YAML)

    Returns
    -------
    Dictionary of fields to plot, model output times, and metadata

    """
    
    # Create list of datetimes to loop over
    times = []
    current = config['time']['start']
    while current <= config['time']['stop']:
        times.append(current)
        current = current + config['time']['step']
        
    # Initialize output dictionary
    out_dict = {}
    for f in config['fields']:
        if 'loc_2d' in config['opt'][f]['stat']:
            out_dict[f] = np.zeros([len(times), len(config['paths']['ens_mem']), 2])
        else:
            out_dict[f] = np.zeros([len(times), len(config['paths']['ens_mem'])])
    
    # Compute statistics for each WRF field
    for i, t in enumerate(times):
        print(f"Reading WRF output for {t.strftime('%Y%m%d %H:%M:%S')}")
        for j, p in enumerate(config['paths']['ens_mem']):
            stats, meta = compute_wrf_field_stat(t.strftime(p), config['fields'], config['opt'])
            for f in config['fields']:
                out_dict[f][i, j] = stats[f]
    
    return out_dict, times, meta


def plot_timeseries(stat_dict, times, opt, meta, plot_dict):
    """
    Make timeseries plots

    Parameters
    ----------
    stat_dict : dictionary
        Statistics from each WRF ensemble member
    times : list of dt.datetimes
        Times corresponding to the statistics in mem_stat_dict
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
    
    for f in mem_stat_dict.keys():
        if opt[f]['ptype'] == 'timeseries':
            avg_mean = np.mean(np.mean(stat_dict[f], axis=1))
            avg_std = np.mean(np.std(stat_dict[f], axis=1))
            print(f"Making timeseries plot for {f} (avg std = {avg_std:.3e})")
            
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
            
            for j in range(stat_dict[f].shape[1]):
                ax.plot(times, stat_dict[f][:, j], c='gray', ls='-', lw=0.5)
            if plot_dict['ctrl']:
                ax.plot(times, stat_dict[f][:, 0], 'k-', lw=1)
            ax.plot(times, np.mean(stat_dict[f], axis=1), 'b-', lw=2, 
                    label=f"mean ({avg_mean:.3e})")
            ax.plot(times, np.std(stat_dict[f], axis=1), 'r--', lw=2, 
                    label=f"std dev ({avg_std:.3e})")
            
            ax.grid(lw=0.25)
            ax.legend()
            ax.set_title(f, size=16)
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=int_min))
            ax.tick_params(axis='x', rotation=45)
            ax.set_ylabel(f"{opt[f]['stat']} {meta[f]['description']} ({meta[f]['units']})", size=12)
            
            plt.savefig(f"{f}_timeseries.png")
            plt.close()


def plot_spatial_2d(stat_dict, times, opt, meta, plot_dict):
    """
    Make 2D spatial plots

    Parameters
    ----------
    stat_dict : dictionary
        Statistics from each WRF ensemble member
    times : list of dt.datetimes
        Times corresponding to the statistics in mem_stat_dict
    opt : dictionary
        Plotting options
    meta : dictionary
        Metadata for variables being plotted
    plot_dict : dictionary
        General plotting options

    Returns
    -------
    None.

    """
    
    for f in mem_stat_dict.keys():
        if opt[f]['ptype'] == 'spatial_2d':
            xstd = np.mean(np.std(stat_dict[f][:, :, 0], axis=1))
            ystd = np.mean(np.std(stat_dict[f][:, :, 1], axis=1))
            print(f"Making 2D spatial plot for {f} (avg stds = {xstd:.3e}, {ystd:.3e})")
            
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
            
            for j in range(stat_dict[f].shape[1]):
                ax.plot(stat_dict[f][:, j, 0], stat_dict[f][:, j, 1], c='gray', ls='-', lw=0.5)
            if plot_dict['ctrl']:
                ax.plot(stat_dict[f][:, 0, 0], stat_dict[f][:, 0, 1], 'k-', lw=1)
            
            ax.set_xlim([np.amin(meta[f]['x_grid']), np.amax(meta[f]['x_grid'])])
            ax.set_ylim([np.amin(meta[f]['y_grid']), np.amax(meta[f]['y_grid'])])
            ax.axes.set_aspect('equal')
            ax.set_title(f"{f}\nTemporally averaged X std dev = {xstd:.3e}\n" +
                         f"Temporally averaged Y std dev = {ystd:.3e}", size=16)
            
            plt.savefig(f"{f}_spatial_2d.png")
            plt.close()


if __name__ == '__main__':

    start = dt.datetime.now()
    print('Starting plot_wrf_ens_spread.py')
    print(f"Time = {start.strftime('%Y%m%d %H:%M:%S')}\n")

    # Read inputs
    param = parse_in_args(sys.argv[1:])
    config = read_yaml(param.in_yaml)
    
    # Compute stats from each WRF ensemble member
    mem_stat_dict, times, metadata = main_time_loop(config)
    print()
    
    # Create plots
    plot_timeseries(mem_stat_dict, times, config['opt'], metadata, config['plot'])
    plot_spatial_2d(mem_stat_dict, times, config['opt'], metadata, config['plot'])
    
    print('\nProgram finished!')
    print(f"Elapsed time = {(dt.datetime.now() - start).total_seconds()} s")


"""
End plot_wrf_ens_spread.py
"""
