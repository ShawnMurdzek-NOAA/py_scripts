"""
Compute Differences Between Two Model Forecasts, Then Convolve Results

Difference fields are convolved with a 2D array of ones with a desired size. Convolution output is
converted to RMSDs. The result is an array of RMSDs for different "patches" of the difference field.
This allows users to isolate small regions of the two forecasts where there are large differences.

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import numpy as np
import scipy.signal as ss
import pandas as pd
import yaml
import sys
import argparse
import datetime as dt
import copy

import pyDA_utils.upp_postprocess as uppp


#---------------------------------------------------------------------------------------------------
# Main Body
#---------------------------------------------------------------------------------------------------

def parse_args(argv):
    """
    Parse command-line arguments

    Parameters
    ----------
    argv : list
        Command-line arguments

    Returns
    -------
    args : dictionary
        Dictionary with parsed command-line arguments

    """

    parser = argparse.ArgumentParser(description='Returns the result of a 2D convolution of the difference between two forecasts.')
    parser.add_argument('YAML_file', help='YAML file with program options')

    args = parser.parse_args(argv)

    return args


def parse_in_yaml(args):
    """
    Parse input YAML file

    Parameters
    ----------
    args : dictionary
        Dictionary with parsed command-line arguments

    Returns
    -------
    param : dictionary
        Program parameters
        
    """

    # Read YAML
    with open(args.YAML_file, 'r') as fptr:
        param = yaml.safe_load(fptr)

    # Expand init time list
    t_start = dt.datetime.strptime(param['model']['init']['start'], '%Y%m%d%H')
    t_stop = dt.datetime.strptime(param['model']['init']['stop'], '%Y%m%d%H')
    param['model']['init'] = pd.date_range(start=t_start,
                                           end=t_stop,
                                           freq=param['model']['init']['step']).tolist()

    # Raise error if rmsd_fhr_diff = True, but fcst_hr does not have 2 entries
    if param['model']['rmsd_fhr_diff']:
        if len(param['model']['fcst_hr']) != 2:
            raise ValueError('model/fcst_hr must have 2 entries if model/rmsd_fhr_diff = True')

    return param


def read_fcst_output(param, itime, fhr):
    """
    Read forecast output files

    """

    fnames = {}
    for key in ['file1', 'file2']:
        dum = copy.deepcopy(param['model'][key])
        dum = dum.replace('{iYYYY}', itime.strftime('%Y'))
        dum = dum.replace('{iMM}', itime.strftime('%m'))
        dum = dum.replace('{iDD}', itime.strftime('%d'))
        dum = dum.replace('{iHH}', itime.strftime('%H'))
        dum = dum.replace('{vYYYY}', (itime + dt.timedelta(hours=fhr)).strftime('%Y'))
        dum = dum.replace('{vMM}', (itime + dt.timedelta(hours=fhr)).strftime('%m'))
        dum = dum.replace('{vDD}', (itime + dt.timedelta(hours=fhr)).strftime('%d'))
        dum = dum.replace('{vHH}', (itime + dt.timedelta(hours=fhr)).strftime('%H'))
        dum = dum.replace('{FFF}', f'{fhr:03d}')
        fnames[key] = dum

    ds1 = xr.open_dataset(fnames['file1'], engine='pynio')
    ds2 = xr.open_dataset(fnames['file2'], engine='pynio')

    return ds1, ds2


def extract_field(ds, param):

    name = param['field']['name']

    # Compute derived fields
    if name in ['WSPD', 'WDIR']:
        ds = uppp.compute_wspd_wdir(ds)
    if name == 'CEIL':
        ds = uppp.compute_ceil_agl(ds, no_ceil=2e4)
        ds['CEIL'] = ds[param['field']['ceil_name']]

    # Extract field
    field = ds[name].values
    if len(field.shape) == 3:
        field = field[param['field']['level'], :, :]

    # Make field binary, if desired
    if param['field']['binary']:
        field = np.int64(field > param['field']['binary_val'])

    return field


def extract_lat_lon(ds, param):

    c_size = param['convolve']['size']
    subset = param['convolve']['subset']

    idx = int(c_size / 2)
    if c_size % 2 == 0:
        lat = ds['gridlat_0'][idx:-idx+1:subset, idx:-idx+1:subset].values
        lon = ds['gridlon_0'][idx:-idx+1:subset, idx:-idx+1:subset].values
    else:
        lat = ds['gridlat_0'][idx:-idx:subset, idx:-idx:subset].values
        lon = ds['gridlon_0'][idx:-idx:subset, idx:-idx:subset].values

    return lat, lon


def perform_diff_2d_convolve(field1, field2, param):

    diff = (field1 - field2)**2
    kernel = np.ones([param['convolve']['size']]*2)
    rmsd = np.sqrt(ss.convolve2d(diff, kernel, mode='valid') / np.sum(kernel))

    # Subset RMSDs
    subset = param['convolve']['subset']
    rmsd = rmsd[::subset, ::subset]

    return rmsd


def rmsd_diff_field_names(param):

    fhr1 = param['model']['fcst_hr'][0]
    fhr2 = param['model']['fcst_hr'][1]

    return f"rmsd_{fhr1}", f"rmsd_{fhr2}"


def init_out_dict(param):

    out = {'lat':[],
           'lon':[],
           'init':[]}

    if param['model']['rmsd_fhr_diff']:
        n1, n2 = rmsd_diff_field_names(param)
        out[n1] = []
        out[n2] = []
    else:
        out['rmsd'] = []
        out['fhr'] = []

    return out


def save_output_1iter(param, rmsd, lat, lon, init, fhr, out):

    size = lat.size
    out['lat'] = out['lat'] + list(lat.flatten())
    out['lon'] = out['lon'] + list(lon.flatten())
    out['init'] = out['init'] + [init]*size

    if param['model']['rmsd_fhr_diff']:
        n1, n2 = rmsd_diff_field_names(param)
        out[n1] = out[n1] + list(rmsd[0].flatten())
        out[n2] = out[n2] + list(rmsd[1].flatten())
    else:
        out['rmsd'] = out['rmsd'] + list(rmsd.flatten())
        out['fhr'] = out['fhr'] + [fhr]*size

    return out


def save_output_fname(out, param):

    # Create DataFrame
    out_df = pd.DataFrame.from_dict(out)

    # Determine output file name
    tag = param['out']['tag']
    var = param['field']['name'].split('_')[0]
    lvl = param['field']['level']
    if param['field']['binary']:
        fname = f"{var}_{lvl}_binary{param['field']['binary_val']}{tag}.csv"
    else:
        fname = f"{var}_{lvl}_{tag}.csv"

    # Sort by RMSD
    if param['model']['rmsd_fhr_diff']:
        n1, n2 = rmsd_diff_field_names(param)
        out_df['rmsd_diff'] = out_df[n2] - out_df[n1]
        out_df.sort_values('rmsd_diff', ascending=False, inplace=True)
    else:
        out_df.sort_values('rmsd', ascending=False, inplace=True)

    # Save to CSV
    out_df.to_csv(fname)

    return out_df


def compute_rmsds_1init(init, param, out_dict, verbose=1):

    if param['model']['rmsd_fhr_diff']:
        rmsd = []

    for fhr in param['model']['fcst_hr']:
        if verbose > 0: print(f"init = {t}, fhr = {fhr}")
        ds1, ds2 = read_fcst_output(param, t, fhr)
        field1 = extract_field(ds1, param)
        field2 = extract_field(ds2, param)
        lat, lon = extract_lat_lon(ds1, param)
        if param['model']['rmsd_fhr_diff']:
            rmsd.append(perform_diff_2d_convolve(field1, field2, param))
        else:
            rmsd = perform_diff_2d_convolve(field1, field2, param)
            out_dict = save_output_1iter(param, rmsd, lat, lon, t, fhr, out_dict)

    if param['model']['rmsd_fhr_diff']:
        out_dict = save_output_1iter(param, rmsd, lat, lon, t, fhr, out_dict)

    return out_dict


if __name__ == '__main__':

    # Record start time
    start = dt.datetime.now()
    print()
    print('Starting diff_2d_convolve')
    print(f"Start time = {start.strftime('%Y%m%d %H:%M:%S')}")

    # Parse command-line args
    args = parse_args(sys.argv[1:])
    param = parse_in_yaml(args)
    print(f"Compute RMSD Diffs? {param['model']['rmsd_fhr_diff']}")
    print()

    # Dictionary to save output
    out_dict = init_out_dict(param)

    # Main loop
    for t in param['model']['init']:
        out_dict = compute_rmsds_1init(t, param, out_dict, verbose=1)

    # Save output to CSV file
    out_df = save_output_fname(out_dict, param)

    # Finish
    print()
    print(f"Done. Elapsed time = {(dt.datetime.now() - start).total_seconds()} s")


"""
End diff_2d_convolve.py
"""
