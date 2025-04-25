"""
Plot RMSE Differences Between Two Simulations

Command Line Arguments
----------------------
sys.argv[1] : Input YAML file name

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import sys
import yaml
import datetime as dt
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import copy
import cartopy.crs as ccrs

import pyDA_utils.plot_model_data as pmd


#---------------------------------------------------------------------------------------------------
# Main Program
#---------------------------------------------------------------------------------------------------

def read_inputs(fname):
    """
    Read program input parameters
    """

    with open(fname, 'r') as fptr:
        param = yaml.safe_load(fptr)

    param['init_times'] = [param['start_time'] + dt.timedelta(hours=i) for i in range(param['nhrs'])]

    return param


def open_ds(tmpl, time, offset=0, debug=1):
    """
    Open GRIB file as a DataSet
    """

    t = time + dt.timedelta(hours=offset)
    name = t.strftime(tmpl)
    if debug > 0: print(f"fname = {name}")
    try:
        ds = xr.open_dataset(name, engine='pynio')
    except:
        print(f"Issue opening files for {t.strftime('%Y%m%d%H')}")
        ds = None

    return ds


def get_empty_da(ds, rmse_field, rmse_zlvl):
    """
    Get empty DataArray that can be used for the RMSEs
    """

    # Extract DataArray
    if np.isnan(rmse_zlvl):
        da = copy.deepcopy(ds[rmse_field])
    else:
        da = copy.deepcopy(ds[rmse_field][rmse_zlvl, ::])

    # Set values to 0 and change name
    da.values = np.zeros(da.shape)
    da.attrs['long_name'] = f"RMSE {da.attrs['long_name']}"

    return da


def init_out_ds(tmpl, time, rmse_field, rmse_zlvl):
    """
    Initialize empty Datasets to hold RMSEs
    """

    ds = open_ds(tmpl, time)
    da = get_empty_da(ds, rmse_field, rmse_zlvl)

    out_ds = []
    for i in range(2):
        out_ds.append(xr.Dataset({'RMSE':copy.deepcopy(da)}))

    return out_ds


def compute_diff(ds, ds_truth, rmse_field, rmse_zlvl, rmse_zlvl_truth):
    """
    Compute 2D difference between a forecast simulation and a "truth" simulation (assumed to 
    be the 1-km NR)
    """

    # Compute difference
    if np.isnan(rmse_zlvl):
        diff = ds[rmse_field].values - ds_truth[rmse_field][2::3, 2::3].values
    else:
        diff = (ds[rmse_field][rmse_zlvl, :, :].values - 
                ds_truth[rmse_field][rmse_zlvl_truth, 2::3, 2::3].values)

    return diff


def compute_rmses(sims, truth, init_times,  rmse_field, rmse_zlvl, rmse_zlvl_truth, fcst_hr, debug=1):
    """
    Compute RMSEs over all times for both simulations
    """

    # Initialize output Datasets
    out_ds = init_out_ds(sims[list(sims.keys())[0]], init_times[0], rmse_field, rmse_zlvl)

    # Compute RMSEs
    ntot = [0, 0]
    for time in init_times:
        if debug > 0: print(f"Looping over {time.strftime('%Y%m%d%H')}")
        ds_truth = ds = open_ds(truth, time, offset=fcst_hr)
        for j, key in enumerate(list(sims.keys())):
            ds = open_ds(sims[key], time)
            if ds is not None:
                diff = compute_diff(ds, ds_truth, rmse_field, rmse_zlvl, rmse_zlvl_truth)
                out_ds[j]['RMSE'].values = out_ds[j]['RMSE'].values + (diff**2)
                ntot[j] = ntot[j] + 1
    for i in range(2):
        out_ds[i]['RMSE'].values = np.sqrt(out_ds[i]['RMSE'].values / ntot[i])

    return out_ds


if __name__ == '__main__':

    start = dt.datetime.now()
    print('\nStarting Program')
    print(f"{start.strftime('%Y%m%d %H:%M:%S')}\n")

    param = read_inputs(sys.argv[1])
    sims = param['sims']
    init_times = param['init_times']

    out_ds = compute_rmses(sims, 
                           param['truth'], 
                           init_times, 
                           param['rmse_field'], 
                           param['rmse_zlvl'], 
                           param['rmse_zlvl_truth'],
                           param['fcst_hr'])
    print(f"\nTime to compute RMSEs = {(dt.datetime.now() - start).total_seconds()}")

    # Plot output
    print('\nMaking Plots')
    fig = plt.figure(figsize=(7, 10))
    for i, ttl in enumerate(list(sims.keys())):
        # Note: Cannot use ccrs.LambertConformal() for map projection, otherwise program hangs 
        # when trying to make the contour plot
        out = pmd.PlotOutput([out_ds[i]], 'upp', fig, 3, 1, i+1, proj=ccrs.PlateCarree())
        out.contourf('RMSE', cntf_kw={'cmap':'inferno_r', 'extend':'max', 'levels':np.arange(0, 10.1, 0.5)})
        out.config_ax(grid=False)
        out.ax.set_title(ttl, size=16)

    diff = pmd.PlotOutput(out_ds, 'upp', fig, 3, 1, 3, proj=ccrs.PlateCarree())
    diff.plot_diff('RMSE', auto=False, cntf_kw={'cmap':'bwr', 'extend':'both',
                                                'levels':np.arange(-5, 5.1, 0.25)})
    diff.config_ax(grid=False)
    diff.ax.set_title('difference (top $-$ bottom)', size=16)
    
    plt.subplots_adjust(left=0.02, bottom=0.03, right=0.98, top=0.93, hspace=0.1)
    plt.suptitle(f"{init_times[0].strftime('%Y%m%d%H')} $-$ {init_times[-1].strftime('%Y%m%d%H')}", size=20)

    print('Saving plot')
    plt.savefig(param['out_fname'])
    plt.close()
    
    print('\nFinished!')
    print(f"Elapsed time = {(dt.datetime.now() - start).total_seconds()} s")


"""
End plot_hcrsxn_rmse_diff.py
"""
