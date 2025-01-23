"""
Histograms of Lapse Rates from FV3 NetCDF Output

Inputs
------
sys.argv[1] : Input YAML file

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import yaml
import datetime as dt


#---------------------------------------------------------------------------------------------------
# Main Code
#---------------------------------------------------------------------------------------------------

def parse_inputs(fname):
    """
    Parse inputs from a YAML file

    Parameters
    ----------
    fname : string
        Input YAML file
    
    Returns
    -------
    param : dictionary
        Input parameters

    """

    with open(fname, 'r') as fptr:
        param = yaml.safe_load(fptr)
    
    return param


def compute_fv3_z(ds):
    """
    Compute model level heights AGL

    This code mimics fv_diag_column.F90 from the FV3 source code

    Parameters
    ----------
    ds : xr.Dataset
        FV3 output from fv_core.res.tile1.nc

    Returns
    -------
    ds : xr.Dataset
        FV3 output with a height AGL field

    """
    
    rgrav = 1 / 9.8
    ds['hgt'] = ds['DZ'].copy()
    ds['hgt'].values[0, -1, :, :] = ds['phis'].values[0, :, :]*rgrav - 0.5*ds['DZ'].values[0, -1, :, :]
    for i in range(ds['hgt'].shape[1]-2, -1, -1):
        ds['hgt'].values[0, i, :, :] = ds['hgt'].values[0, i+1, :, :] - 0.5*(ds['DZ'].values[0, i, :, :] + ds['DZ'].values[0, i+1, :, :])
    ds['hgt'].attrs = {'long_name': 'height AGL', 'units':'m'}

    return ds


def read_ens_mem(param):
    """
    Read in all ensemble members

    Parameters
    ----------
    param : dictionary
        Input parameters

    Returns
    -------
    ens_dict : dictionary
        Ensemble output

    """

    ens_dict = {}
    for key in param['data']:
        ens_dict[key] = {}
        for n in range(1, param['gen_param']['nens']+1):
            fname = param['data'][key]['fmt'].format(num=n)
            if os.path.isfile(fname):
                ens_dict[key][n] = xr.open_dataset(fname)
                ens_dict[key][n] = compute_fv3_z(ens_dict[key][n])
            else:
                print(f"File does not exist, skipping: {fname}")
    
    return ens_dict


def compute_lapse_rates_FV3(ds):
    """
    Compute lapse rates for a single FV3 output dataset

    Parameters
    ----------
    ds : xr.Dataset
        FV3 output
    
    Returns
    -------
    ds : xr.Dataset
        FV3 output with the field 'lapse_rate' added

    """

    # Compute lapse rate
    lr = ((ds['T'].values[:, 1:, :, :] - ds['T'].values[:, :-1, :, :]) /
          (1e-3 * (ds['hgt'].values[:, :-1, :, :] - ds['hgt'].values[:, 1:, :, :])))
    
    # Add lapse rate to Dataset
    lr_dict = {'data':lr, 
               'dims':['Time', 'zaxis_2', 'yaxis_2', 'xaxis_1'],
               'coords':{'zaxis_2':{'dims':'zaxis_2', 'data':np.arange(1, 65)}}}
    da = xr.DataArray.from_dict(lr_dict)
    ds['lapse_rate'] = da

    return ds


def compute_lapse_rate_all(ens_dict):
    """
    Compute lapse rates for all ensemble output

    Parameters
    ----------
    ens_dict : dictionary
        Ensemble output. Key is dataset name, value is another dictionary of ensemble output
    
    Returns
    -------
    ens_dict : dictionary
        Ensemble output with lapse rates included

    """

    for key in ens_dict.keys():
        for n in ens_dict[key].keys():
            ens_dict[key][n] = compute_lapse_rates_FV3(ens_dict[key][n])
    
    return ens_dict


def plot_lr_hist_1d(ens_dict, param):
    """
    Plot histograms of lapse rates

    Parameters
    ----------
    ens_dict : dictionary
        Ensemble output
    param : dictionary
        Plotting parameters

    Returns
    -------
    fig : matplotlib.figure
        Figure containing histograms of lapse rates

    """

    # Create figure and subplots
    nrows = param['gen_param']['subplot_kw']['nrows']
    ncols = param['gen_param']['subplot_kw']['ncols']
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                             figsize=param['gen_param']['figsize'],
                             sharex=True, sharey=True)
    plt.subplots_adjust(**param['gen_param']['subplots_adjust'])
    
    # Plot histogram
    bins = np.array(param['gen_param']['lr_bins'])
    bin_plt = 0.5*(bins[1:] + bins[:-1])
    for n in range(1, param['gen_param']['nens']+1):
        ax = axes[int((n-1) / ncols), (n-1)%ncols]
        for key in param['data']:
            if n in ens_dict[key].keys():
                hist, _ = np.histogram(ens_dict[key][n]['lapse_rate'].values, bins=bins)
                ax.plot(bin_plt, hist, label=key, **param['data'][key]['plot_kw'])
        ax.set_yscale('log')
        ax.set_xlim([bins[0], bins[-1]])
        ax.axvline(9.8, ls='--', c='k')
        ax.set_title(f"mem{n:04}", size=12)
        if param['gen_param']['grid']:
            ax.grid()
    
    # Add annotations
    for i in range(nrows): 
        axes[i, 0].set_ylabel('counts', size=12)
    for i in range(ncols): 
        axes[nrows-1, i].set_xlabel('lapse rate (K km$^{-1}$)', size=12)
    #axes[0, 0].legend()
    plt.suptitle(param['gen_param']['title'], size=16)
    
    return fig


def plot_lr_hgt_hist_2d(ens_dict, param):
    """
    Plot 2D histograms of lapse rates and heights AGL

    Parameters
    ----------
    ens_dict : dictionary
        Ensemble output
    param : dictionary
        Plotting parameters

    Returns
    -------
    fig : matplotlib.figure
        Figure containing histograms of lapse rates

    """

    # Create figure and subplots
    nrows = param['gen_param']['subplot_kw']['nrows']
    ncols = param['gen_param']['subplot_kw']['ncols']
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                             figsize=param['gen_param']['figsize'],
                             sharex=True, sharey=True)
    plt.subplots_adjust(**param['gen_param']['subplots_adjust'])
    
    # Plot histogram
    lr_bins = np.array(param['gen_param']['lr_bins'])
    z_bins = np.arange(param['gen_param']['z_bins']['start'],
                       param['gen_param']['z_bins']['stop'],
                       param['gen_param']['z_bins']['step'])
    for n in range(1, param['gen_param']['nens']+1):
        ax = axes[int((n-1) / ncols), (n-1)%ncols]
        hist = {}
        sims = list(param['data'].keys())
        for key in sims:
            if n in ens_dict[key].keys():
                avg_hgt = 1e-3*0.5*(ens_dict[key][n]['hgt'].values[0, 1:, :, :] + ens_dict[key][n]['hgt'].values[0, :-1, :, :])
                hist[key], _, _ = np.histogram2d(ens_dict[key][n]['lapse_rate'].values.flatten(),
                                                 avg_hgt.flatten(),
                                                 bins=[lr_bins, z_bins])
        if len(hist) == 2:
            cax = ax.pcolormesh(lr_bins, z_bins, (hist[sims[1]] - hist[sims[0]]).T, 
                                shading='flat', **param['gen_param']['hist2d_pcolormesh_kw'])
        ax.set_xlim([lr_bins[0], lr_bins[-1]])
        ax.set_ylim([z_bins[0], z_bins[-1]])
        ax.axvline(9.8, ls='--', c='k')
        ax.set_title(f"mem{n:04}", size=12)
    
    # Add annotations
    for i in range(nrows): 
        axes[i, 0].set_ylabel('height (km AGL)', size=12)
    for i in range(ncols): 
        axes[nrows-1, i].set_xlabel('lapse rate (K km$^{-1}$)', size=12)
    plt.suptitle(f"{param['gen_param']['title']}", size=16)
    cbar = plt.colorbar(cax, ax=axes, orientation='horizontal', pad=0.12)
    cbar.set_label(f"difference in counts ({sims[1]} - {sims[0]})", size=12)
    
    return fig


if __name__ == '__main__':

    start = dt.datetime.now()

    # Read in parameters
    print('Reading in input parameters...')
    param = parse_inputs(sys.argv[1])

    # Read in ensemble members
    print('Reading in ensemble members...')
    ens_dict = read_ens_mem(param)

    # Compute lapse rates
    print('Computing lapse rates...')
    ens_dict = compute_lapse_rate_all(ens_dict)

    # Plot 1D histogram
    print('Plotting 1D histograms...')
    fig = plot_lr_hist_1d(ens_dict, param)
    fname = f"lr_1D_hist_{param['gen_param']['out_tag']}.png"
    plt.savefig(fname)

    # Plot 2D histogram
    print('Plotting 2D histograms...')
    fig = plot_lr_hgt_hist_2d(ens_dict, param)
    fname = f"lr_hgt_2D_hist_{param['gen_param']['out_tag']}.png"
    plt.savefig(fname)

    print()
    print(f'total elapsed time = {(dt.datetime.now() - start).total_seconds()} s')


"""
End plot_fv3_lapse_rate_histograms.py
"""