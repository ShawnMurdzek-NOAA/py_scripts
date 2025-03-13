"""
Plot 2D Horizontal Cross Sections of Ensemble Statistics

Specifically, the min, 25th percentile, median, mean, 75th percentile, and max are plotted

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import sys
import argparse
import numpy as np
import copy
import xarray as xr
import matplotlib.pyplot as plt
import datetime as dt

import pyDA_utils.ensemble_utils as eu
import pyDA_utils.plot_model_data as pmd
import pyDA_utils.upp_postprocess as uppp


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
    Dictionary with parsed input arguments

    """

    parser = argparse.ArgumentParser(description='Script that plots horizontal cross section of \
                                                  ensemble statistics (min, 25th percentile, median, \
                                                  mean, 75th percentile, max) using UPP output')
    
    # Positional arguments
    parser.add_argument('file_tmpl', 
                        help='File template (including path to the file). Member number should be \
                              represented with the {mem} placeholder.',
                        type=str)
    
    parser.add_argument('n_ens', 
                        help='Number of ensemble members. Numbering is assumed to start at 1 and \
                              increase to n_ens',
                        type=int)

    # Optional arguments
    parser.add_argument('-f',
                        dest='field',
                        default='REFC_P0_L200_GLC0',
                        help='Field to plot. Can be 2D or 3D',
                        type=str)
    
    parser.add_argument('-l',
                        dest='level',
                        default=0,
                        help='Vertical level to plot if the field is 3D',
                        type=int)
    
    parser.add_argument('-c',
                        dest='colormap',
                        default='plasma',
                        help='Colormap used for filled contour plots',
                        type=str)
    
    parser.add_argument('-o',
                        dest='out_dir',
                        default='',
                        help='Output directory',
                        type=str)

    parser.add_argument('--outtag',
                        dest='out_tag',
                        default='',
                        help='Output file tag',
                        type=str)
    
    parser.add_argument('--vmin',
                        dest='vmin',
                        default=0,
                        help='vmin value for plotting',
                        type=float)
    
    parser.add_argument('--vmax',
                        dest='vmax',
                        default=70,
                        help='vmax value for plotting',
                        type=float)
    
    parser.add_argument('--minlat',
                        dest='minlat',
                        default=10,
                        help='Minimum latitude (deg N)',
                        type=float)
    
    parser.add_argument('--maxlat',
                        dest='maxlat',
                        default=53,
                        help='Maximum latitude (deg N)',
                        type=float)
    
    parser.add_argument('--minlon',
                        dest='minlon',
                        default=-134.1,
                        help='Minimum longitude (deg E)',
                        type=float)
    
    parser.add_argument('--maxlon',
                        dest='maxlon',
                        default=-60.9,
                        help='Maximum longitude (deg E)',
                        type=float)
    
    return parser.parse_args(argv)


def read_ens_members(param, verbose=False):
    """
    Read ensemble members

    Parameters
    ----------
    param : Namespace datastructure
        Input parameters
    verbose : boolean, optional
        Whether or not creation of eu.ensemble object is verbose
    
    Returns
    -------
    ens_obj : eu.ensemble() object
        Ensemble member output

    """

    fnames = {}
    for i in range(1, param.n_ens):
        fnames[f"mem{i:04d}"] = param.file_tmpl.format(mem=i)

    ens_obj = eu.ensemble(fnames,
                          verbose=verbose,
                          lat_limits=[param.minlat, param.maxlat],
                          lon_limits=[param.minlon, param.maxlon], 
                          zfield=None)

    return ens_obj


def preprocess_ens(ens_obj, field):
    """
    Compute derived fields, if necessary

    Parameters
    ----------
    ens_obj : eu.ensemble() object
        Ensemble member output
    field : string
        Field name

    Returns
    -------
    ens_obj : eu.ensemble() object
        Ensemble member output with derived fields computed

    """

    # Compute ceiling AGL
    if field in ['CEIL_LEGACY', 'CEIL_EXP1', 'CEIL_EXP2']:
        for m in ens_obj.mem_names:
            ens_obj.subset_ds[m] = uppp.compute_ceil_agl(ens_obj.subset_ds[m])

    # Compute wind speed and direction
    if field in ['WSPD', 'WDIR']:
        for m in ens_obj.mem_names:
            ens_obj.subset_ds[m] = uppp.compute_wspd_wdir(ens_obj.subset_ds[m])

    return ens_obj


def compute_ens_stats(ens_obj, param):
    """
    Compute required ensemble statistics

    Parameters
    ----------
    ens_obj : eu.ensemble() object
        Ensemble member output
    param : Namespace datastructure
        Input parameters
    
    Returns
    -------
    stat_ds : xr.Dataset
        Required statistics in a separate DataSet
    
    """

    # Extract template DataArray
    tmpl = copy.deepcopy(ens_obj.subset_ds[ens_obj.mem_names[0]][param.field])

    # Define functions for percentiles
    def p25(x, axis=0):
        return np.percentile(x, 25, axis=axis)
    def p75(x, axis=0):
        return np.percentile(x, 75, axis=axis)

    # Compute stats
    stat_dict = {}
    for name, fct in zip(['min', 'p25', 'mean', 'median', 'p75', 'max'],
                         [np.amin, p25, np.mean, np.median, p75, np.amax]):
        tmpl_copy = copy.deepcopy(tmpl)
        tmpl_copy.values = ens_obj.ens_stats(param.field, stat_fct=fct)
        stat_dict[name] = tmpl_copy

    # Save output to DataSet
    stat_ds = xr.Dataset(stat_dict)

    return stat_ds


def make_plot(stat_ds, param):
    """
    Make final plot

    Parameters
    ----------
    stat_ds : xr.Dataset
        Dataset containing ensemble statistics
    param : Namespace datastructure
        Input parameters
    
    Returns
    -------
    fig : matplotlib.figure object

    """

    # Create figures and get some preliminary data
    fig = plt.figure(figsize=(12, 8))
    stats = [k for k in stat_ds.keys()]
    lvls = np.linspace(param.vmin, param.vmax, 20)

    # Plot statistics
    for i, name in enumerate(stats):
        out = pmd.PlotOutput([stat_ds], 'upp', fig, 2, 3, i+1)
        if len(stat_ds[name].shape) == 2:
            out.contourf(name, 
                         cntf_kw={'cmap':param.colormap, 'levels':lvls},
                         cbar=False)
        else:
            out.contourf(name, 
                         cntf_kw={'cmap':param.colormap, 'levels':lvls},
                         ingest_kw={'zind':param.level},
                         cbar=False)
        out.config_ax(grid=False)
        out.ax_title(txt=name, size=12)

    # Add colorbar
    cb_ax = fig.add_axes([0.05, 0.07, 0.9, 0.03])
    cbar = plt.colorbar(out.cax, cax=cb_ax, orientation='horizontal', aspect=35, pad=0.1)
    cbar.set_label(f"{stat_ds[stats[0]].attrs['long_name']} ({stat_ds[stats[0]].attrs['units']})", size=12)
    
    # Format axes spacing
    plt.subplots_adjust(left=0.03, bottom=0.1, right=0.97, top=0.95, wspace=0.02)

    return fig


if __name__ == '__main__':

    start = dt.datetime.now()
    print('Starting plot_ens_stats_hcrsxn.py')
    print(f"Time = {start.strftime('%Y%m%d %H:%M:%S')}")

    param = parse_in_args(sys.argv[1:])

    print('Reading in ensemble members and computing derived fields...')
    ens_obj = read_ens_members(param)
    ens_obj = preprocess_ens(ens_obj, param.field)

    print('Computing statistics...')
    stat_ds = compute_ens_stats(ens_obj, param)

    print('Creating plot...')
    fig = make_plot(stat_ds, param)
    plt.savefig(f"{param.out_dir}/{param.field.split('_')[0]}_{param.level}_{param.out_tag}.png")

    print('Program Finished!')
    print(f"Elapsed time = {(dt.datetime.now() - start).total_seconds()} s")


"""
End plot_ens_stats_hcrsxn.py
"""
