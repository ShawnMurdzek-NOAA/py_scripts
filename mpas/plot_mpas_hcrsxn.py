"""
Plot an MPAS Horizontal Cross Section

This program should only be used to preliminary analysis. plots are created using the scatterplot
function, so it is possible that some small scale features could be washed out.

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import datetime as dt
import sys
import argparse
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import cartopy.crs as ccrs


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

    parser = argparse.ArgumentParser(description='Script that plots 2D cross sections from native \
                                                  MPAS output netCDF files using \
                                                  matplotlib.pyplot.scatter. Should only be used \
                                                  for quick plots of MPAS native output.')
    
    # Positional arguments
    parser.add_argument('atm_file', 
                        help='MPAS netCDF file containing atmospheric output',
                        type=str)
    
    parser.add_argument('fix_file', 
                        help='MPAS netCDF file containing coordinates of each cell',
                        type=str)

    # Optional arguments
    parser.add_argument('-f',
                        dest='field',
                        default='cldfrac',
                        help='Field to plot. Can be 2D or 3D',
                        type=str)
    
    parser.add_argument('-l',
                        dest='level',
                        default=0,
                        help='Vertical level to plot if the field is 3D',
                        type=int)
    
    parser.add_argument('-r',
                        dest='reduction',
                        default='none',
                        help='Function for reducing data in vertical dimension. \
                              Options: min, max, mean, median',
                        type=str)
    
    parser.add_argument('-c',
                        dest='colormap',
                        default='plasma',
                        help='Colormap used for filled contour plots',
                        type=str)
    
    parser.add_argument('-o',
                        dest='out_dir',
                        default='.',
                        help='Output directory',
                        type=str)
    
    parser.add_argument('-s',
                        dest='mark_size',
                        default=0.05,
                        help='Size of scatterplot dots',
                        type=float)

    parser.add_argument('--file2',
                        dest='atm_file2', 
                        default=np.nan,
                        help='Second MPAS netCDF file containing atmospheric output. \
                              Used for differences. Set to NaN to not use.',
                        type=str)

    parser.add_argument('--outtag',
                        dest='out_tag',
                        default='',
                        help='Output file tag',
                        type=str)
    
    parser.add_argument('--vmin',
                        dest='vmin',
                        default=None,
                        help='vmin value for plotting. Set to None to use defaults',
                        type=float)
    
    parser.add_argument('--vmax',
                        dest='vmax',
                        default=None,
                        help='vmax value for plotting. Set to None to use defaults',
                        type=float)
    
    parser.add_argument('--minlat',
                        dest='minlat',
                        default=20,
                        help='Minimum latitude (deg N)',
                        type=float)
    
    parser.add_argument('--maxlat',
                        dest='maxlat',
                        default=50,
                        help='Maximum latitude (deg N)',
                        type=float)
    
    parser.add_argument('--minlon',
                        dest='minlon',
                        default=-130,
                        help='Minimum longitude (deg E)',
                        type=float)
    
    parser.add_argument('--maxlon',
                        dest='maxlon',
                        default=-65,
                        help='Maximum longitude (deg E)',
                        type=float)
    
    return parser.parse_args(argv)


def read_process_data(param):
    """
    Read and process MPAS input files

    Parameters
    ----------
    param : argparse.Namespace
        Command-line arguments
    
    Returns
    -------
    data : 1D np.array
        Data to plot
    lat : 1D np.array
        Latitudes (deg N)
    lon : 1D np.array
        Longitudes (deg E)
    cbar_label : string
        Label for colorbar
    
    """

    # Read in netCDF files
    atm_ds = [xr.open_dataset(param.atm_file)]
    fix_ds = xr.open_dataset(param.fix_file)
    diff = False
    if type(param.atm_file2) == str:
        diff = True
        atm_ds.append(xr.open_dataset(param.atm_file2))

    # Extract atmospheric field
    # Be sure to remove the time dimension, which should always be 1
    data_ls = []
    for ds in atm_ds:
        if len(ds[param.field].shape) == 2:
            data_ls.append(ds[param.field].values[0, :])
        else:
            if param.reduction == 'none':
                data_ls.append(ds[param.field].values[0, :, param.level])
            elif param.reduction == 'max':
                data_ls.append(np.amax(ds[param.field].values[0, :, :], axis=1))
            elif param.reduction == 'min':
                data_ls.append(np.amin(ds[param.field].values[0, :, :], axis=1))
            elif param.reduction == 'mean':
                data_ls.append(np.mean(ds[param.field].values[0, :, :], axis=1))
            elif param.reduction == 'median':
                data_ls.append(np.median(ds[param.field].values[0, :, :], axis=1))
            else:
                raise ValueError(f"reduction method {param.reduction} is not supported")
        
    # Compute differences
    if diff:
        data = data_ls[0] - data_ls[1]
    else:
        data = data_ls[0]

    # Create colorbar label
    cbar_label = f"{atm_ds[0][param.field].long_name} ({atm_ds[0][param.field].units})"
    if diff: cbar_label = f"diff {cbar_label}"

    # Cell coordinate information
    lat = np.rad2deg(fix_ds['latCell'].values)
    lon = np.rad2deg(fix_ds['lonCell'].values)

    return data, lat, lon, cbar_label


def make_plot(data, lat, lon, cbar_label, param):
    """
    Create plot

    Parameters
    ----------
    data : 1D np.array
        Data to plot
    lat : 1D np.array
        Latitudes (deg N)
    lon : 1D np.array
        Longitudes (deg E)
    cbar_label : string
        Label for colorbar
    param : argparse.Namespace
        Command-line arguments
    
    Returns
    -------
    fig : plt.figure()
       Final plot

    """

    # Create figure and axes
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())

    # Make plot
    cax = ax.scatter(lon, lat, c=data, transform=ccrs.PlateCarree(),
                     vmin=param.vmin, vmax=param.vmax, cmap=param.colormap,
                     s=param.mark_size)
    cbar = plt.colorbar(cax, ax=ax, orientation='horizontal')
    cbar.set_label(cbar_label, size=12)
    
    # Set plotting extent
    ax.set_extent([param.minlon, param.maxlon, param.minlat, param.maxlat])

    # Add coastlines and state boundaries
    ax.coastlines('50m', linewidth=0.5, edgecolor='k')
    borders = cfeature.NaturalEarthFeature(category='cultural',
                                           scale='50m',
                                           facecolor='none',
                                           name='admin_1_states_provinces')
    ax.add_feature(borders, linewidth=0.25, edgecolor='gray')

    return fig


if __name__ == '__main__':

    start = dt.datetime.now()
    print('Starting plot_ens_stats_hcrsxn.py')
    print(f"Time = {start.strftime('%Y%m%d %H:%M:%S')}")

    param = parse_in_args(sys.argv[1:])

    print('Reading in MPAS data...')
    data, lat, lon, cbar_label = read_process_data(param)

    print('Creating plot...')
    fig = make_plot(data, lat, lon, cbar_label, param)
    if param.reduction == 'none':
        lvl_str = str(param.level)
    else:
        lvl_str = param.reduction
    plt.savefig(f"{param.out_dir}/{param.field}_{lvl_str}_{param.out_tag}.png")

    print('Program Finished!')
    print(f"Elapsed time = {(dt.datetime.now() - start).total_seconds()} s")


"""
End plot_mpas_hcrsxn.py
"""