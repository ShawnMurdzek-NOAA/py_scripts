"""
Plot output from NCAR's 3-km MPAS nature run using pyDAmonitor

Before running, you must load the pyDAmonitor environment:

source pyDAmonitor/ush/load_pyDAmonitor.sh

Shawn Murdzek
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import sys
import argparse
import datetime as dt
import uxarray as ux
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import copy


#---------------------------------------------------------------------------------------------------
# Main Program
#---------------------------------------------------------------------------------------------------

# Dictionary with plotting parameters
plot_param = {'refl10cm_max' : {'cmin': 5, 'cmax': 80, 'cmap': 'turbo'},
              'ceil' : {'cmin': 0, 'cmax': 5000, 'cmap': 'plasma'},
              'cldfrac_max' : {'cmin': 0, 'cmax': 1, 'cmap': 'plasma_r'}}

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
                                                  MPAS output using uxarray and pyDAmonitor. \
                                                  Plotted fields are chosen specifically for the \
                                                  NCAR nature run for the cloud DA project')
    
    # Positional arguments
    parser.add_argument('atm_file', 
                        help='MPAS netCDF file containing atmospheric output',
                        type=str)
    
    parser.add_argument('fix_file', 
                        help='MPAS netCDF file containing coordinates of each cell',
                        type=str)

    # Optional arguments
    parser.add_argument('--tag',
                        dest='tag',
                        default='',
                        help='Tag used for title and output file name',
                        type=str)

    return parser.parse_args(argv)


def compute_hgt_agl(fix_fname):
    """
    Compute heights AGL using the MPAS fix file
    """

    ds = xr.open_dataset(fix_fname)
    z = 0.5 * (ds['zgrid'].values[:, 1:] + ds['zgrid'].values[:, :-1])

    return z - ds['ter'].values[:, np.newaxis]


def compute_ceil(ds, z_agl):
    """
    Compute cloud ceiling field
    """

    # Compute ceiling
    thres = 0.5
    tmp = copy.deepcopy(z_agl)
    tmp[ds['cldfrac'][0, :, :] < thres] = np.nan
    ceil = np.nanmin(tmp, axis=1)

    # Add ceiling to dataset
    ds['ceil'] = copy.deepcopy(ds['t2m'])
    ds['ceil'].values[0, :] = ceil
    ds['ceil'].attrs.update({'units': 'm AGL',
                             'long_name': 'cloud ceiling (cldfrac >= 0.5'})

    return ds


def compute_max_cldfrac(ds, z_agl, thres=3657):
    """
    Compute maximum cldfrac below thres (m)
    """

    # Compute max cloud fraction below threshold elevation AGL
    tmp = ds['cldfrac'][0, :, :].values
    tmp[z_agl > thres] = np.nan
    cldfrac_max = np.nanmax(tmp, axis=1)

    # Add max cldfrac to dataset
    ds['cldfrac_max'] = copy.deepcopy(ds['t2m'])
    ds['cldfrac_max'].values[0, :] = cldfrac_max
    ds['cldfrac_max'].attrs.update({'units': 'unitless',
                                    'long_name': f"max cldfrac below {thres} m AGL"})

    return ds


def plot_mpas_2d_raster(ds, param, tag=''):
    """
    Plot 2D MPAS fields as using uxarray's to_raster method
    """

    # Define lakes for base map
    lakes = cfeature.NaturalEarthFeature(category='physical',
                                         name='lakes',
                                         scale='50m',
                                         facecolor='none')

    # Loop over each field to plot
    for key in param:

        # Configure figure and axes
        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.LambertConformal()},
                               figsize=(8, 6),
                               constrained_layout=True)
        ax.set_extent([-112, -66, 16, 56])

        # Rasterize and plot with imshow
        uxvar = ds[key].isel(Time=0)
        raster, pixel_mapping = uxvar.to_raster(ax=ax, pixel_ratio=3, return_pixel_mapping=True)
        raster[raster < param[key]['cmin']] = np.nan
        raster[raster > param[key]['cmax']] = np.nan
        cax = ax.imshow(raster, 
                        cmap=param[key]['cmap'], vmin=param[key]['cmin'], vmax=param[key]['cmax'],
                        origin='lower', extent=ax.get_xlim() + ax.get_ylim())

        # Add underlying base map
        ax.add_feature(cfeature.COASTLINE, linewidth=0.75, edgecolor='k')
        ax.add_feature(cfeature.STATES, linewidth=0.25, edgecolor='gray')
        ax.add_feature(lakes, linewidth=0.25, edgecolor='gray')

        # Add colorbar
        cbar = plt.colorbar(cax, ax=ax, orientation='horizontal', aspect=30)
        cbar.set_label(f"{uxvar.attrs['long_name']} ({uxvar.attrs['units']})", size=12)

        # Add title
        ax.set_title(f"{tag} {key}", size=16)

        # Save to png
        if tag == '':
            fname = f"{key}.png"
        else:
            fname = f"{key}_{tag}.png"

        plt.savefig(fname)

    return None


if __name__ == '__main__':

    start = dt.datetime.now()
    print('Starting plot_mpas_nr.py')
    print(f"Time = {start.strftime('%Y%m%d %H:%M:%S')}")

    in_param = parse_in_args(sys.argv[1:])

    # Open file using uxarray
    uxds = ux.open_dataset(in_param.fix_file, in_param.atm_file)

    # Compute height AGL
    hgt_agl = compute_hgt_agl(in_param.fix_file)

    # Compute cloud fields
    print('Computing cloud fields')
    uxds = compute_ceil(uxds, hgt_agl)
    uxds = compute_max_cldfrac(uxds, hgt_agl)

    # Make plots
    print('Making plots')
    plot_mpas_2d_raster(uxds, plot_param, tag=in_param.tag)

    print('Program Finished!')
    print(f"Elapsed time = {(dt.datetime.now() - start).total_seconds()} s")


"""
End plot_mpas_nr.py
"""
