"""
Plot Synthetic MRMS Output

Plots are intentionally similar to pygraf for easier comparisons to model data

Shawn Murdzek
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import sys
import argparse
import datetime as dt
import glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature


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

    parser = argparse.ArgumentParser(description='Script that plots MRMS composite reflectivity \
                                                  using a series of MRMS output files at fixed \
                                                  heights ASL.')
    
    # Positional arguments
    parser.add_argument('mrms_path', 
                        help='Directory containing MRMS GRIB2 output',
                        type=str)

    parser.add_argument('mrms_time', 
                        help='MRMS analysis time in YYYYMMDDHHMMSS format',
                        type=str)
    
    # Optional arguments
    parser.add_argument('--out_fname',
                        dest='out_fname',
                        default='mrms_out.png',
                        help='Output file name',
                        type=str)

    return parser.parse_args(argv)


def compute_mrms_cref(path, time):
    """
    Compute MRMS composite reflectivity using files in a given path with a given timestamp
    """

    # Create list of files
    files = glob.glob(f"{path}/*{time[:8]}-{time[8:]}_processed.grib2")
    files.sort()

    # Open a sample file
    field = 'CONUSMergedReflectivityQC_P0_L102_GLL0'
    ds = xr.open_dataset(files[0], engine='pynio')

    # Save lat, lon, and reflectivity
    lat = ds['lat_0'].values
    lon = ds['lon_0'].values
    refl_3d = np.zeros([len(files), len(lat), len(lon)])
    refl_3d[0, :, :] = ds[field].values
    ds.close()

    # Loop over remaining files
    for i, f in enumerate(files[1:]):
        ds = xr.open_dataset(f, engine='pynio')
        refl_3d[i+1, :, :] = ds[field].values
        ds.close()

    # Compute composite reflectivity
    cref = np.amax(refl_3d, axis=0)

    return cref, lat, lon


def cref_colors(clevs):
    """
    Reflectivity colormap from pygraf/adb_graphics/specs.py
    """

    rgb_vals = np.array([[1.,         1.,         1.,         1.        ],
                         [0.,         0.9254902,  0.9254902,  1.        ],
                         [0.00392157, 0.62745098, 0.96470588, 1.        ],
                         [0.,         0.,         0.96470588, 1.        ],
                         [0.,         1.,         0.,         1.        ],
                         [0.,         0.78431373, 0.,         1.        ],
                         [0.,         0.56470588, 0.,         1.        ],
                         [1.,         1.,         0.,         1.        ],
                         [0.90588235, 0.75294118, 0.,         1.        ],
                         [1.,         0.56470588, 0.,         1.        ],
                         [1.,         0.,         0.,         1.        ],
                         [0.83921569, 0.,         0.,         1.        ],
                         [0.75294118, 0.,         0.,         1.        ],
                         [1.,         0.,         1.,         1.        ],
                         [0.6,        0.33333333, 0.78823529, 1.        ],
                         [1.,         1.,         1.,         1.        ]])

    return ListedColormap(rgb_vals)


def plot_field_pcolormesh_conus(lat, lon, field, clevs=np.arange(5, 76, 5), title='', cmap='gist_ncar'):
    """
    Plot field over CONUS using CartoPy
    """

    # Create figure and axes
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
    #ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    # Plot data
    #cax = ax.contourf(lon, lat, field, clevs, transform=ccrs.PlateCarree(), cmap=cmap)
    norm = BoundaryNorm(clevs, cmap.N)
    cax = ax.pcolormesh(lon, lat, field, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
    cbar = fig.colorbar(cax, ax=ax, orientation='horizontal')

    # Add annotations
    ax.set_extent([-120, -72, 22, 50])
    ax.coastlines(linewidth=0.5, edgecolor='k')
    ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='k', facecolor='none')
    ax.add_feature(cfeature.LAKES, linewidth=0.5, edgecolor='k', facecolor='none')
    ax.set_title(title, size=18)

    plt.tight_layout()

    return fig


if __name__ == '__main__':

    start = dt.datetime.now()
    print('Starting plot_syn_mrms_data.py')
    print(f"Time = {start.strftime('%Y%m%d %H:%M:%S')}")

    in_param = parse_in_args(sys.argv[1:])

    # Compute composite reflectivity
    mrms_cref, mrms_lat, mrms_lon = compute_mrms_cref(in_param.mrms_path, in_param.mrms_time)

    # Debugging: Only plot REF from a single level
    #fname = '/work2/noaa/wrfruc/RRFS2_RETRO_DATA/May2024/reflectivity/MergedReflectivityQC_05.00_20240508-210042.grib2'
    #title = 'Real MRMS: 5 km ASL, 20240508-210042'
    #fname = '/work/noaa/wrfruc/jdduda/RadarNext/processed_reflectivity/20220501-010000/MergedReflectivityQC_05.00_20220501-010000_processed.grib2'
    #title = 'Syn MRMS: 5 km ASL, 20220501-010000'
    #ds = xr.open_dataset(fname, engine='pynio')
    #mrms_cref = ds['CONUSMergedReflectivityQC_P0_L102_GLL0'].values
    #mrms_lat = ds['lat_0'].values
    #mrms_lon = ds['lon_0'].values

    # Create plot
    contours = np.arange(0, 81, 5)
    cmap = cref_colors(contours)
    t_str = f"{in_param.mrms_time[:8]} {in_param.mrms_time[8:12]} UTC"
    title = f"MRMS Composite Reflectivity (dBZ, shaded)\n{t_str}"
    fig = plot_field_pcolormesh_conus(mrms_lat, mrms_lon, mrms_cref, 
                                      clevs=contours,
                                      cmap=cmap,
                                      title=title)
    plt.savefig(in_param.out_fname)
    plt.close()

    print('Program Finished!')
    print(f"Elapsed time = {(dt.datetime.now() - start).total_seconds()} s")


"""
End plot_syn_mrms_data.py
"""
