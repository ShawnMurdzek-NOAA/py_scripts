"""
Plot JEDI obs locations

Gemini Pro 3.1 with edits from Shawn Murdzek
"""

import argparse
import sys
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

def extract_locations(filepath):
    """Extract longitude and latitude from an IODA v2+ netCDF file."""
    try:
        # IODA v2+ stores metadata variables in the 'MetaData' group
        with xr.open_dataset(filepath, group='MetaData') as ds:
            lons = ds['longitude'].values
            lats = ds['latitude'].values
            return lons, lats
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return np.array([]), np.array([])

def main():
    parser = argparse.ArgumentParser(description="Plot observation locations from JEDI IODA netCDF files (v2+).")
    parser.add_argument("file_list", help="Text file containing a list of netCDF file paths (one per line)")
    parser.add_argument("-o", "--output", default="obs_locations_plot.png", 
                        help="Output image file name (default: obs_locations_plot.png)")
    parser.add_argument("-t", "--title", default="JEDI Assimilated Observation Locations", 
                        help="Image title")
    args = parser.parse_args()

    with open(args.file_list, 'r') as f:
        file_paths = [line.strip() for line in f if line.strip()]

    all_lons = []
    all_lats = []

    for filepath in file_paths:
        lons, lats = extract_locations(filepath)
        if lons.size > 0 and lats.size > 0:
            all_lons.extend(lons)
            all_lats.extend(lats)

    total_obs = len(all_lons)
    
    # Convert lat and lon lists to arrays
    all_lons = np.array(all_lons)
    all_lats = np.array(all_lats)

    if total_obs == 0:
        print("No observation locations found in the provided files.")
        sys.exit(1)

    # Setup the Cartopy plot with Lambert Conformal projection
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())

    # Dynamically set plotting bounds based on data
    #pad = 1
    #ax.set_extent([np.amin(all_lons) - pad,
    #               np.amax(all_lons) + pad,
    #               np.amin(all_lats) - pad,
    #               np.amin(all_lats) + pad], crs=ccrs.PlateCarree())
    ax.set_extent([-124, -68, 21, 54], crs=ccrs.PlateCarree())

    # Add map features
    ax.add_feature(cfeature.LAND, facecolor='whitesmoke')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linewidth=0.8, linestyle='--')
    ax.add_feature(cfeature.LAKES, facecolor='lightblue', edgecolor='black', linewidth=0.5)

    # Add state borders
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none')
    ax.add_feature(states_provinces, linewidth=0.5, linestyle=':')

    # Plot the observations (Data is in lat/lon, so transform is PlateCarree)
    ax.scatter(all_lons, all_lats, s=15, color='red', alpha=0.7, edgecolor='k',
               transform=ccrs.PlateCarree(), zorder=7)

    ax.set_title(f"{args.title}\nTotal Observations: {total_obs:,}", 
                 fontsize=14, pad=10)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"Successfully plotted {total_obs:,} observations. Saved to {args.output}.")

if __name__ == "__main__":
    main()


"""
End plot_jedi_ob_locs.py
"""
