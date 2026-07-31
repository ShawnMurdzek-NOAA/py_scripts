"""
Plot UAS Vertical Grid of Superobs in Height AGL

Used to create a plot in my extra slides for AMS Madison 2026

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

from pyDA_utils import bufr


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

bogus_file = '/work/noaa/wrfruc/murdzek/nature_run_spring/obs/uas_obs_150km/bogus_uas_csv/202204291900.rap.prepbufr.csv'
raw_file = '/work/noaa/wrfruc/murdzek/nature_run_spring/obs/uas_obs_150km/perf_uas_csv/202204291900.rap.fake.prepbufr.csv'
superob_file = '/work/noaa/wrfruc/murdzek/nature_run_spring/obs/uas_obs_150km/superob_uas/202204291900.rap.fake.prepbufr.csv'

superob_grid = '/work2/noaa/wrfruc/murdzek/src/osse_ob_creator/fix_data/RRFS_grid_mean_twice_gspacing.nc'

med_out_file = 'uas_vgrid.png'
lowest_superob_hist_file = 'lowest_superob_hist.png'


#---------------------------------------------------------------------------------------------------
# Create Plot
#---------------------------------------------------------------------------------------------------

# Read in grid for superobbing... this is just used for debugging
rrfs_ds = xr.open_dataset(superob_grid)
rrfs_lat = rrfs_ds['lat'].values
rrfs_lon = rrfs_ds['lon'].values
rrfs_sfc = rrfs_ds['HGT_SFC'].values

# Read in data
bogus_df = bufr.bufrCSV(bogus_file).df
bogus_df = bogus_df.loc[bogus_df['TYP'] == 136]
raw_df = bufr.bufrCSV(raw_file).df
raw_df = raw_df.loc[raw_df['TYP'] == 136]
superob_df = bufr.bufrCSV(superob_file).df
superob_df = superob_df.loc[superob_df['TYP'] == 136]

# Compute height agl
all_sid = np.unique(bogus_df['SID'].values)
nz = np.sum(superob_df['SID'] == all_sid[0]) + 2  # Add 2 as a buffer
hgt_agl = np.zeros([len(all_sid), nz]) * np.nan
nlev = np.zeros([3, len(all_sid)])
NR_sfc = np.zeros(len(all_sid))
superob_grid_sfc = np.zeros(len(all_sid))
for i, s in enumerate(all_sid):

    # Just trying to determine how many obs exist in each dataset
    # Bogus CSVs: 667 / UAS
    # Raw CSVs: 665 / UAS
    # Superob CSVs: 9-13 / UAS
    nlev[0, i] = np.sum(bogus_df['SID'] == s)
    nlev[1, i] = np.sum(raw_df['SID'] == s)
    nlev[2, i] = np.sum(superob_df['SID'] == s)

    # Determine height AGL of lowest UAS observation in raw_df
    # This height will likely be 6 m, which is the third vertical level in bogus_df
    # I think the first two levels from bogus_df are ignored b/c these are below the lowest model
    # vertical level
    min_ntb = raw_df.loc[raw_df['SID'] == s, 'ntb'].values[0]
    lowest_uas_agl = bogus_df.loc[np.logical_and(bogus_df['SID'] == s, bogus_df['ntb'] == min_ntb), 'ZOB'].values[0]

    # Determine height MSL of surface
    sfc_msl = raw_df.loc[raw_df['SID'] == s, 'ZOB'].values[0] - lowest_uas_agl
    NR_sfc[i] = sfc_msl

    # Determine heights AGL of superobs
    superob_z = superob_df.loc[superob_df['SID'] == s, 'ZOB'].values - sfc_msl
    hgt_agl[i, :len(superob_z)] = superob_z

    # Determine height MSL of surface in RRFSv1 using nearest neighbor interpolation
    uas_lon = superob_df.loc[superob_df['SID'] == s, 'XOB'].values[0]
    uas_lat = superob_df.loc[superob_df['SID'] == s, 'YOB'].values[0]
    i_rrfs, j_rrfs = np.unravel_index(np.argmin((uas_lon - rrfs_lon)**2 + (uas_lat - rrfs_lat)**2), rrfs_lat.shape)
    superob_grid_sfc[i] = rrfs_sfc[i_rrfs, j_rrfs]

med_hgt_agl = np.nanmedian(hgt_agl, axis=0)

# Plot median height AGL
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
for z in med_hgt_agl:
    ax.axhline(z, color='b', lw=2, ls='-')
ax.set_xlim([0, 1])
ax.set_ylim([0, np.nanmax(med_hgt_agl) + 5])
ax.set_ylabel('height AGL (m)', size=14)
ax.set_xticks([])
ax.set_title(f"Vertical Location of UAS Superobs\nLowest = {np.nanmin(med_hgt_agl):.1f} m, Highest = {np.nanmax(med_hgt_agl):.1f} m", size=18)
plt.savefig(med_out_file)

# Plot histogram of lowest UAS superob
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
ax.hist(hgt_agl[:, 0], bins=np.arange(0, 105, 5))
ax.grid()
ax.set_xlabel('height (m AGL)', size=15)
ax.set_ylabel('count', size=15)
ax.tick_params(axis='both', labelsize=12)
ax.set_title(f"Lowest UAS Superob in the 150-km Network\nMin = {np.nanmin(hgt_agl[:, 0]):.1f} m, Max = {np.nanmax(hgt_agl[:, 0]):.1f} m", size=18)
plt.savefig(lowest_superob_hist_file, dpi=500)


"""
End plot_uas_vertical_grid.py
"""
