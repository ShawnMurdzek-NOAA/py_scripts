"""
Frequency Histograms Comparing the Nature Run to Observed MRMS Output

This script uses "extracted" netCDF files created by extract_wrfnat_fields.py

6-hr precip amounts are computed by ../misc/compute_precip6hr.py

shawn.s.murdzek@noaa.gov
Date Created: 8 March 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import glob
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import sys
import pickle


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Model (options are 'NR' or 'HRRR')
#model = 'NR'
model = sys.argv[1]

# Dates to use
#eval_dates = [dt.datetime(2022, 2, 1) + dt.timedelta(days=i) for i in range(8)]
eval_dates = [dt.datetime(2022, 4, 29) + dt.timedelta(days=i) for i in range(8)]

# NR data file path and subdirectories to use
if model == 'NR':
    NR_path = '/work2/noaa/wrfruc/murdzek/nature_run_spring/UPP'
elif model == 'HRRR':
    NR_path = '/work2/noaa/wrfruc/murdzek/HRRR_data'

# MRMS data file path
MRMS_path = '/work2/noaa/wrfruc/murdzek/real_obs/mrms'

# MRMS years to use
MRMS_years = np.arange(2016, 2024)

# Option to use days from MRMS before and after eval_dates. Offsets here are days relative to the 
# first date in eval_dates
MRMS_offset = [-7, 0, 7]

# Times to evaluate (strings, HHMM)
eval_times = ['0000', '0600', '1200', '1800']

# Field to evaluate (options: 'cref', 'precip1hr', 'precip6hr')
#field = 'precip1hr'
field = sys.argv[2]

# Domain (options: 'all', 'easternUS')
#domain = 'all'
domain = sys.argv[3]

# Optional: Link to MRMS and NR masks (these should be .npy objects). Set to None if not being used
# Masks are created by the make_NR_MRMS_coverage_mask.py program
if model == 'NR':
    NR_mask_file = './NR_mask.npy'
elif model == 'HRRR':
    NR_mask_file = './HRRR_mask.npy'
MRMS_mask_file = './MRMS_mask.npy'

# Option to zoom into the smallest precip rates for precip1hr or the highest reflectivities for cref
#zoom = 0
zoom = bool(int(sys.argv[4]))

# Option to save/use output from a pickle file
# If use_pickle is True, then the script will attempt to read the pickle file specified. If the file
# is not found, that file will be written to.
use_pickle = True
pickle_fname = './%s_%s_%s_spring.pkl' % (model, field, domain)

# Output file
#out_file = './NR_precip1hr_eval_all.png'
out_file = sys.argv[5]


#---------------------------------------------------------------------------------------------------
# Create Histograms
#---------------------------------------------------------------------------------------------------

# Useful functions for NR variable transforms
def precip_kgpm2_to_mm(x):
    return x * 1e3 / 997.

start_time = dt.datetime.now()

# Try to read from pickle file
if use_pickle:
    try:
        with open(pickle_fname, 'rb') as handle:
            all_data = pickle.load(handle)
        NR_total_counts = all_data['NR_total_counts']
        NR_total_pts = all_data['NR_total_pts']
        MRMS_freq = all_data['MRMS_freq']
        yscale = all_data['yscale']
        xlabel = all_data['xlabel']
    except FileNotFoundError:
        pickle_avail = False
else:
    pickle_avail = False

if not pickle_avail:

    # Define necessary variables for each input field
    if field == 'cref':
        if model == 'NR':
            NR_var = 'REFC_P0_L200_GLC0'
        elif model == 'HRRR':
            NR_var = 'REFC_P0_L10_GLC0'
        NR_transform = None
        MRMS_var = ['MergedReflectivityQCComposite_P0_L102_GLL0']
        MRMS_fname = ['MRMS_MergedReflectivityQCComposite']
        MRMS_no_coverage = -999.
        if zoom:
            bins = np.arange(30, 76, 5)
        else:
            bins = np.arange(5, 76, 5)
        xlabel = 'composite reflectivity (dBZ)'
        yscale = 'linear'
    elif field == 'precip1hr':
        NR_var = 'APCP_P8_L1_GLC0_acc'
        NR_transform = precip_kgpm2_to_mm
        MRMS_var = ['VAR_209_6_37_P0_L102_GLL0', 'GaugeCorrQPE01H_P0_L102_GLL0']
        MRMS_fname = ['MRMS_MultiSensor_QPE_01H_Pass2', 'MRMS_GaugeCorr_QPE_01H']
        MRMS_no_coverage = -3.
        if zoom:
            bins = np.arange(0.5, 15, 0.5)
        else:
            bins = np.arange(1, 150, 5)
        xlabel = '1-hr total precip (mm)'
        yscale = 'log'
    elif field == 'precip6hr':
        NR_var = 'APCP_P8_L1_GLC0_acc'
        NR_transform = precip_kgpm2_to_mm
        MRMS_var = ['VAR_209_6_39_P0_L102_GLL0', 'GaugeCorrQPE06H_P0_L102_GLL0']
        MRMS_fname = ['MRMS_MultiSensor_QPE_06H_Pass2', 'MRMS_GaugeCorr_QPE_06H']
        MRMS_no_coverage = -3.
        if zoom:
            bins = np.arange(0.5, 15, 0.5)
        else:
            bins = np.arange(1, 200, 5)
        xlabel = '6-hr total precip (mm)'
        yscale = 'log'

    # Add 0 to MRMS_offset if empty
    if len(MRMS_offset) == 0:
        MRMS_offset = [0]

    # Spatial domain limits
    if domain == 'all':
        lat_lim = [5, 70]
        lon_lim = [-150, -40]
    elif domain == 'easternUS':
        lat_lim = [5, 70]
        lon_lim = [-100, -40]

    # Load NR and MRMS masks
    if NR_mask_file != None:
        NR_mask_external = np.load(NR_mask_file)
    else:
        NR_mask_external = 1
    if MRMS_mask_file != None:
        MRMS_mask_external = np.load(MRMS_mask_file)
    else:
        MRMS_mask_external = 1

    # Initialize dictionaries
    NR_counts = {}
    NR_total_counts = {}
    NR_total_pts = {}
    for t in eval_times:
        NR_total_counts[t] = np.zeros(bins.size-1)
        NR_total_pts[t] = 0

    MRMS_total_counts = {}
    MRMS_total_pts = {}
    MRMS_freq = {}
    n_MRMS = len(MRMS_years) * len(MRMS_offset)
    MRMS_years_all, MRMS_offset_all = np.meshgrid(MRMS_years, MRMS_offset)
    MRMS_years_all = MRMS_years_all.ravel()
    MRMS_offset_all = MRMS_offset_all.ravel()
    for t in eval_times:
        MRMS_total_counts[t] = np.zeros([bins.size-1, n_MRMS])
        MRMS_freq[t] = np.zeros([bins.size-1, n_MRMS])
        MRMS_total_pts[t] = np.zeros(n_MRMS)

    # Extract NR data
    NR_mask = np.array([[np.nan]])
    for d in eval_dates:
        d_str = d.strftime('%Y%m%d')
        try:
            ds = xr.open_dataset('%s/%s/%s_%s.nc' % (NR_path, d_str, field, d_str))
        except FileNotFoundError:
            continue

        # Create mask based on desired domain
        if np.isnan(NR_mask[0, 0]):
            NR_lat = ds['gridlat_0'].values
            NR_lon = ds['gridlon_0'].values
            NR_mask = NR_mask_external * ((NR_lat >= lat_lim[0]) * (NR_lat <= lat_lim[1]) * 
                                          (NR_lon >= lon_lim[0]) * (NR_lon <= lon_lim[1]))

        # Unfortunately, the numpy datetime objects in ds and the datetime datetime objects in 
        # eval_times cannot be easily compared, so we'll convert both of them to pd.Timestamp objects
        ds_timestamps = np.empty(ds['time'].size, dtype=object)
        for j, t in enumerate(ds['time'].values):
            ds_timestamps[j] = pd.Timestamp(t)

        print('NR extracting data for %s' % d_str)
        for t in eval_times:
            full_time = dt.datetime.strptime(d_str + t, '%Y%m%d%H%M')
            try:
                time_idx = np.where(ds_timestamps == pd.Timestamp(full_time))[0][0]
            except IndexError:
                continue
            if NR_transform != None:
                NR_data = NR_mask * NR_transform(ds[NR_var][time_idx, :, :].values)
            else:    
                NR_data = NR_mask * ds[NR_var][time_idx, :, :].values
            NR_counts[full_time] = np.histogram(NR_data, bins=bins)[0]
            NR_total_counts[t] = NR_total_counts[t] + NR_counts[full_time]
            NR_total_pts[t] = NR_total_pts[t] + np.sum(NR_mask)

    # Extract MRMS data
    MRMS_mask = np.array([[np.nan]])
    for i, (y, o) in enumerate(zip(MRMS_years_all, MRMS_offset_all)):
        print()
        print('extracting MRMS data for year = %d, offset = %d' % (y, o))
        for t in NR_counts.keys():
            MRMS_time = t + dt.timedelta(days=float(o))
            print('extracting MRMS data for %s' % MRMS_time.strftime('%m %d %H:%M'))
            for n, f in enumerate(MRMS_fname):
                fname_list = glob.glob('%s/%d/%d%s*%s*' % (MRMS_path, y, y, MRMS_time.strftime('%m%d-%H%M'), f))
                if len(fname_list) > 0:
                    break
            if len(fname_list) == 0:
                print('MRMS data for %d-%s is missing!' % (y, MRMS_time.strftime('%m-%d %H:%M')))
                continue
            ds = xr.open_dataset(fname_list[0], engine='pynio')
 
            # Create mask for MRMS data
            if np.isnan(MRMS_mask[0, 0]):
                MRMS_lon, MRMS_lat = np.meshgrid(ds['lon_0'].values - 360., ds['lat_0'].values)
                MRMS_mask = MRMS_mask_external * ((MRMS_lat >= lat_lim[0]) * (MRMS_lat <= lat_lim[1]) * 
                                                  (MRMS_lon >= lon_lim[0]) * (MRMS_lon <= lon_lim[1]))

            hhmm = t.strftime('%H%M')
            MRMS_total_counts[hhmm][:, i] = (np.histogram(MRMS_mask * ds[MRMS_var[n]].values, bins=bins)[0] + 
                                             MRMS_total_counts[hhmm][:, i])
            MRMS_total_pts[hhmm][i] = MRMS_total_pts[hhmm][i] + np.sum(np.logical_or(MRMS_mask, ds[MRMS_var[n]].values > MRMS_no_coverage))

        for hhmm in eval_times:
            MRMS_freq[hhmm][:, i] = MRMS_total_counts[hhmm][:, i] / MRMS_total_pts[hhmm][i]

    # Save output to pickle file for use later
    if use_pickle:
        all_data = {}
        all_data['NR_total_counts'] = NR_total_counts
        all_data['NR_total_pts'] = NR_total_pts
        all_data['MRMS_freq'] = MRMS_freq
        all_data['yscale'] = yscale
        all_data['xlabel'] = xlabel
        with open(pickle_fname, 'wb') as handle:
            pickle.dump(all_data, handle)


#---------------------------------------------------------------------------------------------------   
# Plot results
#---------------------------------------------------------------------------------------------------

nrows = int(np.floor(np.sqrt(len(eval_times))))
ncols = int(np.ceil(len(eval_times) / nrows))
figsize = (3 + 2*ncols, 3 + 2*nrows) 
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True, figsize=figsize)
bin_ctrs = bins[:-1] + 0.5 * (bins[1] - bins[0])
for i, t in enumerate(eval_times):
    ax = axes[int(i / ncols), i % ncols]

    ax.plot(bin_ctrs, NR_total_counts[t] / NR_total_pts[t], 'k-', linewidth=2.5) 

    MRMS_freq_pct = {}
    for pct in [0, 10, 25, 50, 75, 90, 100]:
        MRMS_freq_pct[pct] = np.nanpercentile(MRMS_freq[t], pct, axis=1)
    
    ax.plot(bin_ctrs, MRMS_freq_pct[50], 'r-', linewidth=2.5)
    ax.fill_between(bin_ctrs, MRMS_freq_pct[25], MRMS_freq_pct[75], color='r', alpha=0.35)
    ax.fill_between(bin_ctrs, MRMS_freq_pct[10], MRMS_freq_pct[90], color='r', alpha=0.15)
    ax.plot(bin_ctrs, MRMS_freq_pct[0], 'r-', linewidth=0.75)
    ax.plot(bin_ctrs, MRMS_freq_pct[100], 'r-', linewidth=0.75)
 
    ax.set_title('%s UTC' % t, size=14)
    ax.set_yscale(yscale)
    ax.grid() 

for i in range(ncols):
    axes[-1, i].set_xlabel(xlabel, size=12)
for i in range(nrows):
    axes[i, 0].set_ylabel('fraction of gridpoints in domain', size=12)

plt.suptitle('%s Frequencies (black) and MRMS Frequencies (red)\n7-Day Composites' % model, size=16)
plt.savefig(out_file)
plt.close()

print('elapsed time = %.2f min' % ((dt.datetime.now() - start_time).total_seconds() / 60))


"""
End frequency_histograms.py
"""
