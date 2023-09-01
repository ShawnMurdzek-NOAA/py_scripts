"""
Object-Based Comparison of MRMS and NR Composite Reflectivity Objects 

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
import scipy.ndimage as sn
import matplotlib.pyplot as plt
import sys
import glob


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Model (options are 'NR' or 'HRRR')
model = 'NR'

# Dates to use
eval_dates = [dt.datetime(2022, 2, 1) + dt.timedelta(days=i) for i in range(8)]

# NR data file path and subdirectories to use
if model == 'NR':
    NR_path = '/work2/noaa/wrfruc/murdzek/nature_run_winter/UPP'
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

# Domain (options: 'all', 'easternUS')
domain = 'all'

# Reflectivity contour used to define objects (dBZ)
ref_thres = 25

# Minimum reflectivity object size (number of 1 X 1 km^2 gridboxes)
min_size = 9

# Optional: Link to MRMS and NR masks (these should be .npy objects). Set to None if not being used
# Masks are created by the make_NR_MRMS_coverage_mask.py program
if model == 'NR':
    NR_mask_file = './NR_mask.npy'
elif model == 'HRRR':
    NR_mask_file = './HRRR_mask.npy'
MRMS_mask_file = './MRMS_mask.npy'

# Output file
out_file = './NR_cref_obj_%sdbz_%sminsize_%s_winter.png' % (ref_thres, min_size, domain)


#---------------------------------------------------------------------------------------------------
# Extract Data and Create Objects
#---------------------------------------------------------------------------------------------------

# Define necessary variables for each input field
if model == 'NR':
    NR_var = 'REFC_P0_L200_GLC0'
elif model == 'HRRR':
    NR_var = 'REFC_P0_L10_GLC0'
MRMS_var = 'MergedReflectivityQCComposite_P0_L102_GLL0'
MRMS_fname = 'MRMS_MergedReflectivityQCComposite'
MRMS_no_coverage = -999.

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
NR_obj = {}
MRMS_obj = {}
nMRMS = len(MRMS_years) * len(MRMS_offset)
for t in eval_times:
    NR_obj[t] = {}
    MRMS_obj[t] = {}
    for key in ['size', 'max_dbz']:
        NR_obj[t][key] = []
        MRMS_obj[t][key] = []
        for k in range(nMRMS):
            MRMS_obj[t][key].append([])

# Extract NR data
NR_mask = np.array([[np.nan]])
full_times = []
for d in eval_dates:
    d_str = d.strftime('%Y%m%d')
    try:
        ds = xr.open_dataset('%s/%s/cref_%s.nc' % (NR_path, d_str, d_str))
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
        full_t = dt.datetime.strptime(d_str + t, '%Y%m%d%H%M')
        try:
            time_idx = np.where(ds_timestamps == pd.Timestamp(full_t))[0][0]
        except IndexError:
            continue
        full_times.append(full_t)
        NR_data = NR_mask * ds[NR_var][time_idx, :, :].values
        NR_data_labeled, nlabels = sn.label(NR_data >= ref_thres)
        for label in range(1, nlabels):
            size = np.sum(NR_data_labeled == label)
            if size >= min_size:
                NR_obj[t]['size'].append(size)
                NR_obj[t]['max_dbz'].append(np.amax(NR_data * (NR_data_labeled == label)))

# Extract MRMS data
MRMS_mask = np.array([[np.nan]])
MRMS_years_all, MRMS_offset_all = np.meshgrid(MRMS_years, MRMS_offset)
MRMS_years_all = MRMS_years_all.ravel()
MRMS_offset_all = MRMS_offset_all.ravel()
for i, (y, o) in enumerate(zip(MRMS_years_all, MRMS_offset_all)):
    print()
    print('extracting MRMS data for year = %d, offset = %d' % (y, o))
 
    for t in full_times:
        MRMS_time = t + dt.timedelta(days=float(o))
        print('extracting MRMS data for %s' % MRMS_time.strftime('%m %d %H:%M'))
        fname_list = glob.glob('%s/%d/%d%s*%s*' % (MRMS_path, y, y, MRMS_time.strftime('%m%d-%H%M'), MRMS_fname))
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
        MRMS_data = MRMS_mask * ds[MRMS_var].values
        MRMS_data_labeled, nlabels = sn.label(MRMS_data >= ref_thres)
        for label in range(1, nlabels):
            size = np.sum(MRMS_data_labeled == label)
            if size >= min_size:
                MRMS_obj[hhmm]['size'][i].append(size)
                MRMS_obj[hhmm]['max_dbz'][i].append(np.amax(MRMS_data * (MRMS_data_labeled == label)))


#---------------------------------------------------------------------------------------------------
# Create Plots
#---------------------------------------------------------------------------------------------------

print()
print('Making plots...')

ncols = len(eval_times)
fig, axes = plt.subplots(nrows=3, ncols=ncols, height_ratios=[1, 2, 2], figsize=(12, 10), 
                         sharex='row', sharey='row')
plt.subplots_adjust(left=0.07, bottom=0.07, right=0.98, top=0.9, hspace=0.25, wspace=0.1)
labelsize = 14

for i, t in enumerate(eval_times):

    # First subplot: Number of objects
    ax = axes[0, i]
    nobj_MRMS = [len(sublist) for sublist in MRMS_obj[t]['size']]
    bplot = ax.boxplot(nobj_MRMS, vert=False, patch_artist=True, medianprops=dict(color='black'),
                       boxprops=dict(facecolor='lightcoral'))
    ax.plot(len(NR_obj[t]['size']), 1, 'ko')
    ax.set_xlabel('number of objects', size=labelsize)
    ax.grid(axis='x')

    # Second and third subplots: Histograms
    for j, (var, xlabel, bins, xscale) in enumerate(zip(['size', 'max_dbz'], 
                                                        ['object size (gridboxes)', 'max reflectivity (dBZ)'],
                                                        [np.linspace(4, 5000, 200), np.arange(25, 85, 5)],
                                                        ['log', 'linear'])):
        ax = axes[j+1, i]
        bin_ctrs = 0.5 * (bins[1:] + bins[:-1])

        NR_fraction = np.histogram(NR_obj[t][var], bins=bins)[0] / len(NR_obj[t][var]) 
        ax.plot(bin_ctrs, NR_fraction, 'k-', linewidth=2.5)

        MRMS_fraction = np.zeros([nMRMS, len(bin_ctrs)])
        for k in range(nMRMS):
            MRMS_fraction[k, :] = (np.histogram(MRMS_obj[t][var][k], bins=bins)[0] / 
                                   len(MRMS_obj[t][var][k]))
        MRMS_frac_pct = {}
        for pct in [0, 10, 25, 50, 75, 90, 100]:
            MRMS_frac_pct[pct] = np.nanpercentile(MRMS_fraction, pct, axis=0)
        ax.plot(bin_ctrs, MRMS_frac_pct[50], 'r-', linewidth=2.5)
        ax.fill_between(bin_ctrs, MRMS_frac_pct[25], MRMS_frac_pct[75], color='r', alpha=0.35)
        ax.fill_between(bin_ctrs, MRMS_frac_pct[10], MRMS_frac_pct[90], color='r', alpha=0.15)
        ax.plot(bin_ctrs, MRMS_frac_pct[0], 'r-', linewidth=0.75)
        ax.plot(bin_ctrs, MRMS_frac_pct[100], 'r-', linewidth=0.75)

        ax.set_xscale(xscale)
        ax.set_xlabel(xlabel, size=labelsize)
        ax.grid() 

    axes[0, i].set_title('%s UTC' % t, size=labelsize)

for i in range(2):
    axes[i+1, 0].set_ylabel('fraction', size=labelsize)

plt.suptitle('Reflectivity Objects ($Z_{thres}$ = %.1f dBZ, min size = %d gridboxes)' %  
             (ref_thres, min_size), size=18)
plt.savefig(out_file)
plt.close()


"""
End obj_based_mrms_cref_compare.py
"""
