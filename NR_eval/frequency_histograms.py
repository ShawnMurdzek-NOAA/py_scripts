"""
Frequency Histograms Comparing the Nature Run to Observed MRMS Output

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


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Dates to use
eval_dates = [dt.datetime(2022, 4, 29) + dt.timedelta(days=i) for i in range(8)]

# NR data file path and subdirectories to use
NR_path = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_spring/UPP'

# MRMS data file path
MRMS_path = '/mnt/lfs4/BMC/wrfruc/murdzek/real_obs/mrms'

# MRMS years to use
MRMS_years = np.arange(2018, 2023)

# Times to evaluate (strings, HHMM)
eval_times = ['0000', '0600', '1200', '1800']

# Field to evaluate (options: 'cref', 'precip1hr')
field = 'precip1hr'

# Domain (options: 'all', 'easternUS')
# This option does not work yet!
domain = 'all'

# Output file
out_file = './NR_precip1hr_eval_lower_prate.png'


#---------------------------------------------------------------------------------------------------
# Create Histograms
#---------------------------------------------------------------------------------------------------

# Useful functions for NR variable transforms
def precip_kgpm2_to_mm(x):
    return x * 1e3 / 997.

# Define necessary variables for each input field
if field == 'cref':
    NR_var = 'REFC_P0_L200_GLC0'
    NR_transform = None
    MRMS_var = ['MergedReflectivityQCComposite_P0_L102_GLL0']
    MRMS_fname = ['MRMS_MergedReflectivityQCComposite']
    MRMS_no_coverage = -999.
    bins = np.arange(5, 76, 5)
    xlabel = 'composite reflectivity (dBZ)'
    yscale = 'linear'
elif field == 'precip1hr':
    NR_var = 'APCP_P8_L1_GLC0_acc'
    NR_transform = precip_kgpm2_to_mm
    MRMS_var = ['VAR_209_6_37_P0_L102_GLL0', 'GaugeCorrQPE01H_P0_L102_GLL0']
    MRMS_fname = ['MRMS_MultiSensor_QPE_01H_Pass2', 'MRMS_GaugeCorr_QPE_01H']
    MRMS_no_coverage = -3.
    bins = np.arange(0.5, 15, 0.5)
    xlabel = '1-hr total precip (mm)'
    yscale = 'log'

# Spatial domain limits
if domain == 'all':
    lat_lim = []
    lon_lim = []
elif domain == 'easternUS':
    lat_lim = [30, 45]
    lon_lim = [-100, 80]

# Initialize dictionaries
NR_counts = {}
NR_total_counts = {}
NR_total_pts = {}
for t in eval_times:
    NR_total_counts[t] = np.zeros(bins.size-1)
    NR_total_pts[t] = 0

MRMS_counts = {}
MRMS_total_counts = {}
MRMS_total_pts = {}
for y in MRMS_years:
    MRMS_counts[y] = {}
    MRMS_total_counts[y] = {}
    MRMS_total_pts[y] = {}
    for t in eval_times:
        MRMS_total_counts[y][t] = np.zeros(bins.size-1)
        MRMS_total_pts[y][t] = 0

# Extract NR data
for d in eval_dates:
    d_str = d.strftime('%Y%m%d')
    try:
        ds = xr.open_dataset('%s/%s/%s_%s.nc' % (NR_path, d_str, field, d_str))
    except FileNotFoundError:
        continue

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
            NR_data = NR_transform(ds[NR_var][time_idx, :, :].values)
        else:    
            NR_data = ds[NR_var][time_idx, :, :].values
        NR_counts[full_time] = np.histogram(NR_data, bins=bins)[0]
        NR_total_counts[t] = NR_total_counts[t] + NR_counts[full_time]
        NR_total_pts[t] = NR_total_pts[t] + NR_data.size

# Extract MRMS data
for y in MRMS_years:
    print()
    print('extracting MRMS data for %d' % y)
    for t in NR_counts.keys():
        print('extracting MRMS data for %s' % t.strftime('%m %d %H:%M'))
        for n, f in enumerate(MRMS_fname):
            fname_list = glob.glob('%s/%d/%d%s*%s*' % (MRMS_path, y, y, t.strftime('%m%d-%H%M'), f))
            if len(fname_list) > 0:
                break
        if len(fname_list) == 0:
            print('MRMS data for %d-%s is missing!' % (y, t.strftime('%m-%d %H:%M')))
            continue
        ds = xr.open_dataset(fname_list[0], engine='pynio')
        MRMS_counts[y][t] = np.histogram(ds[MRMS_var[n]].values, bins=bins)[0]
        hhmm = t.strftime('%H%M')
        MRMS_total_counts[y][hhmm] = MRMS_total_counts[y][hhmm] + MRMS_counts[y][t]
        MRMS_total_pts[y][hhmm] = MRMS_total_pts[y][hhmm] + np.sum(ds[MRMS_var[n]].values > MRMS_no_coverage)
   
# Plot results
nrows = int(np.floor(np.sqrt(len(eval_times))))
ncols = int(np.ceil(len(eval_times) / nrows))
figsize = (3 + 2*ncols, 3 + 2*nrows) 
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True, figsize=figsize)
bin_ctrs = bins[:-1] + 0.5 * (bins[1] - bins[0])
for i, t in enumerate(eval_times):
    ax = axes[int(i / ncols), i % ncols]

    ax.bar(bin_ctrs, NR_total_counts[t] / NR_total_pts[t], width=(bins[1] - bins[0]), color='b',
           edgecolor='k', linewidth=0.25)    
    for y in MRMS_years:
        ax.plot(bin_ctrs, MRMS_total_counts[y][t] / MRMS_total_pts[y][t], 'r-')
    
    ax.set_title('%s UTC' % t, size=14)
    ax.set_yscale(yscale)
    ax.grid() 

for i in range(ncols):
    axes[-1, i].set_xlabel(xlabel, size=12)
for i in range(nrows):
    axes[i, 0].set_ylabel('fraction of gridpoints in domain', size=12)

plt.suptitle('NR Frequencies (bars) and MRMS Frequencies (lines)\n7-Day Composites', size=16)
plt.savefig(out_file)
plt.close()


"""
End frequency_histograms.py
"""
