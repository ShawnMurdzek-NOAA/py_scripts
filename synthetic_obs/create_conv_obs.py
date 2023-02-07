"""
Create Synthetic Conventional Observations from WRF Nature Run Output

Conventional obs are written into a new CSV file. UAS obs can then be added later.

Prepbufr CSV files for real observations can be created using GSI-utils/bin/prepbufr_decode_csv.x.
To use this utility, do the following:

1. Copy the prepbufr file you wish to decode to GSI-utils/bin and rename it 'prepbufr'
2. Run ./prepbufr_decode_csv.x
3. Move resulting CSV file to wherever you choose

This code uses UPP output on WRF native levels (wrfnat). This is likely more efficient than using 
the wrfred files, which can take some time (and a lot of memory) to create.

Resources:

Hera - An interactive session on a regular node is sufficient
Jet - 10 GB of memory on xjet is not enough

shawn.s.murdzek@noaa.gov
Date Created: 22 November 2022
Environment: adb_graphics (Jet) or pygraf (Hera)
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import numpy as np
import datetime as dt
import json
import math

import create_ob_utils as cou 
import bufr


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

work = '/mnt/lfs4'

# Directory containing wrfnat output from UPP
wrf_dir = work + '/BMC/wrfruc/murdzek/nature_run_spring_v2/'

# Directory containing real prepbufr CSV output
bufr_dir = './'

# File containing conventional observation error characteristics (use GSI errtable file)
error_fname = work + '/BMC/wrfruc/murdzek/sample_real_obs/errtable.rrfs'

# Observation platforms to use (aka subsets, same ones used by BUFR)
ob_platforms = ['ADPUPA', 'AIRCAR', 'AIRCFT', 'PROFLR', 'ADPSFC', 'SFCSHP', 'MSONET', 'GPSIPW']

# Output directory for synthetic prepbufr CSV output
fake_bufr_dir = work + '/BMC/wrfruc/murdzek/nature_run_spring_v2/synthetic_obs/'

# Start and end times for prepbufrs. Step is in min
bufr_start = dt.datetime(2022, 4, 29, 12)
bufr_end = dt.datetime(2022, 4, 29, 13)
bufr_step = 120

# Start and end times for wrfnat UPP output. Step is in min
wrf_start = dt.datetime(2022, 4, 29, 12, 0)
wrf_end = dt.datetime(2022, 4, 29, 13, 0)
wrf_step = 15

# Option for debugging output (0 = none, 1 = some, 2 = a lot)
debug = 2


#---------------------------------------------------------------------------------------------------
# Create Synthetic Observations
#---------------------------------------------------------------------------------------------------

# Open up obs error file
errors = cou.read_ob_errors(error_fname)

# List to save time spent on each BUFR entry
if debug > 1:
    entry_times = []

ntimes = int((bufr_end - bufr_start) / dt.timedelta(minutes=bufr_step) + 1)
for i in range(ntimes):
    t = bufr_start + dt.timedelta(minutes=(i*bufr_step))
    start = dt.datetime.now()
    print()
    print('t = %s' % t.strftime('%Y-%m-%d %H%M%S'))
    print('start time = %s' % start.strftime('%Y%m%d %H:%M:%S'))
    bufr_fname = bufr_dir + t.strftime('/%Y%m%d%H%M.rap.prepbufr.csv')
    bufr_csv = bufr.bufrCSV(bufr_fname)

    # Only keep platforms if we are creating synthetic obs for them
    obs = bufr_csv.df['subset'].unique()
    if debug > 0:
        print('total # of BUFR entries = %d' % len(bufr_csv.df))
        print('available BUFR obs =', obs)
    for s in obs:
        if s not in ob_platforms:
            bufr_csv.df.drop(index=np.where(bufr_csv.df['subset'] == s)[0], inplace=True)
            bufr_csv.df.reset_index(drop=True, inplace=True)

    # Remove obs outside of the range of wrfnat files
    hr_min = ((wrf_start - t).days * 24) + ((wrf_start - t).seconds / 3600.)
    hr_max = ((wrf_end - t).days * 24) + ((wrf_end - t).seconds / 3600.)
    bufr_csv.df.drop(index=np.where(np.logical_or(bufr_csv.df['DHR'] < hr_min, 
                                                  bufr_csv.df['DHR'] > hr_max))[0], inplace=True)
    bufr_csv.df.reset_index(drop=True, inplace=True)

    # Open wrfnat files
    # Ideally this would be done using xr.open_mfdataset(), but the wrfnat files lack a time 
    # dimension, which makes Xarray unable to combine them (and using concatenate requires a lot
    # of memory, so that is not a viable option)
    hr_start = math.floor(bufr_csv.df['DHR'].min()*4) / 4
    hr_end = math.ceil(bufr_csv.df['DHR'].max()*4) / 4 + (wrf_step / 60.)
    wrf_ds = {}
    wrf_hr = []
    for hr in np.arange(hr_start, hr_end + (wrf_step / 60.), wrf_step / 60.):
        wrf_t = t + dt.timedelta(hours=hr)
        print(wrf_dir + wrf_t.strftime('wrfnat_%Y%m%d%H%M.grib2'))
        wrf_ds[hr] = xr.open_dataset(wrf_dir + wrf_t.strftime('wrfnat_%Y%m%d%H%M.grib2'),
                                     engine='pynio')
        wrf_hr.append(hr)
    
    print('time to open GRIB files = %.2f s' % (dt.datetime.now() - start).total_seconds())

    # Remove obs outside of the spatial domain of the wrfnat files
    latmin = wrf_ds[wrf_hr[0]]['gridlat_0'].attrs['corners'].min()
    latmax = wrf_ds[wrf_hr[0]]['gridlat_0'].attrs['corners'].max()
    lonmin = wrf_ds[wrf_hr[0]]['gridlon_0'].attrs['corners'].min() + 360.
    lonmax = wrf_ds[wrf_hr[0]]['gridlon_0'].attrs['corners'].max() + 360.
    bufr_csv.df.drop(index=np.where(np.logical_or(bufr_csv.df['XOB'] < lonmin, 
                                                  bufr_csv.df['XOB'] > lonmax))[0], inplace=True)
    bufr_csv.df.reset_index(drop=True, inplace=True)
    bufr_csv.df.drop(index=np.where(np.logical_or(bufr_csv.df['YOB'] < latmin, 
                                                  bufr_csv.df['YOB'] > latmax))[0], inplace=True)
    bufr_csv.df.reset_index(drop=True, inplace=True)

    if debug > 0:
        print('# BUFR entries remaining = %d' % len(bufr_csv.df))

    # Create output DataFrame
    out_df = bufr_csv.df.copy()
    drop_idx = []

    # Loop over each observation
    wrf_hr = np.array(wrf_hr)
    for j in range(len(out_df)):
    #for j in range(14):
        
        if debug > 1:
            time1 = dt.datetime.now()
            print()

        subset = out_df.iloc[j]
     
        # Determine WRF hour right before observation and weight for temporal interpolation
        ihr = np.where((wrf_hr - subset['DHR']) <= 0)[0][-1]
        twgt = (subset['DHR'] - wrf_hr[ihr]) / (wrf_hr[ihr+1] - wrf_hr[ihr]) 

        if debug > 1:
            time2 = dt.datetime.now()
            print('done determining twgt (%.6f s)' % (time2 - time1).total_seconds())

        # Interpolate horizontally
        xi, yi = cou.wrf_coords(subset['YOB'], subset['XOB'] - 360., wrf_ds[wrf_hr[ihr]], ftype='UPP')   
        if np.isnan(xi):
            # ob location lies outside model domain
            drop_idx.append(j)
            continue
     
        if debug > 1:
            time3 = dt.datetime.now()
            print('done determining xi, yi (%.6f s)' % (time3 - time2).total_seconds())

        if subset['subset'] in ['ADPSFC', 'SFCSHP', 'MSONET']:
            # Surface obs. Don't need to interpolate vertically

            if debug > 1:
                print('2D field: nmsg = %d, subset = %s, TYP = %d' % (subset['nmsg'], 
                                                                      subset['subset'], 
                                                                      subset['TYP']))
                print('twgt = %.3f, xi = %.3f, yi = %.3f' % (twgt, xi, yi))
                print('DHR = %.3f, XOB = %.3f, YOB = %.3f' % (subset['DHR'], subset['XOB'], 
                                                              subset['YOB']))

            # Interpolate temporally (winds must be done separately b/c both 10-m and 80-m winds
            # are reported in UPP). Also adjust ELV and ZOB to match Nature Run terrain
            obs_name = ['POB', 'QOB', 'TOB', 'PRSS', 'ZOB', 'ELV']
            wrf_name = ['PRES_P0_L1_GLC0', 'SPFH_P0_L103_GLC0', 'TMP_P0_L103_GLC0', 
                        'PRES_P0_L1_GLC0', 'HGT_P0_L1_GLC0', 'HGT_P0_L1_GLC0']
            for o, m in zip(obs_name, wrf_name):
                if not np.isnan(subset[o]):
                    if debug > 1:
                        time5 = dt.datetime.now()
                    out_df.loc[j, o] = (twgt * wrf_ds[wrf_hr[ihr]][m].interp(ygrid_0=yi, xgrid_0=xi) +
                                        (1-twgt) * wrf_ds[wrf_hr[ihr+1]][m].interp(ygrid_0=yi, xgrid_0=xi))
                    if debug > 1:
                        time6 = dt.datetime.now()
                        print('finished interp for %s (%.6f s)' % (o, (time6 - time5).total_seconds()))

            obs_name = ['UOB', 'VOB']
            wrf_name = ['UGRD_P0_L103_GLC0', 'VGRD_P0_L103_GLC0']
            for o, m in zip(obs_name, wrf_name):
                if not np.isnan(subset[o]):
                    if debug > 1:
                        time5 = dt.datetime.now()
                    out_df.loc[j, o] = (twgt * wrf_ds[wrf_hr[ihr]][m][0, :, :].interp(ygrid_0=yi, xgrid_0=xi) +
                                        (1-twgt) * wrf_ds[wrf_hr[ihr+1]][m][0, :, :].interp(ygrid_0=yi, xgrid_0=xi))
                    if debug > 1:
                        time6 = dt.datetime.now()
                        print('finished interp for %s (%.6f s)' % (o, (time6 - time5).total_seconds()))

        elif subset['subset'] == 'GPSIPW':
            # GPS-Derived precipitable water

            if debug > 1:
                print('PW field: nmsg = %d, subset = %s, TYP = %d' % (subset['nmsg'], 
                                                                      subset['subset'], 
                                                                      subset['TYP']))
                print('twgt = %.3f, xi = %.3f, yi = %.3f' % (twgt, xi, yi))

            # Interpolate temporally
            m = 'PWAT_P0_L200_GLC0'
            if not np.isnan(subset['PWO']):
                if debug > 1:
                    time5 = dt.datetime.now()
                out_df.loc[j, 'PWO'] = (twgt * wrf_ds[wrf_hr[ihr]][m].interp(ygrid_0=yi, xgrid_0=xi) +
                                        (1-twgt) * wrf_ds[wrf_hr[ihr+1]][m].interp(ygrid_0=yi, xgrid_0=xi))
                if debug > 1:
                    time6 = dt.datetime.now()
                    print('finished interp for PWO (%.6f s)' % (time6 - time5).total_seconds())

        else:
            # 3D fields. Need to interpolate in vertical
            if debug > 1:
                print('3D field: nmsg = %d, subset = %s, TYP = %d' % (subset['nmsg'], 
                                                                      subset['subset'], 
                                                                      subset['TYP']))
                print('twgt = %.3f, xi = %.3f, yi = %.3f' % (twgt, xi, yi))

            zi = np.zeros(2)
            for k in [ihr, ihr+1]:
                z1d = wrf_ds[wrf_hr[k]]['HGT_P0_L105_GLC0'].interp(ygrid_0=yi, xgrid_0=xi)
                z0i = np.where((z1d - subset['ZOB']) <= 0)[0][-1]
                zi[k] = z0i + (subset['ZOB'] - z1d[z0i]) / (z1d[z0i+1] - z1d[z0i]) 

            if debug > 1:
                time4 = dt.datetime.now()
                print('done determining zi [%.2f, %.2f] (%.6f s)' % (zi[0], zi[1], (time4 - time3).total_seconds()))

            # Interpolate temporally
            obs_name = ['POB', 'QOB', 'TOB', 'UOB', 'VOB']
            wrf_name = ['PRES_P0_L105_GLC0', 'SPFH_P0_L105_GLC0', 'TMP_P0_L105_GLC0', 'UGRD_P0_L105_GLC0',
                        'VGRD_P0_L105_GLC0']
            for o, m in zip(obs_name, wrf_name):
                if not np.isnan(subset[o]):
                    if debug > 1:
                        time5 = dt.datetime.now()
                    out_df.loc[j, o] = (twgt * wrf_ds[wrf_hr[ihr]][m].interp(lv_HYBL1=zi[0], ygrid_0=yi, xgrid_0=xi) +
                                        (1-twgt) * wrf_ds[wrf_hr[ihr+1]][m].interp(lv_HYBL1=zi[1], ygrid_0=yi, xgrid_0=xi))
                    if debug > 1:
                        time6 = dt.datetime.now()
                        print('finished interp for %s (%.6f s)' % (o, (time6 - time5).total_seconds()))

        # Add errors

        if debug > 1:
            print('total time = %.6f s' % (dt.datetime.now() - time1).total_seconds())
            entry_times.append((dt.datetime.now() - time1).total_seconds())

    # Drop rows that we skipped
    out_df.drop(index=drop_idx, inplace=True)
    out_df.reset_index(drop=True, inplace=True)
    bufr_csv.df.drop(index=drop_idx, inplace=True)
    bufr_csv.df.reset_index(drop=True, inplace=True)

    # Convert to proper units
    out_df['QOB'] = out_df['QOB'] * 1e6
    out_df['POB'] = out_df['POB'] * 1e-2
    out_df['TOB'] = out_df['TOB'] - 273.15
    out_df['PWO'] = (out_df['PWO'] / 997.) * 1000.
    out_df['ELV'] = np.int64(out_df['ELV'])

    # Write output DataFrame to a CSV file
    bufr.df_to_csv(out_df, fake_bufr_dir + t.strftime('/%Y%m%d%H%M.fake.prepbufr.csv'))
    bufr.df_to_csv(bufr_csv.df, fake_bufr_dir + t.strftime('/%Y%m%d%H%M.real_red.prepbufr.csv'))

if debug > 1:
    entry_times = np.array(entry_times)
    print()
    print('avg time per entry = %.6f s' % np.mean(entry_times))


"""
End create_conv_obs.py
"""
