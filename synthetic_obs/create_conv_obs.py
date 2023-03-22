"""
Create Perfect Simulated Conventional Observations from WRF Nature Run Output

Conventional obs are written into a new CSV file. Instrument errors and other obs (e.g.,  UAS) can 
then be added later. Conventional observations include the following: ADPSFC, SFCSHP, MSONET, 
GPSIPW, ADPUPA, AIRCAR, AIRCFT. Radiosonde drift is NOT accounted for in this program. Use
create_adpupa_obs.py to create realistic radiosonde obs with drift included. Obs in this script
are created used pure interpolation (linear in x, y, and t and logaritmic in p).

Prepbufr CSV files for real observations can be created using GSI-utils/bin/prepbufr_decode_csv.x.
To use this utility, do the following:

1. Copy the prepbufr file you wish to decode to GSI-utils/bin and rename it 'prepbufr'
2. Run ./prepbufr_decode_csv.x
3. Move resulting CSV file to wherever you choose

Passed Arguments:
    argv[1] = Directory containing WRF UPP output
    argv[2] = Directory containing real prepBUFR CSV files
    argv[3] = Output directory for simulated prepBUFR CSV files
    argv[4] = Time for first prepbufr file (YYYYMMDDHH)
    argv[5] = Time for first UPP file (YYYYMMDDHH)
    argv[6] = Time for last UPP file (YYYYMMDDHH)

shawn.s.murdzek@noaa.gov
Date Created: 22 November 2022
Environment: adb_graphics (Jet) or pygraf (Hera, Orion, Hercules)
RAM Required: 25 GB
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import numpy as np
import datetime as dt
import math
import os
import sys

import create_ob_utils as cou 
import bufr
import map_proj as mp


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Directory containing wrfnat output from UPP
#wrf_dir = '/work2/noaa/wrfruc/murdzek/nature_run_spring/UPP/'
wrf_dir = sys.argv[1]

# Directory containing real prepbufr CSV output
#bufr_dir = '/work2/noaa/wrfruc/murdzek/real_obs/obs_rap_csv/'
bufr_dir = sys.argv[2]

# Observation platforms to use (aka subsets, same ones used by BUFR)
ob_platforms = ['ADPUPA', 'AIRCAR', 'AIRCFT', 'ADPSFC', 'SFCSHP', 'MSONET', 'GPSIPW']

# Output directory for synthetic prepbufr CSV output
#fake_bufr_dir = '/work2/noaa/wrfruc/murdzek/nature_run_spring/synthetic_obs_csv/conv/'
fake_bufr_dir = sys.argv[3]

# Start and end times for prepbufrs. Step is in min
bufr_start = dt.datetime.strptime(sys.argv[4], '%Y%m%d%H')
bufr_end = bufr_start + dt.timedelta(hours=1)
#bufr_start = dt.datetime(2022, 4, 29, 12)
#bufr_end = dt.datetime(2022, 4, 29, 13)
bufr_step = 120

# Start and end times for wrfnat UPP output. Step is in min
wrf_start = dt.datetime.strptime(sys.argv[5], '%Y%m%d%H')
wrf_end = dt.datetime.strptime(sys.argv[6], '%Y%m%d%H')
#wrf_start = dt.datetime(2022, 4, 29, 12, 0)
#wrf_end = dt.datetime(2022, 4, 29, 15, 0)
wrf_step = 15

# Option to interpolate height obs (ZOB) for AIRCAR and AIRCFT platforms
# Heights reported by aircraft are calculated by integrating the hydrostatic balance eqn assuming
# the US Standard Atmosphere. Therefore, the heights for the simulated obs should be the same as the
# real obs b/c the pressures are the same.
interp_z_aircft = False

# Option for debugging output (0 = none, 1 = some, 2 = a lot)
debug = 2

# Option to interpolate (lat, lon) coordinates for surface obs (ADPSFC, SFCSHP, MSONET)
# Helpful for debugging, but should usually be set to False b/c it increases runtime
interp_latlon = False


#---------------------------------------------------------------------------------------------------
# Create Synthetic Observations
#---------------------------------------------------------------------------------------------------

# Timing
begin = dt.datetime.now()

# List to save time spent on each BUFR entry
if debug > 1:
    entry_times2d = []
    entry_timesp1 = []
    entry_timesp2 = []
    entry_times3d = []

# Convert wrf_step to decimal hours for ease of use
wrf_step_dec = wrf_step / 60.

ntimes = int((bufr_end - bufr_start) / dt.timedelta(minutes=bufr_step) + 1)
for i in range(ntimes):
    t = bufr_start + dt.timedelta(minutes=(i*bufr_step))
    start_loop = dt.datetime.now()
    print()
    print('t = %s' % t.strftime('%Y-%m-%d %H:%M:%S'))
    print('start time = %s' % start_loop.strftime('%Y%m%d %H:%M:%S'))
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
    hr_start = math.floor(bufr_csv.df['DHR'].min()*4) / 4
    hr_end = math.ceil(bufr_csv.df['DHR'].max()*4) / 4 + (2*wrf_step_dec)
    print('min/max WRF hours = %.2f, %.2f' % (hr_min, hr_max))
    print('min/max BUFR hours = %.2f, %.2f' % (bufr_csv.df['DHR'].min(), bufr_csv.df['DHR'].max()))
    wrf_ds = {}
    wrf_hr = np.arange(hr_start, hr_end, wrf_step_dec)
    for hr in wrf_hr:
        wrf_t = t + dt.timedelta(hours=hr)
        print(wrf_dir + wrf_t.strftime('/%Y%m%d/wrfnat_%Y%m%d%H%M.grib2'))
        wrf_ds[hr] = xr.open_dataset(wrf_dir + wrf_t.strftime('/%Y%m%d/wrfnat_%Y%m%d%H%M.grib2'), 
                                     engine='pynio')

    print('time to open GRIB files = %.2f s' % (dt.datetime.now() - start_loop).total_seconds())
    
    # Extract size of latitude and longitude grids
    shape = wrf_ds[0]['gridlat_0'].shape
    imax = shape[0] - 2
    jmax = shape[1] - 2
    
    # Compute (x, y) coordinates of obs using a Lambert Conformal projection
    # Remove obs outside of the wrfnat domain
    if debug > 0:
        start_map_proj = dt.datetime.now()
        print('Performing map projection with obs...')

    bufr_csv.df['xlc'], bufr_csv.df['ylc'] = mp.ll_to_xy_lc(bufr_csv.df['YOB'], bufr_csv.df['XOB'] - 360.)
    bufr_csv.df['i0'] = np.int32(np.floor(bufr_csv.df['ylc']))
    bufr_csv.df['j0'] = np.int32(np.floor(bufr_csv.df['xlc']))

    bufr_csv.df.drop(index=np.where(np.logical_or(bufr_csv.df['i0'] < 0, 
                                                  bufr_csv.df['i0'] > imax))[0], inplace=True)
    bufr_csv.df.reset_index(drop=True, inplace=True)
    bufr_csv.df.drop(index=np.where(np.logical_or(bufr_csv.df['j0'] < 0, 
                                                  bufr_csv.df['j0'] > jmax))[0], inplace=True)
    bufr_csv.df.reset_index(drop=True, inplace=True)

    # Compute interpolation weights
    bufr_csv.df['iwgt'] = 1. - (bufr_csv.df['ylc'] - bufr_csv.df['i0'])
    bufr_csv.df['jwgt'] = 1. - (bufr_csv.df['xlc'] - bufr_csv.df['j0'])

    if debug > 0:
        print('Finished with map projection and computing horiz interp weights (time = %.3f s)' % 
              (dt.datetime.now() - start_map_proj).total_seconds())
        print('# BUFR entries remaining = %d' % len(bufr_csv.df))

    # Create output DataFrame
    out_df = bufr_csv.df.copy()
    drop_idx = []

    # Determine row indices for each ob type
    ob_idx = {}
    obs = ['ADPSFC', 'SFCSHP', 'MSONET', 'GPSIPW', 'ADPUPA', 'AIRCAR', 'AIRCFT']
    for o in obs:
        ob_idx[o] = list(np.where(out_df['subset'] == o)[0]) 

    # Initialize variables for 3D obs as zeros
    for v in ['QOB', 'TOB', 'ZOB', 'UOB', 'VOB']:
        rows = np.intersect1d(np.array(ob_idx['ADPUPA'] + ob_idx['AIRCAR'] + ob_idx['AIRCFT']), 
                              np.where(np.logical_not(np.isnan(out_df[v]))))
        out_df.loc[rows, v] = 0.

    # Set all ZOB values to NaN if we don't wish to interpolate ZOBs for AIRCAR and AIRCFT
    if not interp_z_aircft:
        out_df.loc[ob_idx['AIRCAR'], 'ZOB'] = np.nan
        out_df.loc[ob_idx['AIRCFT'], 'ZOB'] = np.nan

    # Add some DataFrame columns
    nrow = len(out_df)
    extra_col_int = ['pi0']
    extra_col_float = ['twgt', 'pwgt']
    for c in extra_col_int:
        out_df[c] = np.zeros(nrow, dtype=int)
    for c in extra_col_float:
        out_df[c] = np.zeros(nrow, dtype=float)
    extra_col_int = extra_col_int + ['i0', 'j0']
    extra_col_float = extra_col_float + ['xlc', 'ylc', 'iwgt', 'jwgt']


    #-----------------------------------------------------------------------------------------------
    # Create Obs Based on 2-D Fields (ADPSFC, SFCSHP, MSONET, GPSIPW)
    #-----------------------------------------------------------------------------------------------

    # Extract 2D fields ONLY
    fields2D = ['PRES_P0_L1_GLC0', 'SPFH_P0_L103_GLC0', 'TMP_P0_L103_GLC0', 'HGT_P0_L1_GLC0', 
                'PWAT_P0_L200_GLC0']
    if interp_latlon:
        fields2D.append('gridlat_0')
        fields2D.append('gridlon_0')
    fields2D1 = ['UGRD_P0_L103_GLC0', 'VGRD_P0_L103_GLC0']
    wrf_data = {}
    for hr in wrf_hr:
        wrf_data[hr] = {}
        for f in fields2D:
            wrf_data[hr][f] = wrf_ds[hr][f][:, :].values
        for f in fields2D1:
            wrf_data[hr][f] = wrf_ds[hr][f][0, :, :].values

    # Loop over each ADPSFC, SFCSHP, MSONET, and GPSIPW observation
    for j in (ob_idx['ADPSFC'] + ob_idx['SFCSHP'] + ob_idx['MSONET'] + ob_idx['GPSIPW']):
        
        if debug > 1:
            time_jstart = dt.datetime.now()
            print()

        # Determine WRF hour right before observation and weight for temporal interpolation
        ihr, twgt = cou.determine_twgt(wrf_hr, out_df.loc[j, 'DHR'])

        # Determine surface height above sea level and surface pressure
        sfch = cou.interp_x_y_t(wrf_data, wrf_hr, 'HGT_P0_L1_GLC0', out_df.loc[j], ihr, twgt)
        sfcp = cou.interp_x_y_t(wrf_data, wrf_hr, 'PRES_P0_L1_GLC0', out_df.loc[j], ihr, twgt) * 1e-2

        if debug > 1:
            time_sfc = dt.datetime.now()
            print('2D field: nmsg = %d, subset = %s, TYP = %d' % 
                  (out_df.loc[j, 'nmsg'], out_df.loc[j, 'subset'], out_df.loc[j, 'TYP']))
            print('done determining twgt, sfch, and sfcp (%.6f s)' % 
                  (time_sfc - time_jstart).total_seconds())
            print('twgt = %.3f, i0 = %d, j0 = %d' % 
                  (twgt, out_df.loc[j, 'i0'], out_df.loc[j, 'j0']))
            print('DHR = %.3f, XOB = %.6f, YOB = %.6f' % 
                  (out_df.loc[j, 'DHR'], out_df.loc[j, 'XOB'], out_df.loc[j, 'YOB']))
            print('sfch = %.6f, sfcp = %.6f' % (sfch, sfcp))

        # Surface obs only: Reset surface values to match NR and assign surface pressure values
        if out_df.loc[j, 'subset'] in ['ADPSFC', 'SFCSHP', 'MSONET']:
            for o, v in zip(['ZOB', 'ELV', 'POB', 'PRSS'], [sfch, sfch, sfcp, sfcp]):
                if not np.isnan(out_df.loc[j, o]):
                    out_df.loc[j, o] = v

        # Interpolate temporally
        obs_name = ['QOB', 'TOB', 'UOB', 'VOB', 'PWO']
        wrf_name = ['SPFH_P0_L103_GLC0', 'TMP_P0_L103_GLC0', 'UGRD_P0_L103_GLC0', 
                    'VGRD_P0_L103_GLC0', 'PWAT_P0_L200_GLC0']
        if interp_latlon:
            obs_name = obs_name + ['XOB', 'YOB']
            wrf_name = wrf_name + ['gridlon_0', 'gridlat_0']
        for o, m in zip(obs_name, wrf_name):
            if not np.isnan(out_df.loc[j, o]):
                if debug > 1:
                    time_before_interp = dt.datetime.now()
                out_df.loc[j, o] = cou.interp_x_y_t(wrf_data, wrf_hr, m, out_df.loc[j], ihr, twgt)
                if debug > 1:
                    time_after_interp = dt.datetime.now()
                    print('finished interp for %s (%.6f s)' % 
                          (o, (time_after_interp - time_before_interp).total_seconds()))
                    entry_times2d.append((dt.datetime.now() - time_jstart).total_seconds())


    #-----------------------------------------------------------------------------------------------
    # Create Obs Based on 3-D Fields (ADPUPA, AIRCAR, AIRCFT)
    #-----------------------------------------------------------------------------------------------

    print()
    print('3D Observations')
    print('---------------')
    print()

    # Remove data from wrf_data that we don't need so we can free up some memory
    for hr in wrf_hr:
        for f in ['SPFH_P0_L103_GLC0', 'TMP_P0_L103_GLC0', 'PWAT_P0_L200_GLC0', 'UGRD_P0_L103_GLC0', 
                  'VGRD_P0_L103_GLC0']:
            del wrf_data[hr][f]

    # Create array to save p1d arrays in
    p1d = np.zeros([100, len(out_df)])
    pdone = np.zeros(len(out_df), dtype=int)

    # We will extract 3D fields one at a time b/c these 3D arrays are massive (~6.5 GB each), so it 
    # is not feasible to load all of them at once. Ideally, only 1 should be loaded at any given
    # time. Therefore, unlike the 2D obs, we will loop over each obs time interval and each WRF
    # field.

    # Another approach that could have been taken here is to lazily load each 3D field using Xarray
    # datasets. While this would address the memory issues, it would create computational issues b/c
    # the dataset will have to be queried for each BUFR entry and each variable, which is time 
    # consuming. Tests comparing the use of Xarray datasets to regular Numpy arrays show that the
    # Numpy array approach is ~23x faster.

    # Loop over each variable (starting with pressure)
    obs_name = ['POB', 'TOB', 'QOB', 'ZOB', 'UOB', 'VOB']
    wrf_name = ['PRES_P0_L105_GLC0', 'TMP_P0_L105_GLC0', 'SPFH_P0_L105_GLC0', 
                'HGT_P0_L105_GLC0', 'UGRD_P0_L105_GLC0', 'VGRD_P0_L105_GLC0']
    for o, f in zip(obs_name, wrf_name):

        # Loop over each WRF time
        for hr in wrf_hr: 

            # Determine indices of obs within wrf_step of this output time
            ind = np.where(np.logical_and(out_df['DHR'] > (hr - wrf_step_dec), 
                                          out_df['DHR'] < (hr + wrf_step_dec)))
            ind = np.intersect1d(ind, np.array(ob_idx['ADPUPA'] + ob_idx['AIRCAR'] + ob_idx['AIRCFT']))

            # If no indices, move to next time
            if ind.size == 0:
                continue

            print()
            print('3D Interp: %s' % f)
            print()

            # Extract field from UPP
            # It seems rather silly to include [:, :, :] before calling .values, but this really
            # helps with memory management. Including these indices allows the program to deallocate
            # wrf3d when setting wrf3d = 0.
            if debug > 0:
                time_3d = dt.datetime.now()
            wrf3d = wrf_ds[hr][f][:, :, :].values
            if debug > 0:
                print('time to extract %s = %.6f s' % (f, (dt.datetime.now() - time_3d).total_seconds()))
                for l in os.popen('free -t -m -h').readlines():
                    print(l)
                print()

            # Loop over each ADPUPA, AIRCAR, and AIRCFT observation within this time interval
            for j in ind:

                # Check to make sure that we didn't already drop this row
                if j in drop_idx:
                    continue

                # Drop row if POB is missing
                if np.isnan(out_df.loc[j, 'POB']):
                    drop_idx.append(j)
                    continue

                if debug > 1:
                    time_jstart = dt.datetime.now()
                    print()

                if o == 'POB': 
   
                    # First half of pressure calculation (P1)
                    if np.isclose(p1d[0, j], 0.):         

                        # Determine weight for temporal interpolation
                        ihr, twgt = cou.determine_twgt(wrf_hr, out_df.loc[j, 'DHR'])
                        out_df.loc[j, 'twgt']  = twgt 

                        # Determine surface height above sea level and surface pressure
                        if hr > out_df.loc[j, 'DHR']:
                            hr1 = hr - wrf_step_dec
                            hr2 = hr
                        else:
                            hr1 = hr
                            hr2 = hr + wrf_step_dec
                        sfch = cou.interp_x_y_t(wrf_data, wrf_hr, 'HGT_P0_L1_GLC0', out_df.loc[j],
                                                ihr, twgt)
                        sfcp = cou.interp_x_y_t(wrf_data, wrf_hr, 'PRES_P0_L1_GLC0', out_df.loc[j], 
                                                ihr, twgt) * 1e-2
            
                        # Check whether ob from BUFR file lies underground
                        if (out_df.loc[j, 'ZOB'] < sfch) or (out_df.loc[j, 'POB'] > sfcp):
                            drop_idx.append(j)
                            continue

                        if debug > 1:
                            time_sfc = dt.datetime.now()
                            print('done determining twgt, sfch, and sfcp (%.6f s)' % 
                                  (time_sfc - time_jstart).total_seconds())
            
                        # Perform vertical interpolation in pressure rather than height b/c some obs don't 
                        # have height. The calculation below only gives use part of the p1d array. The 
                        # rest of this of this array will be determined when we loop over the other 3D
                        # pressure array
                        p1d[:, j] = twgt * 1e-2 * cou.interp_x_y(wrf3d, out_df.loc[j], threeD=True)
 
                        if debug > 1:
                            print('total time for p1 = %.6f s' % 
                                  (dt.datetime.now() - time_jstart).total_seconds())
                            entry_timesp1.append((dt.datetime.now() - time_jstart).total_seconds())

                        # Special case: twgt = 1. In this case, we don't need to interpolate in time, so 
                        # we can skip the P2 section of the conditional
                        if np.isclose(twgt, 1):
                            pdone[j] = 1

                    # Second half of pressure calculation (P2)
                    else:
                        p1d[:, j] = p1d[:, j] + ((1.-out_df.loc[j, 'twgt']) * 1e-2 * 
                                                 cou.interp_x_y(wrf3d, out_df.loc[j], threeD=True))
                        pdone[j] = 1
                        
                        if debug > 1:
                            print('total time for p2 = %.6f s' % 
                                  (dt.datetime.now() - time_jstart).total_seconds())
                            entry_timesp2.append((dt.datetime.now() - time_jstart).total_seconds())

                    if pdone[j]: 
                        # Check for extrapolation
                        if (p1d[0, j] > out_df.loc[j, 'POB']):
                            out_df.loc[j, 'pi0'] = np.where(p1d[:, j] > out_df.loc[j, 'POB'])[0][-1]
                            if debug > 1:
                                print('pi0 = %d' % out_df.loc[j, 'pi0'])
                            if out_df.loc[j, 'pi0'] >= (p1d.shape[0] - 1):
                                # Prevent extrapolation in vertical
                                drop_idx.append(j)
                                continue
                            out_df.loc[j, 'POB'], out_df.loc[j, 'pwgt'] = cou.interp_wrf_p1d(p1d[:, j], out_df.loc[j])
                        else:
                            drop_idx.append(j)
                            continue

                        if debug > 1:
                            print('interpolated P = %.2f' % out_df.loc[j, 'POB'])
                            print('actual P = %.2f' % bufr_csv.df.loc[j, 'POB'])

                # Perform interpolation for other variables (not pressure)
                else:

                    if not np.isnan(out_df.loc[j, o]):
                        if debug > 1:
                            time_before_interp = dt.datetime.now()
                        twgt =  out_df.loc[j, 'twgt'] 
                        if np.isclose(out_df.loc[j, o], 0):
                            out_df.loc[j, o] = twgt * cou.interp_x_y_z(wrf3d, out_df.loc[j])
                        else:
                            out_df.loc[j, o] = ((1.-twgt) * cou.interp_x_y_z(wrf3d, out_df.loc[j]) + 
                                                out_df.loc[j, o])
                        if debug > 1:
                            time_after_interp = dt.datetime.now()
                            print('finished interp for %s (%.6f s)' % 
                                  (o, (time_after_interp - time_before_interp).total_seconds()))

                    if debug > 1:
                        print('total time = %.6f s' % (dt.datetime.now() - time_jstart).total_seconds())
                        entry_times3d.append((dt.datetime.now() - time_jstart).total_seconds())

            # Free up memory (shouldn't have to call garbage collector after this)
            wrf3d = 0.

    # Drop rows that we skipped as well as the extra columns we added
    debug_df = out_df.copy()
    out_df.drop(labels=extra_col_int, axis=1, inplace=True)
    out_df.drop(labels=extra_col_float, axis=1, inplace=True)
    out_df.drop(index=drop_idx, inplace=True)
    out_df.reset_index(drop=True, inplace=True)
    bufr_csv.df.drop(index=drop_idx, inplace=True)
    bufr_csv.df.reset_index(drop=True, inplace=True)
    debug_df.drop(index=drop_idx, inplace=True)
    debug_df.reset_index(drop=True, inplace=True)

    # Convert to proper units
    out_df['QOB'] = out_df['QOB'] * 1e6
    out_df['TOB'] = out_df['TOB'] - 273.15
    out_df['PWO'] = (out_df['PWO'] / 997.) * 1000.
    out_df['ELV'] = np.int64(out_df['ELV'])
    if interp_latlon:
        idx = np.where((out_df['subset'] == 'ADPSFC') | (out_df['subset'] == 'SFCSHP') |
                       (out_df['subset'] == 'MSONET'))[0] 
        out_df.loc[idx, 'XOB'] = out_df.loc[idx, 'XOB'] + 360.

    # If we didn't interpolate ZOB for AIRCAR and AIRCFT, copy values from original BUFR file
    if not interp_z_aircft:
        air_idx = np.where((out_df['subset'] == 'AIRCAR') | (out_df['subset'] == 'AIRCFT'))[0]
        out_df.loc[air_idx, 'ZOB'] = bufr_csv.df.loc[air_idx, 'ZOB']

    # Write output DataFrame to a CSV file
    # real_red.prepbufr.csv file can be used for assessing interpolation accuracy
    bufr.df_to_csv(out_df, fake_bufr_dir + t.strftime('/%Y%m%d%H%M.fake.prepbufr.csv'))
    bufr.df_to_csv(bufr_csv.df, fake_bufr_dir + t.strftime('/%Y%m%d%H%M.real_red.prepbufr.csv'))

    # Timing
    print()
    print('time for this BUFR file = %.6f s' % (dt.datetime.now() - start_loop).total_seconds())

if debug > 1:
    print()
    for a, s in zip([entry_times2d, entry_timesp1, entry_timesp2, entry_times3d],
                    ['2D', 'P1', 'P2', '3D']):
        if len(a) > 0:
            print('avg time per %s entries = %.6f s' % (s, np.mean(np.array(a))))

# Total timing
print()
print('time for entire program = %.6f s' % (dt.datetime.now() - begin).total_seconds())


"""
End create_conv_obs.py
"""
