"""
Create Synthetic ADPUPA Observations from WRF Nature Run Output

ADPUPA obs are created by assuming that the radiosonde balloons drift with the horizontal wind
and have a constant ascent rate.

Prepbufr CSV files for real observations can be created using GSI-utils/bin/prepbufr_decode_csv.x.
To use this utility, do the following:

1. Copy the prepbufr file you wish to decode to GSI-utils/bin and rename it 'prepbufr'
2. Run ./prepbufr_decode_csv.x
3. Move resulting CSV file to wherever you choose

This code uses UPP output on WRF native levels (wrfnat). This is likely more efficient than using 
the wrfred files, which can take some time (and a lot of memory) to create.

Debugging Notes:

1. In IPython, the command `%whos` can be used to view all variables in the namespace, which can be
    used to determine which variables are sucking down a lot of RAM.

Resources:

Hera - An interactive session on a regular node is sufficient
Jet - 4 GB of memory on tjet is sufficient if only interpolating 2D fields
    - 25 GB of memory on sjet is sufficient if interpolating 3D fields

shawn.s.murdzek@noaa.gov
Date Created: 27 February 2023
Environment: adb_graphics (Jet) or pygraf (Hera)
"""


#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import numpy as np
import datetime as dt
import math
import os

import create_ob_utils as cou 
import bufr
import map_proj as mp


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

work = '/mnt/lfs4'

# Directory containing wrfnat output from UPP
wrf_dir = work + '/BMC/wrfruc/murdzek/nature_run_spring/UPP'
#wrf_dir = work + '/BMC/wrfruc/murdzek/nature_run_tests/nature_run_spring_v2/output/202204291200/UPP/'

# Directory containing real prepbufr CSV output
bufr_dir = work + '/BMC/wrfruc/murdzek/sample_real_obs/obs_rap/'
#bufr_dir = work + '/BMC/wrfruc/murdzek/py_scripts/synthetic_obs/'
#bufr_dir = work + '/BMC/wrfruc/murdzek/src/py_scripts/synthetic_obs/'
#bufr_dir = work + '/BMC/wrfruc/murdzek/sample_real_obs/test/'

# Output directory for synthetic prepbufr CSV output
fake_bufr_dir = work + '/BMC/wrfruc/murdzek/nature_run_spring/synthetic_obs/'
#fake_bufr_dir = work + '/BMC/wrfruc/murdzek/nature_run_tests/nature_run_spring_v2/synthetic_obs/'

# Start and end times for prepbufrs. Step is in min
bufr_start = dt.datetime(2022, 4, 29, 12)
bufr_end = dt.datetime(2022, 4, 29, 13)
bufr_step = 120

# Start and end times for wrfnat UPP output. Step is in min
wrf_start = dt.datetime(2022, 4, 29, 12, 0)
wrf_end = dt.datetime(2022, 4, 29, 15, 0)
wrf_step = 15

# Ascent rate for radiosonde. 5 m/s value comes from this NWS report: 
# https://www.weather.gov/media/upperair/Documents/Radiosonde%20Ascent%20Rates.pdf
ascent_spd = 5.

# Allow vertical velocity from nature run to influence radiosonde ascent rate?
# Setting this to try might lead to some wonky results if the radiosonde enters a downdraft.
w_impact = True

# Size of subdomain to extract around sounding launch site (in gridpoints)
subdomain_half_size = 400.

# Option for debugging output (0 = none, 1 = some, 2 = a lot)
debug = 2


#---------------------------------------------------------------------------------------------------
# Create Simulated ADPUPA Observations
#---------------------------------------------------------------------------------------------------

# Timing
begin = dt.datetime.now()

# List to save time spent on each BUFR entry
if debug > 1:
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
    if debug > 0:
        print('total # of BUFR entries = %d' % len(bufr_csv.df))
        print('available BUFR obs =', obs)
    bufr_csv.df.drop(index=np.where(bufr_csv.df['subset'] != 'ADPUPA')[0], inplace=True)
    bufr_csv.df.reset_index(drop=True, inplace=True)

    # Remove obs outside of the range of wrfnat files
    hr_min = ((wrf_start - t).days * 24) + ((wrf_start - t).seconds / 3600.)
    hr_max = ((wrf_end - t).days * 24) + ((wrf_end - t).seconds / 3600.)
    bufr_csv.df.drop(index=np.where(np.logical_or(bufr_csv.df['DHR'] < hr_min, 
                                                  bufr_csv.df['DHR'] > hr_max))[0], inplace=True)
    bufr_csv.df.reset_index(drop=True, inplace=True)

    # Estimate how long it will take the last radiosonde to reach the model top (assumed to be
    # ~27,000 m, which is around 20 hPa)
    last_dhr = bufr_csv.df['DHR'].max()
    max_wrfnat_dhr_needed = last_dhr + ((27000. / ascent_spd) / 3600.)

    # Open wrfnat files
    hr_start = math.floor(bufr_csv.df['DHR'].min()*4) / 4
    hr_end = math.ceil(max_wrfnat_dhr_needed*4) / 4 + (2*wrf_step_dec)
    print('min/max WRF hours = %.2f, %.2f' % (hr_min, hr_max))
    print('min/max BUFR hours = %.2f, %.2f' % (bufr_csv.df['DHR'].min(), bufr_csv.df['DHR'].max()))
    wrf_ds = {}
    wrf_hr = np.arange(hr_start, hr_end, wrf_step_dec)
    for hr in wrf_hr:
        wrf_t = t + dt.timedelta(hours=hr)
        print(wrf_dir + wrf_t.strftime('wrfnat_%Y%m%d%H%M.grib2'))
        wrf_ds[hr] = xr.open_dataset(wrf_dir + wrf_t.strftime('wrfnat_%Y%m%d%H%M.grib2'), 
                                     engine='pynio')
    wrf_sec = wrf_hr * 3600.

    print('time to open GRIB files = %.2f s' % (dt.datetime.now() - start_loop).total_seconds())
    
    # Remove obs outside of the spatial domain of the wrfnat files (first, inprecise pass)
    latmin = wrf_ds[0]['gridlat_0'].min().values
    latmax = wrf_ds[0]['gridlat_0'].max().values
    lonmin = wrf_ds[0]['gridlon_0'].min().values + 360.
    lonmax = wrf_ds[0]['gridlon_0'].max().values + 360.
    bufr_csv.df.drop(index=np.where(np.logical_or(bufr_csv.df['XOB'] < lonmin, 
                                                  bufr_csv.df['XOB'] > lonmax))[0], inplace=True)
    bufr_csv.df.reset_index(drop=True, inplace=True)
    bufr_csv.df.drop(index=np.where(np.logical_or(bufr_csv.df['YOB'] < latmin, 
                                                  bufr_csv.df['YOB'] > latmax))[0], inplace=True)
    bufr_csv.df.reset_index(drop=True, inplace=True)

    # Extract size of latitude and longitude grids
    shape = wrf_ds[0]['gridlat_0'].shape
    imax = shape[0] - 2
    jmax = shape[1] - 2
    
    # Compute (x, y) coordinates of obs using a Lambert Conformal projection
    # Remove obs outside of the wrfnat domain (second, precise pass)
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

    bufr_csv.df['iwgt'] = 1. - (bufr_csv.df['ylc'] - bufr_csv.df['i0'])
    bufr_csv.df['jwgt'] = 1. - (bufr_csv.df['xlc'] - bufr_csv.df['j0'])

    if debug > 0:
        print('Finished performing map projection and computing horiz interp weights (time = %.3f s)' % 
              (dt.datetime.now() - start_map_proj).total_seconds())
        print('# BUFR entries remaining = %d' % len(bufr_csv.df))

    # Create output DataFrame
    out_df = bufr_csv.df.copy()
    drop_idx = []

    # Initialize variables for 3D obs as zeros
    for v in ['QOB', 'TOB', 'ZOB', 'UOB', 'VOB']:
        out_df.loc[v] = 0.

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
    # Create Obs Based on 3-D Fields
    #-----------------------------------------------------------------------------------------------

    print()
    print('Compute ADPUPA Observations')
    print('-----------------------------')
    print()

    all_sid = np.unique(out_df['SID'])
    print('number of radiosondes = %d' % len(all_sid))

    for sid in all_sid: 

        single_df = out_df['SID'].loc[out_df['SID'] == sid].copy()
        single_df.sort_values('POB', ascending=False)
        adpupa_x, adpupa_y = mp.ll_to_xy_lc(single_df['YOB'].iloc[0], 
                                            single_df['XOB'].iloc[0] - 360.)
        adpupa_t_s = single_df['DHR'].iloc[0] * 3600.

        ihr = np.where(wrf_sec <= adpupa_t_s)[0][-1]
        wrf_t1_s = wrf_sec[ihr]
        wrf_t2_s = wrf_sec[ihr+1]

        # Extract WRF grids around the radiosonde launch location
        if debug > 0:
            time_3d = dt.datetime.now()
        istart = max(int(adpupa_y) - subdomain_half_size, 0)
        iend = min(int(adpupa_y) + subdomain_half_size + 1, imax+2)
        jstart = max(int(adpupa_x) - subdomain_half_size, 0)
        jend = min(int(adpupa_x) + subdomain_half_size + 1, jmax+2)
        fields2D = ['PRES_P0_L1_GLC0', 'HGT_P0_L1_GLC0']
        fields3D = ['PRES_P0_L105_GLC0', 'TMP_P0_L105_GLC0', 'SPFH_P0_L105_GLC0',
                    'HGT_P0_L105_GLC0', 'UGRD_P0_L105_GLC0', 'VGRD_P0_L105_GLC0']
        wrf_data = {}
        for hr in wrf_hr:
            wrf_data[hr] = {}
            for f in fields2D:
                wrf_data[hr][f] = wrf_ds[hr][f][istart:iend, jstart:jend].values
            for f in fields3D:
                wrf_data[hr][f] = wrf_ds[hr][f][:, istart:iend, jstart:jend].values
        adpupa_x = adpupa_x - jstart
        adpupa_y = adpupa_y - istart

        if debug > 0:
            print('time to extract %s = %.6f s' % (f, (dt.datetime.now() - time_3d).total_seconds()))
            for l in os.popen('free -t -m -h').readlines():
                print(l) 
            print()

        # Find first pressure level that lies above the surface
        ihr, twgt = cou.determine_twgt(wrf_hr, adpupa_t_s / 3600.)
        sfch = cou.interp_wrf_to_obs(wrf_data, wrf_hr, 'HGT_P0_L1_GLC0', single_df.iloc[0], ihr, 
                                     twgt)
        sfcp = cou.interp_wrf_to_obs(wrf_data, wrf_hr, 'PRES_P0_L1_GLC0', single_df.iloc[0], ihr, 
                                     twgt) * 1e-2
        indices = single_df.index.copy()
        for k in indices:
            if (out_df.loc[k, 'ZOB'] < sfch) or (out_df.loc[k, 'POB'] > sfcp):
                drop_idx.append(k)
                single_df.drop(k)
                continue
            else:
                break

        if debug > 1:
            time_sfc = dt.datetime.now()
            print('done determining sfch and sfcp (%.6f s)' % (time_sfc - time_3d).total_seconds())

        # Now, we're ready to compute the radiosonde obs, starting from the lowest pressure level
        for k in single_df.index:

            # Determine weights for interpolation
            out_df.loc[k, 'i0'] np.int32(np.floor(adpupa_y)) 
            out_df.loc[k, 'j0'] np.int32(np.floor(adpupa_x))
            out_df.loc[k, 'iwgt'] = 1. - (adpupa_y - out_df.loc[k, 'i0']) 
            out_df.loc[k, 'jwgt'] = 1. - (adpupa_x - out_df.loc[k, 'j0']) 
            
            prs_f = 'PRES_P0_L105_GLC0'
            p1d = 1e-2 * (twgt * cou._bilinear_interp_horiz(wrf_data[wrf_hr[ihr]][prs_f], 
                                                            out_df.loc[k, 'iwgt'], 
                                                            out_df.loc[k, 'jwgt'],
                                                            out_df.loc[k, 'i0'], 
                                                            out_df.loc[k, 'j0'], threeD=True) +
                          (1.-twgt) * cou._bilinear_interp_horiz(wrf_data[wrf_hr[ihr+1]][prs_f], 
                                                                 out_df.loc[k, 'iwgt'], 
                                                                 out_df.loc[k, 'jwgt'],
                                                                 out_df.loc[k, 'i0'], 
                                                                 out_df.loc[k, 'j0'], threeD=True))
            
            # Check for extrapolation
            if (p1d[0] > out_df.loc[k, 'POB']):
                out_df.loc[k, 'pi0'] = np.where(p1d > out_df.loc[k, 'POB'])[0][-1]
                if debug > 1:
                    print('pi0 = %d' % out_df.loc[j, 'pi0'])
                if out_df.loc[k, 'pi0'] >= (p1d.size - 1):
                    # Prevent extrapolation in vertical
                    drop_idx.append(k)
                    continue
                out_df.loc[k, 'POB'], out_df.loc[k, 'pwgt'] = cou.interp_wrf_p1d(p1d, out_df.loc[k])
            else:
                drop_idx.append(k)
                continue

            if debug > 1:
                time_wgts = dt.datetime.now()
                print('interpolated P = %.2f' % out_df.loc[k, 'POB'])
                print('actual P = %.2f' % bufr_csv.df.loc[k, 'POB'])
                print('time to determine interpolation wgts = %.6f s' % (time_wgts - time_sfc).total_seconds())

            # Interpolate to observation location
            for ob, f in zip(['TOB', 'QOB', 'ZOB', 'UOB', 'VOB'], fields3d[1:]):
                if not np.isnan(out_df.loc[k, o]):
                    out_df.loc[k, o] = (twgt * cou.interp_wrf_3d(wrf_data[wrf_hr[ihr]][f], out_df.loc[k]) +
                                        (1.-twgt) * cou.interp_wrf_3d(wrf_data[wrf_hr[ihr+1]][f], out_df.loc[k]))
            if debug > 1:
                time_interp = dt.datetime.now()
                print('finished interp for %s (%.6f s)' % (o, (time_interp - time_wgts).total_seconds()))
            
            # Update (x, y) coordinate for next observation
            """
            Start here next time!
                1. Determine distance (in vertical) to next pressure level. If distance = 0, don't change anything
                2. Determine time it will take to reach next pressure level using the radiosonde ascent rate
                3. Use (u, v) winds at current level to advect radiosonde (can get a dx and dy doing this)
                4. Update adpupa_x, adpupa_y, and adpupa_t_s accordingly
            """

            # Update ihr if needed
            if adpupa_t_s > wrf_t2_s:
                ihr = ihr + 1
                wrf_t1_s = wrf_sec[ihr] 
                wrf_t2_s = wrf_sec[ihr+1] 

    # Drop rows that we skipped as well as the extra columns we added
    out_df.drop(labels=extra_col_int, axis=1, inplace=True)
    out_df.drop(labels=extra_col_float, axis=1, inplace=True)
    out_df.drop(index=drop_idx, inplace=True)
    out_df.reset_index(drop=True, inplace=True)
    bufr_csv.df.drop(index=drop_idx, inplace=True)
    bufr_csv.df.reset_index(drop=True, inplace=True)

    # Convert to proper units
    out_df['QOB'] = out_df['QOB'] * 1e6
    out_df['TOB'] = out_df['TOB'] - 273.15
    out_df['PWO'] = (out_df['PWO'] / 997.) * 1000.
    out_df['ELV'] = np.int64(out_df['ELV'])
    if interp_latlon:
        idx = np.where((out_df['subset'] == 'ADPSFC') | (out_df['subset'] == 'SFCSHP') |
                       (out_df['subset'] == 'MSONET'))[0] 
        out_df.loc[idx, 'XOB'] = out_df.loc[idx, 'XOB'] + 360.

    # Write output DataFrame to a CSV file
    # real_red.prepbufr.csv file can be used for assessing interpolation accuracy
    bufr.df_to_csv(out_df, fake_bufr_dir + t.strftime('/%Y%m%d%H%M.fake.adpupa.csv'))
    bufr.df_to_csv(bufr_csv.df, fake_bufr_dir + t.strftime('/%Y%m%d%H%M.real_red.adpupa.csv'))

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
End create_adpupa_obs.py
"""
