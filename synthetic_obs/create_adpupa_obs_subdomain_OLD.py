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

"""
Development Notes:

1. Current version takes ~1800 s (30 min) to extract all the fields for a single radiosonde ob when 
    subdomain_half_size = 400. This is unacceptable b/c there are ~70 radiosonde sites in the US, 
    which means it will take ~35 hrs just to extract all the data for the radiosondes. Based on the 
    timing tests below, a better option might be to (1) only extract the p, u, and v arrays so that
    the radiosonde locations can be determined, then interpolate to these locations in create_conv_obs.py
    and (2) extract the entire p, u, and v arrays for the two times straddling the radiosonde obs, 
    then compute the radiosonde locations for all sites at once before moving to the next time window.
    Note: To make (2) tractable, the p, u, and v arrays will need to be thinned in the horizontal 
    (grabbing every 3rd value should be sufficient). Make this thinning value an input parameter so 
    it can be easily changed in the future if desired.

    a. Timing tests:
        i. Extracting a single value from a UPP dataset that's not loaded takes ~0.00117 s (using
            ds[i, j, k].values)
        ii. Extracting a vertical slice (100 values) from a UPP dataset that's not loaded takes
            ~27.7 s (using ds[:, j, k].values)
        iii. Extracting a vertical slice (10 values) from a UPP dataset that's not loaded takes
            ~2.96 s (using ds[k1:k2, j, k].values)
        iv. Extracting all fields for a single radiosonde ob when subdomain_half_size = 200 takes
            ~1700 s (so essentially the same as when the subdomain_half_size = 400). This agrees 
            with other tests which suggest that extracting a domain that's larger in the x and y 
            directions does not increase the time it takes to extract the data from the DataSet.

    b. Memory considerations:
        i. A 3D array takes up ~6.8 GB in memory 
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
wrf_dir = work + '/BMC/wrfruc/murdzek/nature_run_spring/UPP/20220429/'
#wrf_dir = work + '/BMC/wrfruc/murdzek/nature_run_tests/nature_run_spring_v2/output/202204291200/UPP/'

# Directory containing real prepbufr CSV output
#bufr_dir = work + '/BMC/wrfruc/murdzek/real_obs/obs_rap/'
#bufr_dir = work + '/BMC/wrfruc/murdzek/py_scripts/synthetic_obs/'
#bufr_dir = work + '/BMC/wrfruc/murdzek/src/py_scripts/synthetic_obs/'
bufr_dir = work + '/BMC/wrfruc/murdzek/real_obs/test/'

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

# Size of subdomain to extract around sounding launch site (in gridpoints)
subdomain_half_size = 400.

# Option for debugging output (0 = none, 1 = some, 2 = a lot)
debug = 2

# Save full dataframe for debugging?
save_debug_df = True


#---------------------------------------------------------------------------------------------------
# Create Simulated ADPUPA Observations
#---------------------------------------------------------------------------------------------------

# Timing
begin = dt.datetime.now()

# Constants
g = 9.81
Rd = 287.04

# List to save time spent on each BUFR entry
if debug > 1:
    entry_times = []

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
    bufr_csv.df.drop(index=np.where(bufr_csv.df['subset'] != 'ADPUPA')[0], inplace=True)
    bufr_csv.df.reset_index(drop=True, inplace=True)
    if debug > 0:
        print('total # of BUFR entries = %d' % len(bufr_csv.df))

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
        out_df[v] = 0.

    # Add some DataFrame columns
    nrow = len(out_df)
    extra_col_int = ['pi0']
    extra_col_float = ['twgt', 'pwgt', 'adpupa_x', 'adpupa_y', 'adpupa_t_s']
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

        single_df = out_df.loc[out_df['SID'] == sid].copy()
        single_df.sort_values('POB', ascending=False)
        adpupa_x, adpupa_y = mp.ll_to_xy_lc(single_df['YOB'].iloc[0], 
                                            single_df['XOB'].iloc[0] - 360.)
        adpupa_t_s = single_df['DHR'].iloc[0] * 3600.

        ihr = np.where(wrf_sec <= adpupa_t_s)[0][-1]
        wrf_t1_s = wrf_sec[ihr]
        wrf_t2_s = wrf_sec[ihr+1]

        # Extract WRF grids around the radiosonde launch location
        # Extracting these fields takes A LOT of time. Thus, we will only extract the fields needed
        # to propagate the radiosonde (p, u, v) and interpolate the other fields afterwards using the
        # regular create_conv_obs.py program using the updated radiosonde locations computed here
        if debug > 0:
            time_3d = dt.datetime.now()
        istart = int(max(int(adpupa_y) - subdomain_half_size, 0))
        iend = int(min(int(adpupa_y) + subdomain_half_size + 1, imax+2))
        jstart = int(max(int(adpupa_x) - subdomain_half_size, 0))
        jend = int(min(int(adpupa_x) + subdomain_half_size + 1, jmax+2))
        fields2D = ['PRES_P0_L1_GLC0', 'HGT_P0_L1_GLC0']
        fields3D = ['PRES_P0_L105_GLC0', 'UGRD_P0_L105_GLC0', 'VGRD_P0_L105_GLC0']
        wrf_data = {}
        for hr in wrf_hr:
            print('hr = %.2f' % hr)
            wrf_data[hr] = {}
            for f in fields2D:
                print('extracting %s' % f)
                wrf_data[hr][f] = wrf_ds[hr][f][istart:iend, jstart:jend].values
            for f in fields3D:
                print('extracting %s' % f)
                wrf_data[hr][f] = wrf_ds[hr][f][:, istart:iend, jstart:jend].values
        adpupa_x = adpupa_x - jstart
        adpupa_y = adpupa_y - istart
        single_df['xlc'] = single_df['xlc'] - jstart
        single_df['ylc'] = single_df['ylc'] - istart
        single_df['i0'] = single_df['i0'] - istart
        single_df['j0'] = single_df['j0'] - jstart

        if debug > 0:
            time_extract = dt.datetime.now()
            print('time to extract fields = %.6f s' % (time_extract - time_3d).total_seconds())
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
            print('done determining sfch and sfcp (%.6f s)' % (time_sfc - time_extract).total_seconds())

        # Now, we're ready to compute the radiosonde obs, starting from the lowest pressure level
        sorted_idx = single_df.index
        for n, k in enumerate(sorted_idx):

            if debug > 1:
                kloop_start = dt.datetime.now()

            # Save radiosonde location
            if save_debug_df:
                out_df.loc[k, 'adpupa_x'] = adpupa_x
                out_df.loc[k, 'adpupa_y'] = adpupa_y
                out_df.loc[k, 'adpupa_t_s'] = adpupa_t_s

            # Determine weights for interpolation
            out_df.loc[k, 'i0'] = np.int32(np.floor(adpupa_y)) 
            out_df.loc[k, 'j0'] = np.int32(np.floor(adpupa_x))
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
                    print('pi0 = %d' % out_df.loc[k, 'pi0'])
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
                print('time to determine interpolation wgts = %.6f s' % (time_wgts - kloop_start).total_seconds())

            # Interpolate to observation location
            for ob, f in zip(['UOB', 'VOB'], fields3D[1:]):
                if not np.isnan(out_df.loc[k, ob]):
                    out_df.loc[k, ob] = (twgt * cou.interp_wrf_3d(wrf_data[wrf_hr[ihr]][f], out_df.loc[k]) +
                                         (1.-twgt) * cou.interp_wrf_3d(wrf_data[wrf_hr[ihr+1]][f], out_df.loc[k]))
            if debug > 1:
                time_interp = dt.datetime.now()
                print('finished interp for U and V (%.6f s)' % (time_interp - time_wgts).total_seconds())
            
            # Update (x, y) coordinate for next observation
            if k != sorted_idx[-1]:
                delta_p = single_df['POB'].iloc[n+1] - out_df[k, 'POB']
                if np.isclose(delta_p, 0):
                    continue
                omega = -1e-2 * g * ascent_spd * out_df.loc[k, 'POB'] / (Rd * out_df.loc[k, 'TOB'])
                delta_t = delta_p / omega
                adpupa_x = adpupa_x + (1e3 * out_df.loc[k, 'UOB'] * delta_t)
                adpupa_y = adpupa_y + (1e3 * out_df.loc[k, 'VOB'] * delta_t)
                adpupa_t_s = adpupa_t_s + delta_t
     
                # Switch to next radiosonde if balloon exits the subdomain
                if ((adpupa_x < 0) or (adpupa_y < 0) or (adpupa_x > (2*subdomain_half_size)) or
                    (adpupa_y > (2*subdomain_half_size))):
                    drop_idx = drop_idx + list(sorted_idx[n+1:])
                    break
 
                # Update ihr if needed
                if adpupa_t_s > wrf_t2_s:
                    ihr = ihr + 1
                    wrf_t1_s = wrf_sec[ihr] 
                    wrf_t2_s = wrf_sec[ihr+1]
 
            if debug > 1:
                time_advect = dt.datetime.now()
                entry_times.append(time_advect - kloop_start)
                print('finished advecting radiosonde (%.6f s)' % (time_advect - time_interp).total_seconds())

    # Drop rows that we skipped
    out_df.drop(index=drop_idx, inplace=True)
    out_df.reset_index(drop=True, inplace=True)
    bufr_csv.df.drop(index=drop_idx, inplace=True)
    bufr_csv.df.reset_index(drop=True, inplace=True)
    
    # Convert to proper units
    out_df['QOB'] = out_df['QOB'] * 1e6
    out_df['TOB'] = out_df['TOB'] - 273.15
    out_df['ELV'] = np.int64(out_df['ELV'])

    if save_debug_df:
        bufr.df_to_csv(out_df, fake_bufr_dir + t.strftime('/%Y%m%d%H%M.debug.adpupa.csv'))
    
    # Drop the extra columns we added
    out_df.drop(labels=extra_col_int, axis=1, inplace=True)
    out_df.drop(labels=extra_col_float, axis=1, inplace=True)

    # Write output DataFrame to a CSV file
    # real_red.prepbufr.csv file can be used for assessing interpolation accuracy
    bufr.df_to_csv(out_df, fake_bufr_dir + t.strftime('/%Y%m%d%H%M.fake.adpupa.csv'))
    bufr.df_to_csv(bufr_csv.df, fake_bufr_dir + t.strftime('/%Y%m%d%H%M.real_red.adpupa.csv'))

    # Timing
    print()
    print('time for this BUFR file = %.6f s' % (dt.datetime.now() - start_loop).total_seconds())

if debug > 1:
    print()
    if len(entry_times) > 0:
        print('avg time per entry = %.6f s' % (s, np.mean(np.array(entry_times))))

# Total timing
print()
print('time for entire program = %.6f s' % (dt.datetime.now() - begin).total_seconds())


"""
End create_adpupa_obs.py
"""
