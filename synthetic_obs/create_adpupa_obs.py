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

Resources: Most of the tests I have done with this program use 25-30 GB of RAM

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
wrf_dir = work + '/BMC/wrfruc/murdzek/nature_run_spring/UPP/20220429/'
#wrf_dir = work + '/BMC/wrfruc/murdzek/nature_run_tests/nature_run_spring_v2/output/202204291200/UPP/'

# Directory containing real prepbufr CSV output
bufr_dir = work + '/BMC/wrfruc/murdzek/real_obs/obs_rap/'
#bufr_dir = work + '/BMC/wrfruc/murdzek/py_scripts/synthetic_obs/'
#bufr_dir = work + '/BMC/wrfruc/murdzek/src/py_scripts/synthetic_obs/'
#bufr_dir = work + '/BMC/wrfruc/murdzek/real_obs/test/'

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

# When extracting the p, u, and v arrays to propagate the radiosonde, only retain every nth
# gridpoint, where n = thin. Smaller values consume more memory.
thin = 3

# Spacing between gridpoints after thinning (for the 1-km NR, wrf_dx = thin)
wrf_dx = thin

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

    # Open first wrfnat files
    wrf_hr = [math.floor(bufr_csv.df['DHR'].min()*4) / 4, 0]
    print('min/max WRF hours = %.2f, %.2f' % (hr_min, hr_max))
    print('min/max BUFR hours = %.2f, %.2f' % (bufr_csv.df['DHR'].min(), bufr_csv.df['DHR'].max()))
    wrf_t = [t + dt.timedelta(hours=wrf_hr[0]), 0]
    print(wrf_dir + wrf_t[0].strftime('wrfnat_%Y%m%d%H%M.grib2'))
    wrf_ds = [xr.open_dataset(wrf_dir + wrf_t[0].strftime('wrfnat_%Y%m%d%H%M.grib2'), 
                              engine='pynio'), 0]

    print('time to open first GRIB file = %.2f s' % (dt.datetime.now() - start_loop).total_seconds())
    
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
    imax = (shape[0] - 2) / wrf_dx
    jmax = (shape[1] - 2) / wrf_dx
    
    # Compute (x, y) coordinates of obs using a Lambert Conformal projection
    # Remove obs outside of the wrfnat domain (second, precise pass)
    if debug > 0:
        start_map_proj = dt.datetime.now()
        print('Performing map projection with obs...')

    bufr_csv.df['xlc'], bufr_csv.df['ylc'] = mp.ll_to_xy_lc(bufr_csv.df['YOB'], bufr_csv.df['XOB'] - 360.)
    bufr_csv.df['xlc'] = bufr_csv.df['xlc'] / wrf_dx
    bufr_csv.df['ylc'] = bufr_csv.df['ylc'] / wrf_dx
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

    # Sort DataFrame by descending pressure. Don't drop the original indices, b/c we'll sort the
    # DataFrame based on these indices later
    bufr_csv.df.sort_values(['SID', 'POB'], ascending=False, inplace=True)
    bufr_csv.df.reset_index(inplace=True)

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
    done_sid = []
    nsonde = len(all_sid)
    print('number of radiosondes = %d' % nsonde)

    # Create arrays of radiosonde locations, times, indices, and surface values. These will be used 
    # in the while loop below. Note that the units for the (x, y) coordinates are in gridpoints, NOT
    # km. To get the coordinates in km, multiply by wrf_dx.
    adpupa_x = np.zeros(len(all_sid))
    adpupa_y = np.zeros(len(all_sid))
    adpupa_t_s = np.zeros(len(all_sid))
    adpupa_idx = np.zeros(len(all_sid), dtype=int)
    adpupa_last_idx = np.zeros(len(all_sid), dtype=int)
    adpupa_sfch = np.zeros(len(all_sid))
    adpupa_sfcp = np.zeros(len(all_sid))
    for j, sid in enumerate(all_sid):
        idx_sid = out_df.loc[out_df['SID'] == sid].index[0]
        adpupa_x[j] = out_df['xlc'].iloc[idx_sid]
        adpupa_y[j] = out_df['ylc'].iloc[idx_sid]
        adpupa_t_s[j] = out_df['DHR'].iloc[idx_sid] * 3600.
        adpupa_idx[j] = idx_sid
        adpupa_last_idx[j] = out_df.loc[out_df['SID'] == sid].index[-1]

    # Start loop over each time window between wrfnat times
    fields2D = ['PRES_P0_L1_GLC0', 'HGT_P0_L1_GLC0']
    fields3D = ['PRES_P0_L105_GLC0', 'TMP_P0_L105_GLC0', 'UGRD_P0_L105_GLC0', 'VGRD_P0_L105_GLC0',
                'SPFH_P0_L105_GLC0', 'HGT_P0_L105_GLC0']
    ob_fields = ['POB', 'TOB', 'UOB', 'VOB', 'QOB', 'ZOB']
    itime = 0
    wrf_data = {}
    while len(done_sid) < nsonde:
        time_start_itime = dt.datetime.now()
        print()

        # Update the first wrfnat file to be the second wrfnat file from the previous iteration
        if itime > 0:
            del wrf_data[wrf_hr[0]]
            wrf_hr[0] = wrf_hr[1]
            wrf_t[0] = wrf_t[1]
            wrf_ds[0] = wrf_ds[1]
            print('top of itime loop (itime = %d, wrf_hr[0] = %.2f)' % (itime, wrf_hr[0]))
        else:
            print('top of itime loop (itime = %d, wrf_hr[0] = %.2f)' % (itime, wrf_hr[0]))
            wrf_data[wrf_hr[0]] = {}
            for f in fields2D:
                print('extracting %s' % f)
                wrf_data[wrf_hr[0]][f] = wrf_ds[0][f][::thin, ::thin].values
            for f in fields3D:
                print('extracting %s' % f)
                wrf_data[wrf_hr[0]][f] = wrf_ds[0][f][:, ::thin, ::thin].values
            xmax = wrf_data[wrf_hr[0]][fields2D[0]].shape[1]
            ymax = wrf_data[wrf_hr[0]][fields2D[0]].shape[0]

        # Extract data for the next wrfnat file
        wrf_hr[1] = wrf_hr[0] + (wrf_step / 60.)
        wrf_t[1] = t + dt.timedelta(hours=wrf_hr[1])
        wrf_ds[1] = xr.open_dataset(wrf_dir + wrf_t[1].strftime('wrfnat_%Y%m%d%H%M.grib2'), 
                                    engine='pynio')
        wrf_data[wrf_hr[1]] = {}
        for f in fields2D:
            print('extracting %s' % f)
            wrf_data[wrf_hr[1]][f] = wrf_ds[1][f][::thin, ::thin].values
        for f in fields3D: 
            print('extracting %s' % f)
            wrf_data[wrf_hr[1]][f] = wrf_ds[1][f][:, ::thin, ::thin].values
        
        if debug > 0:
            time_extract = dt.datetime.now()
            print('time to extract fields = %.6f s' % (time_extract - time_start_itime).total_seconds())
            for l in os.popen('free -t -m -h').readlines():
                print(l) 

        # Loop over each sounding station, skipping those that are done
        for idx_sid, sid in enumerate(all_sid): 
            if sid in done_sid:
                continue
            if adpupa_t_s[idx_sid] > (wrf_hr[1] * 3600.):
                continue
            print()
            print('SID = %s' % sid)

            # Find first pressure level that lies above the surface if not done so already
            ihr, twgt = cou.determine_twgt(wrf_hr, adpupa_t_s[idx_sid] / 3600.)
            if np.isclose(adpupa_sfch[idx_sid], 0):
                adpupa_sfch[idx_sid] = cou.interp_x_y_t(wrf_data, wrf_hr, 'HGT_P0_L1_GLC0', 
                                                        out_df.iloc[adpupa_idx[idx_sid]], ihr, 
                                                        twgt)
                adpupa_sfcp[idx_sid] = cou.interp_x_y_t(wrf_data, wrf_hr, 'PRES_P0_L1_GLC0', 
                                                        out_df.iloc[adpupa_idx[idx_sid]], ihr, 
                                                        twgt) * 1e-2
            
                indices = out_df.loc[out_df['SID'] == sid].index.values
                out_df.loc[indices, 'ELV'] = adpupa_sfch[idx_sid]
                for k in indices:
                    if (out_df.loc[k, 'POB'] > adpupa_sfcp[idx_sid]):
                        drop_idx.append(k)
                    else:
                        adpupa_idx[idx_sid] = k
                        break

                if debug > 1:
                    time_sfc = dt.datetime.now()
                    print('done determining sfch and sfcp (%.2f, %.2f, %.6f s)' % 
                          (adpupa_sfch[idx_sid], adpupa_sfcp[idx_sid], 
                           (time_sfc - time_extract).total_seconds()))
    
            # Now, we're ready to compute the radiosonde obs, starting from the lowest pressure level
            k = adpupa_idx[idx_sid]
            while (k <= adpupa_last_idx[idx_sid] and adpupa_t_s[idx_sid] <= (wrf_hr[1] * 3600.)):

                if debug > 1:
                    kloop_start = dt.datetime.now()

                # Save radiosonde location
                if save_debug_df:
                    out_df.loc[k, 'adpupa_x'] = adpupa_x[idx_sid]
                    out_df.loc[k, 'adpupa_y'] = adpupa_y[idx_sid]
                    out_df.loc[k, 'adpupa_t_s'] = adpupa_t_s[idx_sid]

                # Determine weights for interpolation
                out_df.loc[k, 'i0'] = np.int32(np.floor(adpupa_y[idx_sid])) 
                out_df.loc[k, 'j0'] = np.int32(np.floor(adpupa_x[idx_sid]))
                out_df.loc[k, 'iwgt'] = 1. - (adpupa_y[idx_sid] - out_df.loc[k, 'i0']) 
                out_df.loc[k, 'jwgt'] = 1. - (adpupa_x[idx_sid] - out_df.loc[k, 'j0']) 
            
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
                        # Balloon has exited the top of the domain
                        drop_idx = drop_idx + list(range(k, adpupa_last_idx[idx_sid]+1))
                        done_sid.append(sid)
                        print('RAOB %s exited top of domain' % sid)
                        break
                    out_df.loc[k, 'POB'], out_df.loc[k, 'pwgt'] = cou.interp_wrf_p1d(p1d, out_df.loc[k])
                else:
                    drop_idx.append(k)
                    k = k + 1
                    continue

                if debug > 1:
                    time_wgts = dt.datetime.now()
                    print('interpolated P = %.2f' % out_df.loc[k, 'POB'])
                    print('actual P = %.2f' % bufr_csv.df.loc[k, 'POB'])
                    print('time to determine interpolation wgts = %.6f s' % (time_wgts - kloop_start).total_seconds())

                # Interpolate to observation location
                tmp_sonde_vals = {}
                for ob, f in zip(ob_fields[1:], fields3D[1:]):
                    tmp_sonde_vals[ob] = (twgt * cou.interp_x_y_z(wrf_data[wrf_hr[ihr]][f], out_df.loc[k]) +
                                          (1.-twgt) * cou.interp_x_y_z(wrf_data[wrf_hr[ihr+1]][f], out_df.loc[k]))
                    if not np.isnan(out_df.loc[k, ob]):
                        out_df.loc[k, ob] = tmp_sonde_vals[ob]

                if debug > 1:
                    time_interp = dt.datetime.now()
                    print('finished interp for T, U, and V (%.6f s)' % (time_interp - time_wgts).total_seconds())
            
                # Update (x, y) coordinate for next observation
                if k == adpupa_last_idx[idx_sid]:
                    done_sid.append(sid)
                    break
                else:
                    delta_p = out_df.loc[k+1, 'POB'] - out_df.loc[k, 'POB']
                    if np.isclose(delta_p, 0):
                        k = k + 1
                        continue
                    omega = -g * ascent_spd * out_df.loc[k, 'POB'] / (Rd * tmp_sonde_vals['TOB'])
                    delta_t = delta_p / omega
                    adpupa_x[idx_sid] = adpupa_x[idx_sid] + (1e-3 * tmp_sonde_vals['UOB'] * delta_t / wrf_dx)
                    adpupa_y[idx_sid] = adpupa_y[idx_sid] + (1e-3 * tmp_sonde_vals['VOB'] * delta_t / wrf_dx)
                    adpupa_t_s[idx_sid] = adpupa_t_s[idx_sid] + delta_t    
 
                    # Switch to next radiosonde if balloon exits the domain
                    if ((adpupa_x[idx_sid] < 0) or (adpupa_x[idx_sid] > xmax) or 
                        (adpupa_x[idx_sid] < 0) or (adpupa_y[idx_sid] > ymax)):
                        drop_idx = drop_idx + list(range(k+1, adpupa_last_idx[idx_sid]+1))
                        done_sid.append(sid)
                        print('RAOB %s exited lateral boundaries of domain' % sid)
                        break
 
                    # If final time within this window, save final index for next time window 
                    k = k + 1
                    if (adpupa_t_s[idx_sid] > (wrf_hr[1] * 3600.)):
                        adpupa_idx[idx_sid] = k
 
                if debug > 1:
                    time_advect = dt.datetime.now()
                    entry_times.append((time_advect - kloop_start).total_seconds())
                    print('finished advecting radiosonde (%.6f s)' % (time_advect - time_interp).total_seconds())

        itime = itime + 1
    
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
        # Convert (x, y) coords back to km
        out_df['adpupa_x'] = out_df['adpupa_x'] * wrf_dx
        out_df['adpupa_y'] = out_df['adpupa_y'] * wrf_dx
        bufr.df_to_csv(out_df, fake_bufr_dir + t.strftime('/%Y%m%d%H%M.debug.adpupa.csv'))
    
    # Drop the extra columns we added
    out_df.drop(labels=extra_col_int, axis=1, inplace=True)
    out_df.drop(labels=extra_col_float, axis=1, inplace=True)

    # Use original indices to sort DataFrames
    out_df.sort_values('index', inplace=True)
    out_df.drop(labels=['index'], axis=1, inplace=True)
    bufr_csv.df.sort_values('index', inplace=True)
    bufr_csv.df.drop(labels=['index'], axis=1, inplace=True)

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
        print('avg time per entry = %.6f s' % np.mean(np.array(entry_times)))

# Total timing
print()
print('time for entire program = %.6f s' % (dt.datetime.now() - begin).total_seconds())


"""
End create_adpupa_obs.py
"""
