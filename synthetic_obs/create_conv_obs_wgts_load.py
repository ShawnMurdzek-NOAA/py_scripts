"""
Create Synthetic Conventional Observations from WRF Nature Run Output

Conventional obs are written into a new CSV file. UAS obs can then be added later. Conventional
observations include the following: ADPSFC, SFCSHP, MSONET, GPSIPW, ADPUPA, AIRCAR, AIRCFT.

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
    - 60 GB of memory on xjet is sufficient if interpolating 3D fields

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
import gc

import create_ob_utils as cou 
import bufr


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

work = '/mnt/lfs4'

# Directory containing wrfnat output from UPP
wrf_dir = work + '/BMC/wrfruc/murdzek/nature_run_spring_v2/'
#wrf_dir = work + '/BMC/wrfruc/murdzek/nature_run_tests/nature_run_spring_v2/output/202204291200/UPP/'

# Directory containing real prepbufr CSV output
bufr_dir = './'

# File containing conventional observation error characteristics (use GSI errtable file)
error_fname = work + '/BMC/wrfruc/murdzek/sample_real_obs/errtable.rrfs'

# Observation platforms to use (aka subsets, same ones used by BUFR)
ob_platforms = ['ADPUPA', 'AIRCAR', 'AIRCFT', 'PROFLR', 'ADPSFC', 'SFCSHP', 'MSONET', 'GPSIPW']

# Output directory for synthetic prepbufr CSV output
fake_bufr_dir = work + '/BMC/wrfruc/murdzek/nature_run_spring_v2/synthetic_obs/'
#fake_bufr_dir = work + '/BMC/wrfruc/murdzek/nature_run_tests/nature_run_spring_v2/synthetic_obs/'

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
    print('t = %s' % t.strftime('%Y-%m-%d %H:%M:%S'))
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
    hr_start = math.floor(bufr_csv.df['DHR'].min()*4) / 4
    hr_end = math.ceil(bufr_csv.df['DHR'].max()*4) / 4 + (wrf_step / 60.)
    wrf_ds = {}
    wrf_hr = np.arange(hr_start, hr_end + (wrf_step / 60.), wrf_step / 60.)
    for hr in wrf_hr:
        wrf_t = t + dt.timedelta(hours=hr)
        print(wrf_dir + wrf_t.strftime('wrfnat_%Y%m%d%H%M.grib2'))
        wrf_ds[hr] = xr.open_dataset(wrf_dir + wrf_t.strftime('wrfnat_%Y%m%d%H%M.grib2'), 
                                     engine='pynio')

    print('time to open GRIB files = %.2f s' % (dt.datetime.now() - start).total_seconds())
    
    # Extract latitude and longitude grids
    wrf_lat = wrf_ds[0]['gridlat_0'].values
    wrf_lon = wrf_ds[0]['gridlon_0'].values

    # Remove obs outside of the spatial domain of the wrfnat files
    latmin = wrf_lat.min()
    latmax = wrf_lat.max()
    lonmin = wrf_lon.min() + 360.
    lonmax = wrf_lon.max() + 360.
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

    # Determine row indices for each ob type
    ob_idx = {}
    obs = ['ADPSFC', 'SFCSHP', 'MSONET', 'GPSIPW', 'ADPUPA', 'AIRCAR', 'AIRCFT']
    for o in obs:
        ob_idx[o] = list(np.where(out_df['subset'] == o)[0]) 

    # Add some DataFrame columns
    nrow = len(out_df)
    extra_col = ['twgt', 'xi0', 'yi0', 'xwgt1', 'xwgt2', 'ywgt1', 'ywgt2', 'pi0', 'pwgt']
    for c in extra_col:
        out_df[c] = np.zeros(nrow)


    #-----------------------------------------------------------------------------------------------
    # Create Obs Based on 2-D Fields (ADPSFC, SFCSHP, MSONET, GPSIPW)
    #-----------------------------------------------------------------------------------------------

    # Extract 2D fields ONLY
    fields2D = ['PRES_P0_L1_GLC0', 'SPFH_P0_L103_GLC0', 'TMP_P0_L103_GLC0', 'HGT_P0_L1_GLC0', 
                'PWAT_P0_L200_GLC0']
    fields2D1 = ['UGRD_P0_L103_GLC0', 'VGRD_P0_L103_GLC0']
    wrf_data = {}
    for hr in wrf_hr:
        wrf_data[hr] = {}
        for f in fields2D:
            wrf_data[hr][f] = wrf_ds[hr][f].values
        for f in fields2D1:
            wrf_data[hr][f] = wrf_ds[hr][f][0, :, :].values
    
    # Loop over each ADPSFC, SFCSHP, MSONET, and GPSIPW observation
    for j in (ob_idx['ADPSFC'] + ob_idx['SFCSHP'] + ob_idx['MSONET'] + ob_idx['GPSIPW']):
        
        if debug > 1:
            time1 = dt.datetime.now()
            print()

        subset = out_df.iloc[j]
     
        # Determine WRF hour right before observation and weight for temporal interpolation
        ihr = np.where((wrf_hr - subset['DHR']) <= 0)[0][-1]
        twgt = (wrf_hr[ihr+1] - subset['DHR']) / (wrf_hr[ihr+1] - wrf_hr[ihr]) 

        if debug > 1:
            time2 = dt.datetime.now()
            print('done determining twgt (%.6f s)' % (time2 - time1).total_seconds())

        # Interpolate horizontally (code here is adapted from cou.wrf_coords)
        lat = subset['YOB']
        lon = subset['XOB'] - 360.
        close = np.unravel_index(np.argmin((wrf_lon - lon)**2 + (wrf_lat - lat)**2),
                                 wrf_lon.shape)
        clat = wrf_lat[close[0], close[1]]
        clon = wrf_lon[close[0], close[1]]

        # Check to make sure that we are not extrapolating by excluding (lat, lon) coordinates that
        # have an edge as their closest (lat, lon) coordinate in model space
        if (close[0] == 0 or close[0] == (wrf_lat.shape[0] - 1) or
            close[1] == 0 or close[1] == (wrf_lon.shape[1] - 1)):
            drop_idx.append(j)
            continue

        # Determine pseudo-bilinear interpolation weights in horizontal
        if lat < clat:
            yi0 = close[0] - 1
        else:
            yi0 = close[0]
        if lon < clon:
            xi0 = close[1] -1
        else:
            xi0 = close[1]
        ywgt1 = (wrf_lat[yi0+1, xi0] - lat) / (wrf_lat[yi0+1, xi0] - wrf_lat[yi0, xi0])
        xwgt1 = (wrf_lon[yi0, xi0+1] - lon) / (wrf_lon[yi0, xi0+1] - wrf_lon[yi0, xi0])
        ywgt2 = 1. - ywgt1
        xwgt2 = 1. - xwgt1

        if debug > 1:
            time3 = dt.datetime.now()
            print('done determining xi0, yi0 (%.6f s)' % (time3 - time2).total_seconds())

        # Determine surface height above sea level and surface pressure
        m = 'HGT_P0_L1_GLC0'
        sfch = (twgt * (xwgt1 * ywgt1 * wrf_data[wrf_hr[ihr]][m][yi0, xi0] +
                        xwgt2 * ywgt1 * wrf_data[wrf_hr[ihr]][m][yi0, xi0+1] +
                        xwgt1 * ywgt2 * wrf_data[wrf_hr[ihr]][m][yi0+1, xi0] +
                        xwgt2 * ywgt2 * wrf_data[wrf_hr[ihr]][m][yi0+1, xi0+1]) +
                (1.-twgt) * (xwgt1 * ywgt1 * wrf_data[wrf_hr[ihr+1]][m][yi0, xi0] +
                             xwgt2 * ywgt1 * wrf_data[wrf_hr[ihr+1]][m][yi0, xi0+1] +
                             xwgt1 * ywgt2 * wrf_data[wrf_hr[ihr+1]][m][yi0+1, xi0] +
                             xwgt2 * ywgt2 * wrf_data[wrf_hr[ihr+1]][m][yi0+1, xi0+1]))
        m = 'PRES_P0_L1_GLC0'
        sfcp = (twgt * (xwgt1 * ywgt1 * wrf_data[wrf_hr[ihr]][m][yi0, xi0] +
                        xwgt2 * ywgt1 * wrf_data[wrf_hr[ihr]][m][yi0, xi0+1] +
                        xwgt1 * ywgt2 * wrf_data[wrf_hr[ihr]][m][yi0+1, xi0] +
                        xwgt2 * ywgt2 * wrf_data[wrf_hr[ihr]][m][yi0+1, xi0+1]) +
                (1.-twgt) * (xwgt1 * ywgt1 * wrf_data[wrf_hr[ihr+1]][m][yi0, xi0] +
                             xwgt2 * ywgt1 * wrf_data[wrf_hr[ihr+1]][m][yi0, xi0+1] +
                             xwgt1 * ywgt2 * wrf_data[wrf_hr[ihr+1]][m][yi0+1, xi0] +
                             xwgt2 * ywgt2 * wrf_data[wrf_hr[ihr+1]][m][yi0+1, xi0+1])) * 1e-2

        if debug > 1:
            time4 = dt.datetime.now()
            print('done determining sfch and sfcp (%.6f s)' % (time4 - time3).total_seconds())

        # Surface obs
        if subset['subset'] in ['ADPSFC', 'SFCSHP', 'MSONET']:

            if debug > 1:
                print('2D field: nmsg = %d, subset = %s, TYP = %d' % (subset['nmsg'], 
                                                                      subset['subset'], 
                                                                      subset['TYP']))
                print('twgt = %.3f, xi0 = %d, yi0 = %d' % (twgt, xi0, yi0))
                print('DHR = %.3f, XOB = %.3f, YOB = %.3f' % (subset['DHR'], subset['XOB'], 
                                                              subset['YOB']))

            # Reset surface values to match NR and assign surface pressure values
            for o, v in zip(['ZOB', 'ELV', 'POB', 'PRSS'], [sfch, sfch, sfcp, sfcp]):
                if not np.isnan(subset[o]):
                    out_df.loc[j, o] = v

            # Interpolate temporally
            obs_name = ['QOB', 'TOB', 'UOB', 'VOB']
            wrf_name = ['SPFH_P0_L103_GLC0', 'TMP_P0_L103_GLC0', 'UGRD_P0_L103_GLC0', 
                        'VGRD_P0_L103_GLC0']
            for o, m in zip(obs_name, wrf_name):
                if not np.isnan(subset[o]):
                    if debug > 1:
                        time5 = dt.datetime.now()
                    out_df.loc[j, o] = (twgt * (xwgt1 * ywgt1 * wrf_data[wrf_hr[ihr]][m][yi0, xi0] +
                                                xwgt2 * ywgt1 * wrf_data[wrf_hr[ihr]][m][yi0, xi0+1] +
                                                xwgt1 * ywgt2 * wrf_data[wrf_hr[ihr]][m][yi0+1, xi0] +
                                                xwgt2 * ywgt2 * wrf_data[wrf_hr[ihr]][m][yi0+1, xi0+1]) +
                                        (1.-twgt) * (xwgt1 * ywgt1 * wrf_data[wrf_hr[ihr+1]][m][yi0, xi0] +
                                                     xwgt2 * ywgt1 * wrf_data[wrf_hr[ihr+1]][m][yi0, xi0+1] +
                                                     xwgt1 * ywgt2 * wrf_data[wrf_hr[ihr+1]][m][yi0+1, xi0] +
                                                     xwgt2 * ywgt2 * wrf_data[wrf_hr[ihr+1]][m][yi0+1, xi0+1]))
                    if debug > 1:
                        time6 = dt.datetime.now()
                        print('finished interp for %s (%.6f s)' % (o, (time6 - time5).total_seconds()))

        # GPS-Derived precipitable water
        elif subset['subset'] == 'GPSIPW':

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
                out_df.loc[j, 'PWO'] = (twgt * (xwgt1 * ywgt1 * wrf_data[wrf_hr[ihr]][m][yi0, xi0] +
                                                xwgt2 * ywgt1 * wrf_data[wrf_hr[ihr]][m][yi0, xi0+1] +
                                                xwgt1 * ywgt2 * wrf_data[wrf_hr[ihr]][m][yi0+1, xi0] +
                                                xwgt2 * ywgt2 * wrf_data[wrf_hr[ihr]][m][yi0+1, xi0+1]) +
                                        (1.-twgt) * (xwgt1 * ywgt1 * wrf_data[wrf_hr[ihr+1]][m][yi0, xi0] +
                                                     xwgt2 * ywgt1 * wrf_data[wrf_hr[ihr+1]][m][yi0, xi0+1] +
                                                     xwgt1 * ywgt2 * wrf_data[wrf_hr[ihr+1]][m][yi0+1, xi0] +
                                                     xwgt2 * ywgt2 * wrf_data[wrf_hr[ihr+1]][m][yi0+1, xi0+1]))
                if debug > 1:
                    time6 = dt.datetime.now()
                    print('finished interp for PWO (%.6f s)' % (time6 - time5).total_seconds())


    #-----------------------------------------------------------------------------------------------
    # Create Obs Based on 3-D Fields (ADPUPA, AIRCAR, AIRCFT)
    #-----------------------------------------------------------------------------------------------

    # Remove data from wrf_data that we don't need so we can free up some memory
    for hr in wrf_hr:
        for f in ['SPFH_P0_L103_GLC0', 'TMP_P0_L103_GLC0', 'PWAT_P0_L200_GLC0', 'UGRD_P0_L103_GLC0', 
                  'VGRD_P0_L103_GLC0']:
            del wrf_data[hr][f]

    # Call garbage collector to free up memory
    gc.collect()

    # We will extract 3D fields one at a time b/c these 3D arrays are massive (~6.5 GB each), so it 
    # is not feasible to load all of them at once. Ideally, only 2 should be loaded at any given
    # time. Therefore, unlike the 2D obs, we will loop over each obs time interval and each WRF
    # field.

    # Another approach that could have been taken here is to lazily load each 3D field using Xarray
    # datasets. While this would address the memory issues, it would create computational issues b/c
    # the dataset will have to be queried for each BUFR entry and each variable, which is time 
    # consuming. Tests comparing the use of Xarray datasets to regular Numpy arrays show that the
    # Numpy array approach is ~23x faster.

    # Loop over each WRF time interval
    for hr1, hr2 in zip(wrf_hr[:-1], wrf_hr[1:]): 

        # Determine indices of obs in this time interval
        ind = np.where(np.logical_and(out_df['DHR'] >= hr1, out_df['DHR'] < hr2))
        ind = np.intersect1d(ind, np.array(ob_idx['ADPUPA'] + ob_idx['AIRCAR'] + ob_idx['AIRCFT']))

        # If no indices, move to next time
        if ind.size == 0:
            continue

        # Loop over each variable (starting with pressure)
        obs_name = ['POB', 'TOB', 'QOB', 'ZOB', 'UOB', 'VOB']
        wrf_name = ['PRES_P0_L105_GLC0', 'TMP_P0_L105_GLC0', 'SPFH_P0_L103_GLC0', 
                    'HGT_P0_L105_GLC0', 'UGRD_P0_L105_GLC0', 'VGRD_P0_L105_GLC0']
        for o, f in zip(obs_name, wrf_name):
            print('3D Interp: %s' % f)

            # Extract field from UPP
            wrf_data[hr1]['v3d'] = wrf_ds[hr1][f].values
            wrf_data[hr2]['v3d'] = wrf_ds[hr2][f].values

            # Loop over each ADPUPA, AIRCAR, and AIRCFT observation within this time interval
            for j in ind:

                # Check to make sure that we didn't already drop this variable
                if j in drop_idx:
                    continue

                if debug > 1:
                    time1 = dt.datetime.now()
                    print()

                if o == 'POB':         

                    # Determine weight for temporal interpolation
                    twgt = (hr2 - out_df.loc[j, 'DHR']) / (hr2 - hr1) 

                    if debug > 1:
                        time2 = dt.datetime.now()
                        print('done determining twgt (%.6f s)' % (time2 - time1).total_seconds())

                    # Interpolate horizontally (code here is adapted from cou.wrf_coords)
                    lat = out_df.loc[j, 'YOB']
                    lon = out_df.loc[j, 'XOB'] - 360.
                    close = np.unravel_index(np.argmin((wrf_lon - lon)**2 + (wrf_lat - lat)**2),
                                             wrf_lon.shape)
                    clat = wrf_lat[close[0], close[1]]
                    clon = wrf_lon[close[0], close[1]]

                    # Check to make sure that we are not extrapolating by excluding (lat, lon) coordinates that
                    # have an edge as their closest (lat, lon) coordinate in model space
                    if (close[0] == 0 or close[0] == (wrf_lat.shape[0] - 1) or
                        close[1] == 0 or close[1] == (wrf_lon.shape[1] - 1)):
                        drop_idx.append(j)
                        continue

                    # Determine pseudo-bilinear interpolation weights in horizontal
                    if lat < clat:
                        yi0 = close[0] - 1
                    else:
                        yi0 = close[0]
                    if lon < clon:
                        xi0 = close[1] -1
                    else:
                        xi0 = close[1]
                    ywgt1 = (wrf_lat[yi0+1, xi0] - lat) / (wrf_lat[yi0+1, xi0] - wrf_lat[yi0, xi0])
                    xwgt1 = (wrf_lon[yi0, xi0+1] - lon) / (wrf_lon[yi0, xi0+1] - wrf_lon[yi0, xi0])
                    ywgt2 = 1. - ywgt1
                    xwgt2 = 1. - xwgt1

                    if debug > 1:
                        time3 = dt.datetime.now()
                        print('done determining xi0, yi0 (%.6f s)' % (time3 - time2).total_seconds())

                    # Determine surface height above sea level and surface pressure
                    m = 'HGT_P0_L1_GLC0'
                    sfch = (twgt * (xwgt1 * ywgt1 * wrf_data[hr1][m][yi0, xi0] +
                                    xwgt2 * ywgt1 * wrf_data[hr1][m][yi0, xi0+1] +
                                    xwgt1 * ywgt2 * wrf_data[hr1][m][yi0+1, xi0] +
                                    xwgt2 * ywgt2 * wrf_data[hr1][m][yi0+1, xi0+1]) +
                            (1.-twgt) * (xwgt1 * ywgt1 * wrf_data[hr2][m][yi0, xi0] +
                                         xwgt2 * ywgt1 * wrf_data[hr2][m][yi0, xi0+1] +
                                         xwgt1 * ywgt2 * wrf_data[hr2][m][yi0+1, xi0] +
                                         xwgt2 * ywgt2 * wrf_data[hr2][m][yi0+1, xi0+1]))
                    m = 'PRES_P0_L1_GLC0'
                    sfcp = (twgt * (xwgt1 * ywgt1 * wrf_data[hr1][m][yi0, xi0] +
                                    xwgt2 * ywgt1 * wrf_data[hr1][m][yi0, xi0+1] +
                                    xwgt1 * ywgt2 * wrf_data[hr1][m][yi0+1, xi0] +
                                    xwgt2 * ywgt2 * wrf_data[hr1][m][yi0+1, xi0+1]) +
                            (1.-twgt) * (xwgt1 * ywgt1 * wrf_data[hr2][m][yi0, xi0] +
                                         xwgt2 * ywgt1 * wrf_data[hr2][m][yi0, xi0+1] +
                                         xwgt1 * ywgt2 * wrf_data[hr2][m][yi0+1, xi0] +
                                         xwgt2 * ywgt2 * wrf_data[hr2][m][yi0+1, xi0+1]))
            
                    # Check whether ob from BUFR file lies underground
                    if (out_df.loc[j, 'ZOB'] < sfch) or (out_df.loc[j, 'POB'] > sfcp):
                        drop_idx.append(j)
                        continue

                    if debug > 1:
                        time4 = dt.datetime.now()
                        print('done determining sfch and sfcp (%.6f s)' % (time4 - time3).total_seconds())
            
                    # Perform vertical interpolation in pressure rather than height b/c some obs don't 
                    # have height 
                    m = 'PRES_P0_L105_GLC0'
                    p1d = (twgt * (xwgt1 * ywgt1 * wrf_data[hr1][m][yi0, xi0] +
                                   xwgt2 * ywgt1 * wrf_data[hr1][m][yi0, xi0+1] +
                                   xwgt1 * ywgt2 * wrf_data[hr1][m][yi0+1, xi0] +
                                   xwgt2 * ywgt2 * wrf_data[hr1][m][yi0+1, xi0+1]) +
                           (1.-twgt) * (xwgt1 * ywgt1 * wrf_data[hr2][m][yi0, xi0] +
                                        xwgt2 * ywgt1 * wrf_data[hr2][m][yi0, xi0+1] +
                                        xwgt1 * ywgt2 * wrf_data[hr2][m][yi0+1, xi0] +
                                        xwgt2 * ywgt2 * wrf_data[hr2][m][yi0+1, xi0+1])) * 1e-2
                    pi0 = np.where(p1d > out_df.loc[j, 'POB'])[0][-1]
                    pwgt = (p1d[pi0+1] - out_df.loc[j, 'POB']) / (p1d[pi0+1] - p1d[pi0])
 
                    if debug > 1:
                        time5 = dt.datetime.now()
                        print('done determining pi0 = %d (%.6f s)' % (pi0, (time5 - time4).total_seconds()))

                else:
                    # Extract indicies and weights saved from POB
                    twgt = out_df.loc[j, 'twgt'] 
                    xi0 = out_df.loc[j, 'xi0'] 
                    yi0 = out_df.loc[j, 'yi0'] 
                    xwgt1 = out_df.loc[j, 'xwgt1'] 
                    xwgt2 = out_df.loc[j, 'xwgt2'] 
                    ywgt1 = out_df.loc[j, 'ywgt1'] 
                    ywgt2 = out_df.loc[j, 'ywgt2'] 
                    pi0 = out_df.loc[j, 'pi0'] 
                    pwgt = out_df.loc[j, 'pwgt'] 

                # Interpolate in 4 dimensions
                if not np.isnan(subset[o]):
                    if debug > 1:
                        time6 = dt.datetime.now()
                    v1 = (twgt * (xwgt1 * ywgt1 * wrf_data[hr1][m][pi0, yi0, xi0] +
                                  xwgt2 * ywgt1 * wrf_data[hr1][m][pi0, yi0, xi0+1] +
                                  xwgt1 * ywgt2 * wrf_data[hr1][m][pi0, yi0+1, xi0] +
                                  xwgt2 * ywgt2 * wrf_data[hr1][m][pi0, yi0+1, xi0+1]) +
                          (1.-twgt) * (xwgt1 * ywgt1 * wrf_data[hr2][m][pi0, yi0, xi0] +
                                       xwgt2 * ywgt1 * wrf_data[hr2][m][pi0, yi0, xi0+1] +
                                       xwgt1 * ywgt2 * wrf_data[hr2][m][pi0, yi0+1, xi0] +
                                       xwgt2 * ywgt2 * wrf_data[hr2][m][pi0, yi0+1, xi0+1]))
                    v2 = (twgt * (xwgt1 * ywgt1 * wrf_data[hr1][m][pi0+1, yi0, xi0] +
                                  xwgt2 * ywgt1 * wrf_data[hr1][m][pi0+1, yi0, xi0+1] +
                                  xwgt1 * ywgt2 * wrf_data[hr1][m][pi0+1, yi0+1, xi0] +
                                  xwgt2 * ywgt2 * wrf_data[hr1][m][pi0+1, yi0+1, xi0+1]) +
                          (1.-twgt) * (xwgt1 * ywgt1 * wrf_data[hr2][m][pi0+1, yi0, xi0] +
                                       xwgt2 * ywgt1 * wrf_data[hr2][m][pi0+1, yi0, xi0+1] +
                                       xwgt1 * ywgt2 * wrf_data[hr2][m][pi0+1, yi0+1, xi0] +
                                       xwgt2 * ywgt2 * wrf_data[hr2][m][pi0+1, yi0+1, xi0+1]))
                    out_df.loc[j, o] = (v1**pwgt) * (v2**(1.-pwgt))
                    if debug > 1:
                        time7 = dt.datetime.now()
                        print('finished interp for %s (%.6f s)' % (o, (time7 - time6).total_seconds()))

                # Add errors

                if debug > 1:
                    print('total time = %.6f s' % (dt.datetime.now() - time1).total_seconds())
                    entry_times.append((dt.datetime.now() - time1).total_seconds())

                # Memory management
                del wrf_data[hr1]['v3d']
                del wrf_data[hr2]['v3d']
                gc.collect()

    # Drop rows that we skipped as well as the extra columns we added
    out_df.drop(labels=extra_col, inplace=True)
    out_df.drop(index=drop_idx, inplace=True)
    out_df.reset_index(drop=True, inplace=True)
    bufr_csv.df.drop(index=drop_idx, inplace=True)
    bufr_csv.df.reset_index(drop=True, inplace=True)

    # Convert to proper units
    out_df['QOB'] = out_df['QOB'] * 1e6
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
