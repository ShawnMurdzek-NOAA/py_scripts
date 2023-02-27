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
    - 25 GB of memory on sjet is sufficient if interpolating 3D fields

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
import os
import scipy.spatial as ss
import pickle

import create_ob_utils as cou 
import bufr
import map_proj as mp


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

work = '/mnt/lfs4'

# Directory containing wrfnat output from UPP
wrf_dir = work + '/BMC/wrfruc/murdzek/nature_run_spring_v2/'
#wrf_dir = work + '/BMC/wrfruc/murdzek/nature_run_tests/nature_run_spring_v2/output/202204291200/UPP/'

# Directory containing real prepbufr CSV output
bufr_dir = work + '/BMC/wrfruc/murdzek/sample_real_obs/obs_rap/'
#bufr_dir = work + '/BMC/wrfruc/murdzek/py_scripts/synthetic_obs/'
#bufr_dir = work + '/BMC/wrfruc/murdzek/src/py_scripts/synthetic_obs/'
#bufr_dir = work + '/BMC/wrfruc/murdzek/sample_real_obs/test/'

# Observation platforms to use (aka subsets, same ones used by BUFR)
ob_platforms = ['ADPUPA', 'AIRCAR', 'AIRCFT', 'ADPSFC', 'SFCSHP', 'MSONET', 'GPSIPW']
#ob_platforms = ['AIRCAR', 'AIRCFT', 'ADPSFC', 'SFCSHP', 'MSONET', 'GPSIPW']

# Output directory for synthetic prepbufr CSV output
fake_bufr_dir = work + '/BMC/wrfruc/murdzek/nature_run_spring_v2/synthetic_obs/'
#fake_bufr_dir = work + '/BMC/wrfruc/murdzek/nature_run_tests/nature_run_spring_v2/synthetic_obs/'

# Start and end times for prepbufrs. Step is in min
bufr_start = dt.datetime(2022, 4, 29, 12)
bufr_end = dt.datetime(2022, 4, 29, 13)
bufr_step = 120

# Start and end times for wrfnat UPP output. Step is in min
wrf_start = dt.datetime(2022, 4, 29, 12, 0)
wrf_end = dt.datetime(2022, 4, 29, 15, 0)
wrf_step = 15

# Option for debugging output (0 = none, 1 = some, 2 = a lot)
debug = 2

# Option to interpolate (lat, lon) coordinates for surface obs (ADPSFC, SFCSHP, MSONET)
# Helpful for debugging, but should usually be set to False b/c it increases runtime
interp_latlon = True


#---------------------------------------------------------------------------------------------------
# Functions (Should move to separate file at some point)
#---------------------------------------------------------------------------------------------------

def determine_twgt(wrf_hr, dhr):
    """
    Compute the weight and index needed for time interpolation using 1-D linear interpolation

    Inputs
    ------
    wrf_hr : array
        Decimal hours for each WRF UPP file
    dhr : float
        Decimal hour for the given observation

    Returns
    -------
    ihr : integer
        Index for wrf_hr associated with twgt 
    twgt : float
        Interpolation weight associated with ihr

    """
       
    ihr = np.where((wrf_hr - dhr) <= 0)[0][-1]
    twgt = (wrf_hr[ihr+1] - dhr) / (wrf_hr[ihr+1] - wrf_hr[ihr])

    return ihr, twgt
 

def _bilinear_interp_horiz(field, iwgt, jwgt, i0, j0, threeD=False):
    """
    Interpolate the given field in the horizontal

    Inputs
    ------
    field : array 
        2D field to interpolate
    iwgt : float
        Interpolation weight for i0
    jwgt : float
        Interpolation weight fro j0
    i0 : integer
        Lower-left index for the first dimension of field
    j0 : integer
        Lower-left index for the second dimension of field
    threeD : boolean, optional
        Is this field actually 3D?

    Returns
    -------
    val : float
        Interpolated value

    """

    iind = np.array([i0,        i0,            i0+1,          i0+1])
    jind = np.array([j0,        j0+1,          j0,            j0+1])
    wgts = np.array([iwgt*jwgt, iwgt*(1-jwgt), (1-iwgt)*jwgt, (1-iwgt)*(1-jwgt)])

    if threeD:
        val = np.zeros(field.shape[0])
        for w, i, j in zip(wgts, iind, jind):
            val = val + w * field[:, i, j]
    else:
        val = 0.
        for w, i, j in zip(wgts, iind, jind):
            val = val + w * field[i, j]

    return val


def _log_interp(val1, val2, wgt1):

    return (val1**wgt1) * (val2**(1.-wgt1))


def _interp_x_y_t(field1, field2, iwgt, jwgt, twgt, i0, j0, threeD=False):
    """
    Interpolate linearly in the horizontal and time dimensions

    Inputs
    ------
    field1 : array
        2D field to interpolate with time weight = twgt
    field2 : array
        2D field to interpoalte with time weight = 1 - twgt
    iwgt, jwgt, twgt : float
        Interpolation weights for the (i0, j0) gridpoint in field1
    i0, j0 : integers
        Lower-left index for field1 and field2
    threeD : boolean, optional
        Is this field actually 3D?

    Returns
    -------
    val : float
        Interpolated value

    """

    val = (twgt * _bilinear_interp_horiz(field1, iwgt, jwgt, i0, j0, threeD=threeD) + 
           (1.-twgt) * _bilinear_interp_horiz(field2, iwgt, jwgt, i0, j0, threeD=threeD))

    return val     


def _interp_x_y_z(field, iwgt, jwgt, pwgt, i0, j0, pi0):
    """
    Interpolate linearly in the horizontal dimensions and logarithmically in pressure

    Inputs
    ------
    field : array
        3D field to interpolate
    iwgt, jwgt, pwgt : float
        Interpolation weights for the (pi0, i0, j0) gridpoint in field
    i0, j0, pi0 : integers
        Lower-left index for field

    Returns
    -------
    val : float
        Interpolated value

    """

    val = _log_interp(_bilinear_interp_horiz(field[pi0, :, :], iwgt, jwgt, i0, j0),
                      _bilinear_interp_horiz(field[pi0+1, :, :], iwgt, jwgt, i0, j0), pwgt)

    return val     


def interp_wrf_to_obs(wrf_data, wrf_hr, var, ob_subset, ihr, twgt, threeD=False):
    """
    Wrapper function for interpolation in x, y, and t

    Inputs
    ------
    wrf_data : dictionary
        Dictionary of Xarray datasets containing UPP output data
    wrf_hr : list
        List of valid times (in decimal hours) for wrf_data. These are the keys for wrf_data
    var : string
        Variable from the UPP dataset to interpolate
    ob_subset : series
        Pandas series containing a single observation
    ihr : integer
        Index for the UPP dataset immediately before the observation time
    twgt : float
        Interpolation weight associated with the ihr UPP dataset
    threeD : boolean, optional
        Is this UPP field 3D?

    Returns
    -------
    val : float
        UPP output linearly interpolated in x, y, and t to the observation location

    """

    val = _interp_x_y_t(wrf_data[wrf_hr[ihr]][var], wrf_data[wrf_hr[ihr+1]][var], twgt, 
                        ob_subset['iwgt'], ob_subset['jwgt'], ob_subset['i0'], ob_subset['j0'], 
                        threeD=threeD)

    return val


def interp_wrf_3d(wrf3d, ob_subset):
    """
    Wrapper function for interpolation in x, y, and p

    Inputs
    ------
    wrf3d : array
        WRF 3D output array to interpolate
    ob_subset : series
        Pandas series containing a single observation

    Returns
    -------
    val : float
        UPP output linearly interpolated in x, y, and p to the observation location

    """

    val = _interp_x_y_z(wrf3d, ob_subset['iwgt'], ob_subset['jwgt'], ob_subset['pwgt'], 
                        ob_subset['i0'], ob_subset['j0'], ob_subset['pi0'])

    return val


def interp_wrf_p1d(p1d, ob_subset):
    """
    Wrapper function for interpolation of pressure in one dimension

    Inputs
    ------
    p1d : array
        1D pressure array to interpolate
    ob_subset : series
        Pandas series containing a single observation

    Returns
    -------
    val : float
        UPP output logarithmically interpolated in p to the observation location
    pwgt : float
        Weight used to logarithmic interpolation

    """

    pi0 = ob_subset['pi0']
    pwgt = (p1d[pi0+1] - ob_subset['POB']) / (p1d[pi0+1] - p1d[pi0])
    val = _log_interp(p1d[pi0], p1d[pi0+1], pwgt)

    return val, pwgt


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
        print(wrf_dir + wrf_t.strftime('wrfnat_%Y%m%d%H%M.grib2'))
        wrf_ds[hr] = xr.open_dataset(wrf_dir + wrf_t.strftime('wrfnat_%Y%m%d%H%M.grib2'), 
                                     engine='pynio')

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
            time1 = dt.datetime.now()
            print()

        subset = out_df.loc[j]
     
        # Determine WRF hour right before observation and weight for temporal interpolation
        ihr, twgt = determine_twgt(wrf_hr, subset['DHR'])

        if debug > 1:
            time2 = dt.datetime.now()
            print('done determining twgt (%.6f s)' % (time2 - time1).total_seconds())

        # Determine surface height above sea level and surface pressure
        sfch = interp_wrf_to_obs(wrf_data, wrf_hr, 'HGT_P0_L1_GLC0', subset, ihr, twgt)
        sfcp = interp_wrf_to_obs(wrf_data, wrf_hr, 'PRES_P0_L1_GLC0', subset, ihr, twgt)

        if debug > 1:
            time4 = dt.datetime.now()
            print('sfch = %.6f' % sfch)
            print('sfcp = %.6f' % sfcp)
            print('done determining sfch and sfcp (%.6f s)' % (time4 - time2).total_seconds())

        # Surface obs
        if subset['subset'] in ['ADPSFC', 'SFCSHP', 'MSONET']:

            if debug > 1:
                print('2D field: nmsg = %d, subset = %s, TYP = %d' % (subset['nmsg'], 
                                                                      subset['subset'], 
                                                                      subset['TYP']))
                print('twgt = %.3f, i0 = %d, j0 = %d' % (twgt, subset['i0'], subset['j0']))
                print('DHR = %.3f, XOB = %.6f, YOB = %.6f' % (subset['DHR'], subset['XOB'], 
                                                              subset['YOB']))

            # Reset surface values to match NR and assign surface pressure values
            for o, v in zip(['ZOB', 'ELV', 'POB', 'PRSS'], [sfch, sfch, sfcp, sfcp]):
                if not np.isnan(subset[o]):
                    out_df.loc[j, o] = v

            # Interpolate temporally
            obs_name = ['QOB', 'TOB', 'UOB', 'VOB']
            wrf_name = ['SPFH_P0_L103_GLC0', 'TMP_P0_L103_GLC0', 'UGRD_P0_L103_GLC0', 
                        'VGRD_P0_L103_GLC0']
            if interp_latlon:
                obs_name.append('XOB')
                obs_name.append('YOB')
                wrf_name.append('gridlon_0')
                wrf_name.append('gridlat_0')
            for o, m in zip(obs_name, wrf_name):
                if not np.isnan(subset[o]):
                    if debug > 1:
                        time5 = dt.datetime.now()
                    out_df.loc[j, o] = interp_wrf_to_obs(wrf_data, wrf_hr, m, subset, ihr, twgt)
                    if debug > 1:
                        time6 = dt.datetime.now()
                        print('finished interp for %s (%.6f s)' % (o, (time6 - time5).total_seconds()))
                        entry_times2d.append((dt.datetime.now() - time1).total_seconds())

        # GPS-Derived precipitable water
        elif subset['subset'] == 'GPSIPW':

            if debug > 1:
                print('PW field: nmsg = %d, subset = %s, TYP = %d' % (subset['nmsg'], 
                                                                      subset['subset'], 
                                                                      subset['TYP']))
                print('twgt = %.3f, xi = %.3f, yi = %.3f' % (twgt, xi, yi))

            # Interpolate temporally
            if not np.isnan(subset['PWO']):
                if debug > 1:
                    time5 = dt.datetime.now()
                out_df.loc[j, 'PWO'] = interp_wrf_to_obs(wrf_data, wrf_hr, 'PWAT_P0_L200_GLC0', 
                                                         subset, ihr, twgt)
                if debug > 1:
                    time6 = dt.datetime.now()
                    print('finished interp for PWO (%.6f s)' % (time6 - time5).total_seconds())
                    entry_times2d.append((dt.datetime.now() - time1).total_seconds())


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
                    time1 = dt.datetime.now()
                    print()

                if (o == 'POB' and np.isclose(p1d[0, j], 0.)):         

                    # Determine weight for temporal interpolation
                    twgt = 1. - (np.abs(hr - out_df.loc[j, 'DHR']) / wrf_step_dec) 

                    if debug > 1:
                        time2 = dt.datetime.now()
                        print('done determining twgt (%.3f, %.6f s)' % (twgt, (time2 - time1).total_seconds()))
        
                    # Determine surface height above sea level and surface pressure
                    if hr > out_df.loc[j, 'DHR']:
                        hr1 = hr - wrf_step_dec
                        hr2 = hr
                    else:
                        hr1 = hr
                        hr2 = hr + wrf_step_dec
                    sfch = interp_wrf_to_obs(wrf_data, wrf_hr, 'HGT_P0_L1_GLC0', out_df.loc[j],
                                             ihr, twgt)
                    sfcp = interp_wrf_to_obs(wrf_data, wrf_hr, 'PRES_P0_L1_GLC0', out_df.loc[j], 
                                             ihr, twgt)
            
                    # Check whether ob from BUFR file lies underground
                    if (out_df.loc[j, 'ZOB'] < sfch) or (out_df.loc[j, 'POB'] > sfcp):
                        drop_idx.append(j)
                        continue

                    if debug > 1:
                        time4 = dt.datetime.now()
                        print('done determining sfch and sfcp (%.6f s)' % (time4 - time2).total_seconds())
            
                    # Perform vertical interpolation in pressure rather than height b/c some obs don't 
                    # have height. The calculation below only gives use part of the p1d array. The 
                    # rest of this of this array will be determined when we loop over the other 3D
                    # pressure array
                    p1d[:, j] = twgt * 1e-2 * _bilinear_interp_horiz(wrf3d, out_df.loc[j, 'iwgt'], 
                                                                     out_df.loc[j, 'jwgt'],
                                                                     out_df.loc[j, 'i0'], out_df.loc[j, 'j0'],
                                                                     threeD=True)
 
                    if debug > 1:
                        time5 = dt.datetime.now()
                        print('done computing first part of p1d (%.6f s)' % (time5 - time4).total_seconds())

                    # Special case: twgt = 1. In this case, we don't need to interpolate in time, so 
                    # we can skip the P2 section of the conditional
                    if np.isclose(twgt, 1):
                        # Check for extrapolation
                        if (p1d[0, j] > out_df.loc[j, 'POB']):
                            out_df.loc[j, 'pi0'] = np.where(p1d[:, j] > out_df.loc[j, 'POB'])[0][-1]
                            if debug > 1:
                                print('pi0 = %d' % out_df.loc[j, 'pi0'])
                            if out_df.loc[j, 'pi0'] >= (p1d.shape[0] - 1):
                                # Prevent extrapolation in vertical
                                drop_idx.append(j)
                                continue
                            out_df.loc[j, 'POB'], out_df.loc[j, 'pwgt'] = interp_wrf_p1d(p1d[:, j], out_df.loc[j])
                        else:
                            drop_idx.append(j)
                            continue

                        if debug > 1: 
                            print('interpolated P = %.2f' % out_df.loc[j, 'POB'])
                            print('actual P = %.2f' % bufr_csv.df.loc[j, 'POB'])

                    # save some stuff
                    out_df.loc[j, 'twgt']  = twgt 

                    if debug > 1:
                        print('total time for p1 = %.6f s' % (dt.datetime.now() - time1).total_seconds())
                        entry_timesp1.append((dt.datetime.now() - time1).total_seconds())

                elif o == 'POB':

                    # Extract indicies and weights saved from first POB
                    twgt =  out_df.loc[j, 'twgt'] 

                    # Finish computing p1d
                    p1d[:, j] = p1d[:, j] + (1.-twgt) * 1e-2 * _bilinear_interp_horiz(wrf3d, out_df.loc[j, 'iwgt'], 
                                                                    out_df.loc[j, 'jwgt'],
                                                                    out_df.loc[j, 'i0'], out_df.loc[j, 'j0'],
                                                                    threeD=True)

                    # Check for extrapolation
                    if (p1d[0, j] > out_df.loc[j, 'POB']):
                        out_df.loc[j, 'pi0'] = np.where(p1d[:, j] > out_df.loc[j, 'POB'])[0][-1]
                        if debug > 1:
                            print('pi0 = %d' % out_df.loc[j, 'pi0'])
                        if out_df.loc[j, 'pi0'] >= (p1d.shape[0] - 1):
                            # Prevent extrapolation in vertical
                            drop_idx.append(j)
                            continue
                        out_df.loc[j, 'POB'], out_df.loc[j, 'pwgt'] = interp_wrf_p1d(p1d[:, j], out_df.loc[j])
                    else:
                        drop_idx.append(j)
                        continue

                    if debug > 1:
                        print('interpolated P = %.2f' % out_df.loc[j, 'POB'])
                        print('actual P = %.2f' % bufr_csv.df.loc[j, 'POB'])
                        print('total time for p2 = %.6f s' % (dt.datetime.now() - time1).total_seconds())
                        entry_timesp2.append((dt.datetime.now() - time1).total_seconds())

                else:

                    # Interpolate
                    if not np.isnan(out_df.loc[j, o]):
                        if debug > 1:
                            time6 = dt.datetime.now()
                        twgt =  out_df.loc[j, 'twgt'] 
                        if np.isclose(out_df.loc[j, o], 0):
                            out_df.loc[j, o] = twgt * interp_wrf_3d(wrf3d, out_df.loc[j])
                        else:
                            out_df.loc[j, o] = ((1.-twgt) * interp_wrf_3d(wrf3d, out_df.loc[j]) + 
                                                out_df.loc[j, o])
                        if debug > 1:
                            time7 = dt.datetime.now()
                            print('finished interp for %s (%.6f s)' % (o, (time7 - time6).total_seconds()))

                    if debug > 1:
                        print('total time = %.6f s' % (dt.datetime.now() - time1).total_seconds())
                        entry_times3d.append((dt.datetime.now() - time1).total_seconds())

            # Free up memory (shouldn't have to call garbage collector after this)
            wrf3d = 0.

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
