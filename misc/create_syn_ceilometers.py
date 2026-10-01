"""
Create synthetic ceilometer observations from GRIB2 model output

Environment: /work2/noaa/wrfruc/murdzek/conda/miniforge_hercules/env/grib2io

Shawn Murdzek
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import sys
import datetime as dt
import argparse
import copy
import numpy as np
import pandas as pd
import grib2io


#---------------------------------------------------------------------------------------------------
# Main Program
#---------------------------------------------------------------------------------------------------

def parse_in_args(argv):
    """
    Parse command-line arguments from sys.argv[0]
    """

    parser = argparse.ArgumentParser(description='Script that creates synthetic ceilometer \
                                                  observations by interpolating model GRIB2 \
                                                  output to user-specific ceilometer locations. \
                                                  Output is saved to a DART-compatible obs_seq \
                                                  file.')

    # Positional arguments
    parser.add_argument('NR_file',
                        help='Nature run native-level UPP output in GRIB2 format',
                        type=str)

    parser.add_argument('ceil_file',
                        help='Ceilometer station locations in CSV format. Must contain 4 columns: \
                              SID, lat, lon, and error.',
                        type=str)

    # Optional arguments
    parser.add_argument('--min_range',
                        dest='min_range',
                        default=0,
                        help='Minimum ceilometer range (m AGL)',
                        type=float)

    parser.add_argument('--max_range',
                        dest='max_range',
                        default=3657.6,
                        help='Maximum ceilometer range (m AGL)',
                        type=float)

    parser.add_argument('--min_cldfrac_thres',
                        dest='min_cldfrac_thres',
                        default=0.05,
                        help='Minimum detectable cloud fraction (decimal)',
                        type=float)

    parser.add_argument('--clear_ob_interval',
                        dest='clear_ob_interval',
                        default=2000,
                        help='Interval between clear observations (m AGL)',
                        type=float)

    parser.add_argument('--out_file',
                        dest='out_file',
                        default='syn_ceil.out',
                        help='Output file name',
                        type=str)

    return parser.parse_args(argv)


def reconstruct_3d_hyb_lvl_field(grib2_obj, field, nz):
    """
    Reconstruct a 3D field from 2D horizontal slices from grib2_obj

    Assumes field has hybrid-level coordinates
    """

    # Construct array to save output
    msg_sample = grib2_obj.select(shortName=field)[0]
    out = np.zeros((nz, msg_sample.ny, msg_sample.nx))

    # Fill output array
    for i in range(nz):
        msg = grib2_obj.select(shortName=field, level=f"{i+1} hybrid level")[0]
        out[i, :, :] = msg.data

    return out


def compute_cld_hgt_amt(grib2_obj, min_rng=0, max_rng=3657.6, min_thres=0.05):
    """
    Return two arrays with the height and amount of the lowest cloud layer
    """

    # Determine number of hybrid levels
    cldfrac_field = 'FRACCC'
    nz = len(grib2_obj.levels_by_var(cldfrac_field))

    # Get 3D fields
    cldfrac_3d = reconstruct_3d_hyb_lvl_field(grib2_obj, cldfrac_field, nz)
    sfc_hgt = grib2_obj.select(shortName='HGT', level='surface')[0].data
    hgt_agl_3d = reconstruct_3d_hyb_lvl_field(grib2_obj, 'HGT', nz) - sfc_hgt

    # Initialize output arrays as "clear" obs (amt = 0, hgt = NaN)
    cld_amt = np.zeros(sfc_hgt.shape)
    cld_hgt = np.ones(sfc_hgt.shape) * np.nan

    # Determine height indices to loop over
    max_hgt_1d = np.amax(hgt_agl_3d, axis=(1, 2))
    min_hgt_1d = np.amin(hgt_agl_3d, axis=(1, 2))
    loop_idx = np.where((min_hgt_1d >= min_rng) & (max_hgt_1d <= max_rng))[0]

    # Set cloud fraction to NaN outside of ceilometer range
    cldfrac_3d[(hgt_agl_3d < min_rng) | (hgt_agl_3d > max_rng)] = np.nan

    # Set cloud fraction to clear below min detectable threshold
    cldfrac_3d[cldfrac_3d <= min_thres] = 0

    # For each vertical level, search for (i, j) locations where cldfrac > 0 and a cloud height
    # has not been assigned yet. Assign cldfrac for these locations
    for k in loop_idx:
        cldfrac_slice = cldfrac_3d[k, :, :]
        hgt_slice = hgt_agl_3d[k, :, :]
        cond = (cldfrac_slice > 0) & np.isnan(cld_hgt)
        cld_amt[cond] = cldfrac_slice[cond]
        cld_hgt[cond] = hgt_slice[cond]

    return cld_amt, cld_hgt


def create_ceil_obs(grib2_obj, sites, min_rng=0, max_rng=3657.6, min_thres=0.05):
    """
    Create ceilometer observations by interpolating NR output to ceilometer locations
    """

    # Compute 2D arrays of cloud height and cloud amount
    cld_amt, cld_hgt = compute_cld_hgt_amt(grib2_obj, min_rng=min_rng, max_rng=max_rng, min_thres=min_thres)

    # Perform interpolation
    # Section 3 in the GRIB2 message contains geospatial information needed for interpolation
    sect3 = grib2_obj.select(shortName='HGT', level='surface')[0].section3
    grid_def = grib2io.Grib2GridDef.from_section3(sect3)
    ceil_lat = sites['lat'].values
    ceil_lon = sites['lon'].values
    ceil_amt = grib2io.interpolate_to_stations(cld_amt, 'neighbor', grid_def, ceil_lat, ceil_lon)
    ceil_hgt = grib2io.interpolate_to_stations(cld_hgt, 'neighbor', grid_def, ceil_lat, ceil_lon)

    # Save output to DataFrame
    ceil_out = copy.deepcopy(sites)
    ceil_out['cld_amt'] = ceil_amt
    ceil_out['cld_hgt'] = ceil_hgt

    return ceil_out


def convert_to_oktas(obs):
    """
    Convert ceilometer obs from fractional cloud amount to oktas

    Source: https://worldweather.wmo.int/oktas.htm

    Okta definitions are a big ambiguous, so I decided to center the bins with the exception of
    0/8 and 8/8, which are only used if there are no clouds or all clouds, respectively
    """

    obs.loc[(obs['cld_amt'] > 0) & (obs['cld_amt'] <= 3/16), 'cld_amt'] = 1/8
    obs.loc[(obs['cld_amt'] > 3/16) & (obs['cld_amt'] <= 5/16), 'cld_amt'] = 2/8
    obs.loc[(obs['cld_amt'] > 5/16) & (obs['cld_amt'] <= 7/16), 'cld_amt'] = 3/8
    obs.loc[(obs['cld_amt'] > 7/16) & (obs['cld_amt'] <= 9/16), 'cld_amt'] = 4/8
    obs.loc[(obs['cld_amt'] > 9/16) & (obs['cld_amt'] <= 11/16), 'cld_amt'] = 5/8
    obs.loc[(obs['cld_amt'] > 11/16) & (obs['cld_amt'] <= 13/16), 'cld_amt'] = 6/8
    obs.loc[(obs['cld_amt'] > 13/16) & (obs['cld_amt'] < 1), 'cld_amt'] = 7/8

    return obs


def add_clear_obs_hgts(obs, clr_ob_interval=500, min_rng=0, max_rng=3657.6):
    """
    If a station is clear, add clear obs at regular intervals

    Do not add clear obs below cloudy obs
    """

    # Define height of clear obs
    clr_ob_hgt = np.arange(min_rng + 0.5*clr_ob_interval, max_rng, clr_ob_interval)

    # Split obs DataFrame into clear and cloudy obs
    cond = obs['cld_amt'] == 0
    obs_clr = copy.deepcopy(obs.loc[cond, :])
    obs_cld = copy.deepcopy(obs.loc[~cond, :])

    # Add clear obs
    if len(clr_ob_hgt) > 0:
        list_df = [obs_cld]
        for hgt in clr_ob_hgt:
            tmp = copy.deepcopy(obs_clr)
            tmp.loc[:, 'cld_hgt'] = hgt
            list_df.append(tmp)
        obs = pd.concat(list_df)
    else:
        obs = obs_cld

    # Sort output
    obs.sort_values('SID', axis=0, inplace=True)
    obs.reset_index(inplace=True)

    return obs


def write_ceil_to_obs_seq(obs, out_fname, valid_dt, clr_ob_interval=500):
    """
    Write ceilometer obs to DART obs_seq file

    Adapted from a function that writes synthetic UAS obs to obs_seq format
    """

    # Compute time in DART format (seconds and days since 0000 UTC 1 Jan 1601)
    time_diff = valid_dt - dt.datetime(1601, 1, 1)
    time_sec = time_diff.seconds
    time_days = time_diff.days

    # Convert lat and lon to radians
    obs['lat'] = np.deg2rad(obs['lat'])
    obs['lon'] = np.deg2rad(obs['lon'])

    # Observation types
    # Key: Column name in DataFrame
    # name: Variable name in DART
    # type: Variable type in DART
    ob_types = {'cld_amt': {'name': 'CEIL_CLOUD_AMOUNT', 'type': 777}}
    
    # Determine the number of observations
    # NaN can be used for missing observations
    nobs = len(obs)
    for v in ob_types:
        nobs = nobs - np.sum(np.isnan(obs[v]))
    
    # Write obs_seq file
    with open(out_fname, 'w') as fptr:
        fptr.write(" obs_sequence\n")
        fptr.write("obs_kind_definitions\n")

        fptr.write(f"    {len(ob_types)} \n")
        for v in ob_types:
            fptr.write(f"    {ob_types[v]['type']}          {ob_types[v]['name']}   \n")
    
        fptr.write("  num_copies:            1  num_qc:            1\n")
        fptr.write(f" num_obs:       {nobs}  max_num_obs:       {nobs}\n")
        fptr.write("observation\n")
        fptr.write("QC value\n")
        fptr.write(f"  first:            1  last:       {nobs}\n")
        
        n = 0
        qc_val = 1
        for v in ob_types:
            kind = ob_types[v]['type']
            for i in range(len(obs)):
                fptr.write(f" OBS            {n+1}\n")
                fptr.write(f"   {obs['cld_amt'][i]:20.14f}\n")
                fptr.write(f"   {qc_val:20.14f}\n")

                if nobs == 1:
                    fptr.write(" -1 -1 -1\n") # Only 1 ob
                elif n+1 == 1: 
                    fptr.write(f" -1 {n+2} -1\n") # First ob
                elif n+1 == nobs:
                    fptr.write(f" {n} -1 -1\n") # Last ob
                else:
                    fptr.write(f" {n} {n+2} -1\n") 
            
                fptr.write("obdef\n")
                fptr.write("loc3d\n")
            
                fptr.write(f"    {obs['lon'][i]:20.14f}          {obs['lat'][i]:20.14f}          {obs['cld_hgt'][i]:20.14f}     3\n")
                fptr.write("kind\n")       
                fptr.write(f"     {kind}     \n")       
                fptr.write(f"    {time_sec}          {time_days}     \n")               
                fptr.write(f"    {obs['error'][i]:20.14f}  \n")
                
                n = n + 1
                
    return None


if __name__ == '__main__':

    start = dt.datetime.now()
    print('Starting create_syn_ceilometers.py')
    print(f"Time = {start.strftime('%Y%m%d %H:%M:%S')}")

    # Parse command-line arguments
    in_param = parse_in_args(sys.argv[1:])

    # Read in input files
    ceil_sites = pd.read_csv(in_param.ceil_file)
    NR_data = grib2io.open(in_param.NR_file, save_index=False)

    # Create ceilometer observations
    ceil_obs = create_ceil_obs(NR_data, ceil_sites, 
                               min_rng=in_param.min_range, 
                               max_rng=in_param.max_range, 
                               min_thres=in_param.min_cldfrac_thres)

    # Convert to oktas and add clear obs
    ceil_obs = convert_to_oktas(ceil_obs)
    ceil_obs = add_clear_obs_hgts(ceil_obs, 
                                  clr_ob_interval=in_param.clear_ob_interval,
                                  min_rng=in_param.min_range,
                                  max_rng=in_param.max_range)

    # Write out to DART obs_seq file
    write_ceil_to_obs_seq(ceil_obs, in_param.out_file, NR_data.read(1).validDate)

    print('Program Finished!')
    print(f"Elapsed time = {(dt.datetime.now() - start).total_seconds()} s")


"""
End create_syn_ceilometers.py
"""
