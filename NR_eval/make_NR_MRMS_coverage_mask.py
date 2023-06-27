"""
Create Masks for NR and MRMS Output

Both masks hide gridpoints that occur in regions with no MRMS coverage (a nearest neighbor
interpolation scheme is used to determine which NR gridpoints lie in MRMS regions with no coverage)
and gridpoints that lie outside the domain of the other produce (i.e., only gridpoints within the
common domain are retained).

This script requires quite a bit of memory. I usually run it with 20 GB of RAM

Timing:
    - Creating the NR mask takes 7-8 min
    - Creating the MRMS mask takes

shawn.s.murdzek@noaa.gov
Date Created: 13 March 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import numpy as np
import scipy.interpolate as si
import datetime as dt

import pyDA_utils.map_proj as mp


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

NR_file = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_spring/UPP/20220429/wrfnat_202204291300.grib2'
#NR_file = '/mnt/lfs4/BMC/wrfruc/murdzek/HRRR_data/20220429/2211917000001.grib2'
MRMS_file = '/mnt/lfs4/BMC/wrfruc/murdzek/real_obs/mrms/2015/20150429-180012.MRMS_MergedReflectivityQCComposite_00.50_20150429-180012.grib2'
NR_mask_file = './NR_mask.npy'
#NR_mask_file = './HRRR_mask.npy'
MRMS_mask_file = './MRMS_mask.npy'

MRMS_no_coverage = -999.
MRMS_field = 'MergedReflectivityQCComposite_P0_L102_GLL0'


#---------------------------------------------------------------------------------------------------
# Create Masks
#---------------------------------------------------------------------------------------------------

NR_ds = xr.open_dataset(NR_file, engine='pynio')
MRMS_ds = xr.open_dataset(MRMS_file, engine='pynio')

print('creating interpolation fct (time = %s)' % dt.datetime.now().strftime('%H:%M:%S'))
MRMS_lon, MRMS_lat = np.meshgrid(MRMS_ds['lon_0'].values - 360., MRMS_ds['lat_0'].values)
MRMS_coverage = np.int16(MRMS_ds[MRMS_field].values > MRMS_no_coverage)
interp_fct = si.NearestNDInterpolator(list(zip(MRMS_lon.ravel(), MRMS_lat.ravel())), 
                                      MRMS_coverage.ravel())

print('performing interpolation (time = %s)' % dt.datetime.now().strftime('%H:%M:%S'))
NR_lon = NR_ds['gridlon_0'].values
NR_lat = NR_ds['gridlat_0'].values
NR_coverage = interp_fct(NR_lon, NR_lat)

print("masking gridpoints that don't lie in both domains (time = %s)" % 
      dt.datetime.now().strftime('%H:%M:%S'))
MRMS_x, MRMS_y = mp.ll_to_xy_lc(MRMS_lat, MRMS_lon)
NR_ny, NR_nx = NR_ds['gridlat_0'].shape
MRMS_mask = (MRMS_coverage * np.logical_and(MRMS_x >= 0, MRMS_x <= NR_nx) * 
             np.logical_and(MRMS_y >= 0, MRMS_y <= NR_ny))
NR_mask = (NR_coverage * np.logical_and(NR_lon >= MRMS_lon[0, 0], NR_lon <= MRMS_lon[-1, -1]) * 
           np.logical_and(NR_lat >= MRMS_lat[-1, -1], NR_lat <= MRMS_lat[0, 0]))

print('saving NR mask (time = %s)' % dt.datetime.now().strftime('%H:%M:%S'))
np.save(NR_mask_file, NR_mask)

print('saving MRMS mask (time = %s)' % dt.datetime.now().strftime('%H:%M:%S'))
np.save(MRMS_mask_file, MRMS_mask)


"""
End make_NR_MRMS_coverage_mask.py
"""
