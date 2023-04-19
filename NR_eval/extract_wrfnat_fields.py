"""
Extract Single Fields from UPP Output Files

This simple script extracts a single 2D field from a bunch of UPP wrfnat files and saves them into
a single, separate netcdf file for easy use. UPP wrfnat files can be from the Nature Run or HRRR.

Command line arguments:
    argv[1] = Date (YYYYMMDD)

shawn.s.murdzek@noaa.gov
Date Created: 7 March 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import datetime as dt
import sys
import os


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

path = '/work2/noaa/wrfruc/murdzek/nature_run_winter/UPP'
date = sys.argv[1]
#date = '20220506'

# Input UPP wrfnat files (file names MUST end in .grib2, or else Xarray will throw an error)

# NR
template = '%s/%s/wrfnat_%s' % (path, date, date) + '%04d.grib2'
upp_files = [template % i for i in range(0, 1801, 600)]

# HRRR
#t1 = dt.datetime.strptime(date, '%Y%m%d')
#t0 = t1 - dt.timedelta(days=1)
#upp_files =  ['%s/%s/%s23000001.grib2' % (path, t0.strftime('%Y%m%d'), t0.strftime('%y%j'))]
#template = '%s/%s/%s' % (path, date, t1.strftime('%y%j')) + '%02d000001.grib2'
#upp_files = upp_files + [template % i for i in range(5, 19, 6)]

# Output netcdf file
out_file = '%s/%s/cref_%s.nc' % (path, date, date)
#out_file = '%s/%s/precip1hr_%s.nc' % (path, date, date)

# Field to extract from each UPP file
field = 'REFC_P0_L200_GLC0'  # Used in NR
#field = 'REFC_P0_L10_GLC0'   # Used in HRRR
#field = 'APCP_P8_L1_GLC0_acc'


#---------------------------------------------------------------------------------------------------
# Extract Fields and Combine
#---------------------------------------------------------------------------------------------------

extracted_field = []
times = []
for f in upp_files:
    print('extracting %s' % f)
    ds = xr.open_dataset(f, engine='pynio')
    extracted_field.append(ds[field].copy())
    times.append(dt.datetime.strptime(ds[field].attrs['initial_time'], '%m/%d/%Y (%H:%M)') +
                 dt.timedelta(hours=float(ds[field].attrs['forecast_time'][0])))
    ds.close()

final_da = xr.concat(extracted_field, dim='time')
final_da = final_da.assign_coords(coords={'time':times})

# Delete some metadata entries that are now useless
del final_da.attrs['initial_time']
del final_da.attrs['forecast_time']
del final_da.attrs['forecast_time_units']

final_da.to_netcdf(out_file)


"""
End extract_wrfnat_fields.py
"""
