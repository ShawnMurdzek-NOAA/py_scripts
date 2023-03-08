"""
Extract Single Fields from UPP Output Files

This simple script extracts a single 2D field from a bunch of UPP wrfnat files and saves them into
a single, separate netcdf file for easy use.

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


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

path = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_spring/UPP'
#date = sys.argv[1]
date = '20220430'

# Input wrfnat files
template = '%s/%s/wrfnat_%s' % (path, date, date) + '%04d.grib2'
upp_files =  [template % i for i in range(0, 1801, 600)]

# Output netcdf file
out_file = '%s/%s/cref_%s.nc' % (path, date, date)
#out_file = '%s/%s/precip1hr_%s.nc' % (path, date, date)

# Field to extract from each UPP file
field = 'REFC_P0_L200_GLC0'
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
