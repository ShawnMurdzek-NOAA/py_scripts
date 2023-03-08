"""
Extract Single Fields from UPP Output Files

This simple script extracts a single 2D field from a bunch of UPP wrfnat files and saves them into
a single, separate netcdf file for easy use.

shawn.s.murdzek@noaa.gov
Date Created: 7 March 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import datetime as dt


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

path = '/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_spring/UPP'
date = '20220430'

# Input wrfnat files
template = '%s/%s/wrfnat_%s' % (path, date, date) + '%04d.grib2'
upp_files =  [template % i for i in range(0, 1801, 600)]

# Output netcdf file
out_file = '%s/%s/cref_%s.nc' % (path, date, date)

# Field to extract from each UPP file
field = 'REFC_P0_L200_GLC0'


#---------------------------------------------------------------------------------------------------
# Extract Fields and Combine
#---------------------------------------------------------------------------------------------------

extracted_field = []
times = []
for f in upp_files:
    ds = xr.open_dataset(f, engine='pynio')
    extracted_field.append(ds[field].copy())
    times.append(dt.datetime.strptime(ds[field].attrs['initial_time'], '%m/%d/%Y (%H:%M)') +
                 dt.timedelta(hours=ds[field].attrs['forecast_time'][0]))
    ds.close()

final_da = xr.concat(extracted_field, dim='time')
final_da.assign_coords(coords={'time':times})

final_da.to_netcdf(out_file)


"""
End extract_wrfnat_fields.py
"""
