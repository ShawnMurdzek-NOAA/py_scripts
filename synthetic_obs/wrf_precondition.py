"""
Precondition WRF NetCDF Output Files Prior to Synthetic Observation Generation

Preconditioning consists of paring down the number of fields and computing derived fields.

shawn.s.murdzek@noaa.gov
Date Created: 16 November 2022
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import netCDF4 as nc
import wrf


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

files = ['/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_10km_io102/WRF/run/wrfout_d01_2021-07-24_23_00_00']
save_files = ['/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_10km_io102/WRF/run/wrfred_2021-07-24_23_00_00']

# Fields to compute using wrf.getvar
compute_fields = ['height', 'height_agl', 'zstag', 'temp', 'td', 'rh', 'pw', 'pressure']

# Fields to save into the preconditioned netCDF file (in addition to compute_fields)
save_fields = ['XLAT', 'XLONG', 'U', 'V', 'W', 'QVAPOR', 'QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 
               'QGRAUP', 'QHAIL', 'REFL_10CM', 'XLAT_U', 'XLONG_U', 'XLAT_V', 'XLONG_V', 'Q2', 'T2', 
               'PSFC', 'U10', 'V10']

# Method ('serial' or 'together'). 'serial' requires ~20%  less memory than 'together'
method = 'serial'


#---------------------------------------------------------------------------------------------------
# Perform Preconditioning
#---------------------------------------------------------------------------------------------------

for f, sf in zip(files, save_files):

    fptr = nc.Dataset(f)

    if method == 'together':

        out_dict = {}
        for field in compute_fields+save_fields:
            out_dict[field] = wrf.getvar(fptr, field).drop_vars('Time')
            out_dict[field].attrs['projection'] = 'LambertConformal'

        new_ds = xr.Dataset(data_vars=out_dict)
        new_ds.to_netcdf(sf)

    elif method == 'serial':

        for j, field in enumerate(compute_fields+save_fields):
            out = wrf.getvar(fptr, field).drop_vars('Time')
            out.attrs['projection'] = 'LambertConformal'
            ds = xr.Dataset(data_vars={field:out})

            if j == 0:
                ds.to_netcdf(sf)
            else:
                ds.to_netcdf(sf, mode='a')
            ds.close()

    fptr.close()


"""
End wrf_precondition.py 
"""
