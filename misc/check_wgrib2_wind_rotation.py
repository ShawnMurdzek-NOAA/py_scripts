"""
Check wgrib2 Wind Rotation

Compare winds rotated from grid-relative to earth-relative using wgrib2 (w/ the new_grid_winds 
option) and using the gridrot_0 field in UPP

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import subprocess


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

grid_rel_fname = '/work2/noaa/wrfruc/murdzek/nature_run_winter/UPP/20220202/tmp/wrfnat_202202021000.grib2'
earth_rel_fname = '/work2/noaa/wrfruc/murdzek/nature_run_winter/UPP/20220202/tmp/out.grib2'


#---------------------------------------------------------------------------------------------------
# Examine Differences
#---------------------------------------------------------------------------------------------------

grid_rel_ds = xr.open_dataset(grid_rel_fname, engine='pynio')
earth_rel_ds = xr.open_dataset(earth_rel_fname, engine='pynio')

# Rotate winds to be earth-relative in grid_rel_ds
print('rotating winds to be earth-relative')
gridrot = grid_rel_ds['gridrot_0'].values
grid_rel_ds['UEAR_P0_L103_GLC0'] = (np.sin(gridrot)*grid_rel_ds['VGRD_P0_L103_GLC0'] +
                                    np.cos(gridrot)*grid_rel_ds['UGRD_P0_L103_GLC0'])
grid_rel_ds['VEAR_P0_L103_GLC0'] = (np.cos(gridrot)*grid_rel_ds['VGRD_P0_L103_GLC0'] -
                                    np.sin(gridrot)*grid_rel_ds['UGRD_P0_L103_GLC0'])
earth_rel_ds['UEAR_P0_L103_GLC0'] = earth_rel_ds['UGRD_P0_L103_GLC0']
earth_rel_ds['VEAR_P0_L103_GLC0'] = earth_rel_ds['VGRD_P0_L103_GLC0']

# Compute some RMSDs
fields = ['TMP_P0_L105_GLC0', 'UGRD_P0_L105_GLC0', 
          'UGRD_P0_L103_GLC0', 'UEAR_P0_L103_GLC0',
          'VGRD_P0_L103_GLC0', 'VEAR_P0_L103_GLC0']
print()
for f in fields:
    maxdiff = np.amax(np.abs(grid_rel_ds[f] - earth_rel_ds[f]))
    rmsd = np.sqrt(np.mean((grid_rel_ds[f] - earth_rel_ds[f])**2))
    print('RMSD for {f} = {v:.3e} (max diff = {d:.3e})'.format(f=f, v=rmsd, d=maxdiff))




"""
End check_wgrib2_wind_rotation.py
"""
