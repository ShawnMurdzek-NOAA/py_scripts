"""
Compute Root-Mean-Squared Differences Between Two Sets of UPP Output

It is assumed that both sets of UPP output are on the same grid

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import xarray as xr
import datetime as dt


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

rrfs1_dir = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/real_red_data/winter_updated/NCO_dirs/ptmp/prod/'
#rrfs2_dir = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/real_red_data/winter_no_sfc/NCO_dirs/ptmp/prod/'
rrfs2_dir = '/work2/noaa/wrfruc/murdzek/RRFS_ORION/NCO_dirs/ptmp/prod/'

# UPP files
init_times = [dt.datetime(2022, 2, 1, 9) + dt.timedelta(hours=i) for i in range(3)]
upp_set1 = ['%s/rrfs.%s/%s/rrfs.t%sz.natlev.f003.conus_3km.grib2' % 
            (rrfs1_dir, t.strftime('%Y%m%d'), t.strftime('%H'), t.strftime('%H')) for t in init_times]
upp_set2 = ['%s/rrfs.%s/%s/rrfs.t%sz.natlev.f003.conus_3km.grib2' % 
            (rrfs2_dir, t.strftime('%Y%m%d'), t.strftime('%H'), t.strftime('%H')) for t in init_times]
labels = [t.strftime('%Y%m%d %H UTC') for t in init_times]

# Variables
upp_vars = ['TMP_P0_L105_GLC0', 'SPFH_P0_L105_GLC0', 'UGRD_P0_L105_GLC0', 'VGRD_P0_L105_GLC0',
            'HGT_P0_L105_GLC0']

# Output text file
out_text = 'RRFS_orion_update_diff_fhr003_UPDATED.txt'


#---------------------------------------------------------------------------------------------------
# Compute RMSDs
#---------------------------------------------------------------------------------------------------

# Open output file
fptr = open(out_text, 'w')

for fname1, fname2, l in zip(upp_set1, upp_set2, labels):
    print()
    print('opening %s' % l)
    ds1 = xr.open_dataset(fname1, engine='pynio')
    ds2 = xr.open_dataset(fname2, engine='pynio')
    fptr.write('\n%s:\n' % l)    

    for var in upp_vars:
        print('computing RMSD for %s' % var)
        rmsd = np.sqrt(np.mean((ds1[var].values - ds2[var].values)**2))
        fptr.write('%s = %.8f\n' % (var, rmsd))
        print(rmsd)

fptr.close()


"""
End model_rmsd.py
"""
