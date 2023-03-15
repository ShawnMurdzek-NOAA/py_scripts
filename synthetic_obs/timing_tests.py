
import xarray as xr
import datetime as dt
import os

ds = xr.open_dataset('/mnt/lfs4/BMC/wrfruc/murdzek/nature_run_spring/UPP/20220430/wrfnat_202204301300.grib2', engine='pynio')

x = ds['PRES_P0_L105_GLC0'].values
for l in os.popen('free -t -m -h').readlines():
    print(l)

print()
print('size of x in bytes = %d' % x.nbytes)

'''
f2d = 'PRES_P0_L1_GLC0'
f3d = 'PRES_P0_L105_GLC0'
for iend in [100, 400]:
    for jend in [100, 400]:
        print()
        print('iend = %d, jend = %d' % (iend, jend))
        t1 = dt.datetime.now()
        x2d = ds[f2d][2:iend, 2:jend].values
        t2 = dt.datetime.now()
        x3d = ds[f3d][0, 2:iend, 2:jend].values
        t3 = dt.datetime.now()
        y3d = ds[f3d][:, 2:iend, 2:jend].values
        t4 = dt.datetime.now()
        print('time to extract 2D field = %.2f s' % (t2 - t1).total_seconds())
        print('time to extract 3D field (1 z level) = %.2f s' % (t3 - t2).total_seconds())
        print('time to extract 3D field (all z levels) = %.2f s' % (t4 - t3).total_seconds())
'''
