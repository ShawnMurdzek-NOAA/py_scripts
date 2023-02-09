"""
Plot UPP Output Fields

shawn.s.murdzek@noaa.gov
Date Created: 27 January 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
import plot_model_data as pmd
import xarray as xr


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

upp_files = []
save_fnames = []
for hr in range(12, 18):
    upp_files.append('/scratch1/BMC/wrfruc/murdzek/nature_run_spring_hrrr/output/20220429%d00/UPP/wrfnat_20220429%d00.grib2' % (hr, hr+1))
    save_fnames.append('/scratch1/BMC/wrfruc/murdzek/nature_run_spring_hrrr/other_plots/20220429%d00/upp_APCP_20220429%d00.png' % (hr, hr+1))

upp_files = ['/scratch1/BMC/wrfruc/murdzek/nature_run_tests/nature_run_spring_esstem/output/202204291300/UPP/wrfnat_202204291400.grib2']
save_fnames = ['/scratch1/BMC/wrfruc/murdzek/src/py_scripts/upp/rain_sum.png']

field = 'RWMR_P0_L105_GLC0'

name = 'NR: esstem'


#---------------------------------------------------------------------------------------------------
# Plot Data
#---------------------------------------------------------------------------------------------------

# Sample colors from plasma colorbar and define plotting levels
plevels = np.array([0.01, 0.1, 0.5, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]) * 1e-3
nlvl = len(plevels)
colors = [mcm.plasma(i / (nlvl-1)) for i in range(nlvl)]

for fname, save_fname in zip(upp_files, save_fnames):
    print('Plotting for %s' % fname)

    fig = plt.figure(figsize=(8, 4))
    plt.subplots_adjust(left=0.02, bottom=0.04, right=0.98, top=0.97)
    ds = xr.open_dataset(fname, engine='pynio')
    out = pmd.PlotOutput([ds], 'upp', fig, 1, 1, 1)
    out.contourf(field, ingest_kw={'red_fct':np.sum}, cntf_kw={'cmap':'plasma', 'extend':'max', 
                                                               'levels':np.arange(1e-3, 60e-3, 3e-3)})
    out.config_ax(grid=False)
    out.ax_title(txt=name, size=14)

    plt.savefig(save_fname)
    plt.close()


"""
End plot_upp.py
""" 
