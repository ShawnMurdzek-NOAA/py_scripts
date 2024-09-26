"""
Plot Difference Between UPP Output Fields

shawn.s.murdzek@noaa.gov
Date Created: 3 February 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
from matplotlib.backends.backend_pdf import PdfPages
import xarray as xr

import pyDA_utils.plot_model_data as pmd


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# UPP output files
upp1_files = []
upp2_files = []
home = '/work2/noaa/wrfruc/murdzek'
upp1_files = [f'{home}/RRFS_OSSE/syn_data_rrfs-workflow_orion/winter_uas_35km/NCO_dirs/ptmp/prod/rrfs.20220201/09/rrfs.t09z.natlev.f001.conus_3km.grib2']
upp2_files = [f'{home}/RRFS_orion_rocky9_test/NCO_dirs/ptmp/prod/rrfs.20220201/09/rrfs.t09z.natlev.f001.conus_3km.grib2']

# Field and level to plot (set level to np.nan for a 2D field)
#field = 'REFC_P0_L200_GLC0'
field = 'TMP_P0_L105_GLC0'
zind = [10, 10]

# Colorbar and contour levels
cmap = 'Reds'
lvls = np.arange(240, 315, 5)
extend = 'both'

# Titles
name1 = 'CentOS 7'
name2 = 'Rocky 9'

# Output file name (must include a %d placeholder)
save_fname = './orion_os_test_tmp10_diff_%d.png'


#---------------------------------------------------------------------------------------------------
# Plot Data
#---------------------------------------------------------------------------------------------------

# Create PDF
#pdf = PdfPages(save_fname)

# Create plots
for i, (fname1, fname2) in enumerate(zip(upp1_files, upp2_files)):
    print('Plotting for %s' % fname1)
    fig = plt.figure(figsize=(8, 11))
    ds = []
    for j, (f, n) in enumerate(zip([fname1, fname2], [name1, name2])):
        ds.append(xr.open_dataset(f, engine='pynio'))
        out = pmd.PlotOutput([ds[-1]], 'upp', fig, 3, 1, j+1)
        out.contourf(field, cntf_kw={'cmap':cmap, 'extend':extend, 'levels':lvls}, 
                     ingest_kw={'zind':zind})
        out.config_ax(grid=False)
        out.ax_title(txt=n, size=14)

    diff = pmd.PlotOutput(ds, 'upp', fig, 3, 1, 3)
    diff.plot_diff(field, auto=False, cntf_kw={'cmap':'bwr', 'extend':'both', 
                                               'levels':np.linspace(-2, 2, 20)},
                                      ingest_kw={'zind':zind})
    diff.config_ax(grid=False)
    diff.ax_title(txt='difference (top - bottom)', size=14)

    plt.subplots_adjust(left=0.02, bottom=0.04, right=0.98, top=0.97, hspace=0.1)
    #pdf.savefig(fig)
    plt.savefig(save_fname % i)
    plt.close()

#pdf.close()


"""
End plot_upp_diffs.py
""" 
