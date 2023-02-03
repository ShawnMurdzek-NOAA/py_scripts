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

import plot_model_data as pmd


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# UPP output files
upp1_files = []
upp2_files = []
home = '/scratch1/BMC/wrfruc/murdzek'
for hr in range(12, 14):
    upp1_files.append('%s/nature_run_spring_mynn/output/20220429%d00/UPP/wrfnat_20220429%d00.grib2' % (home, hr, hr+1))
    upp2_files.append('%s/nature_run_spring_v2/output/20220429%d00/UPP/wrfnat_20220429%d00.grib2' % (home, hr, hr+1))

# Field and level to plot (set level to np.nan for a 2D field)
#field = 'REFC_P0_L200_GLC0'
field = 'UGRD_P0_L103_GLC0'
zind = 0

# Colorbar and contour levels
cmap = 'plasma'
lvls = np.arange(-20, 20, 2)
extend = 'both'

# Titles
name1 = 'NR: MYNN'
name2 = 'NR: Shin-Hong'

# Output file name (must be a PDF)
save_fname = './mynn_ShinHong_U10_diff.pdf'


#---------------------------------------------------------------------------------------------------
# Plot Data
#---------------------------------------------------------------------------------------------------

# Create PDF
pdf = PdfPages(save_fname)

# Create plots
for fname1, fname2 in zip(upp1_files, upp2_files):
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
                                               'levels':np.linspace(-5, 5, 20)},
                                      ingest_kw={'zind':zind})
    diff.config_ax(grid=False)
    diff.ax_title(txt='difference (top - bottom)', size=14)

    plt.subplots_adjust(left=0.02, bottom=0.04, right=0.98, top=0.97, hspace=0.1)
    pdf.savefig(fig)
    plt.close()

pdf.close()


"""
End plot_upp_diffs.py
""" 
