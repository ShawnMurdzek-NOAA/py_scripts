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


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

upp_files = []
save_fnames = []
for hr in range(12, 19):
    upp_files.append('/scratch1/BMC/wrfruc/murdzek/nature_run_spring/output/20220429%d00/UPP/wrfnat_20220429%d00.grib2' % (hr, hr+1))
    save_fnames.append('/scratch1/BMC/wrfruc/murdzek/nature_run_spring/other_plots/20220429%d00/upp_APCP_20220429%d00.png' % (hr, hr+1))

field = 'APCP_P8_L1_GLC0_acc'

name = 'Nature Run'


#---------------------------------------------------------------------------------------------------
# Plot Data
#---------------------------------------------------------------------------------------------------

# Sample colors from plasma colorbar and define plotting levels
plevels = np.array([0.002, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 2.0]) * 0.0254 * 997
nlvl = len(plevels)
colors = [mcm.plasma(i / (nlvl-1)) for i in range(nlvl)]

for fname, save_fname in zip(upp_files, save_fnames):
    print('Plotting for %s' % fname)

    fig = plt.figure(figsize=(8, 4))
    plt.subplots_adjust(left=0.02, bottom=0.04, right=0.98, top=0.97)
    out = pmd.PlotOutput(fname, 'upp', fig, 1, 1, 1)
    out.contourf(field, cntf_kw={'colors':colors, 'extend':'max', 'levels':plevels})
    out.config_ax(grid=False)
    out.ax_title(txt=name, size=14)

    plt.savefig(save_fname)
    plt.close()


"""
End plot_upp.py
""" 
