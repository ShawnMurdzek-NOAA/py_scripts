"""
Plot Stage IV Precipitation Estimates

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

path = '/scratch1/BMC/wrfruc/murdzek/stage4_precip/data'
start = dt.datetime(2022, 4, 29, 12)
end = dt.datetime(2022, 4, 29, 18)
step = 3600     # time between output files, in sec
acc_window = '01h'
domain = 'conus'

save_dir = '/scratch1/BMC/wrfruc/murdzek/stage4_precip/plots'


#---------------------------------------------------------------------------------------------------
# Plot Data
#---------------------------------------------------------------------------------------------------

# Sample colors from plasma colorbar and define plotting levels
plevels = np.array([0.002, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 2.0]) * 0.0254 * 997
nlvl = len(plevels)
colors = [mcm.plasma(i / (nlvl-1)) for i in range(nlvl)]

t = start
while t <= end:
    fname = '%s/st4_%s.%s.%s.grb2' % (path, domain, t.strftime('%Y%m%d%H'), acc_window)
    print('Plotting for %s' % fname)

    fig = plt.figure(figsize=(8, 4))
    plt.subplots_adjust(left=0.02, bottom=0.04, right=0.98, top=0.97)
    out = pmd.PlotOutput(fname, 'stage4', fig, 1, 1, 1)
    out.contourf('APCP_P8_L1_GST0_acc', cntf_kw={'colors':colors, 'extend':'max', 
                                                 'levels':plevels})
    out.config_ax(grid=False)
    out.ax_title(size=14)

    plt.savefig('%s/st4_%s_%s_%s.png' % (save_dir, domain, t.strftime('%Y%m%d%H'), acc_window))
    plt.close()

    t = t + dt.timedelta(seconds=step)


"""
End plot_stage4_precip.py
""" 
