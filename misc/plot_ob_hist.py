"""
Plot Histograms of Observations Values for Multiple BUFR Files

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

from pyDA_utils import bufr


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# BUFR file names
bufr_fnames = {'perfect':'/work2/noaa/wrfruc/murdzek/nature_run_spring/obs/uas_obs_150km/perf_uas_csv/202205042100.rap.fake.prepbufr.csv',
               'error':'/work2/noaa/wrfruc/murdzek/nature_run_spring/obs/uas_obs_150km/err_uas_csv/202205042100.rap.fake.prepbufr.csv',
               'superob':'/work2/noaa/wrfruc/murdzek/nature_run_spring/obs/uas_obs_150km/superob_uas/202205042100.rap.fake.prepbufr.csv'}

# Observation type(s) to plot
ob_types = [136]

# Variable to plot
ob_var = 'RHOB'

# Histogram bin edges
bins = np.arange(-0.5, 101, 1)

# X limit to plot
xlim = [92, 100.5]

# Output file name
save_fname = 'UAS_RH_hist.png'


#---------------------------------------------------------------------------------------------------
# Create Histograms
#---------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
ax.set_xlabel(ob_var, size=14)
ax.set_ylabel('density', size=14)

for key, c in zip(bufr_fnames.keys(), ['k', 'b', 'r', 'g']):
    print(f'Plotting histogram for {key}')
    bufr_obj = bufr.bufrCSV(bufr_fnames[key])
    bufr_obj.select_obtypes(ob_types)

    # Compute RH
    if ob_var == 'RHOB':
        bufr_obj.df = bufr.compute_RH(bufr_obj.df)

    bufr_obj.ob_hist(ob_var, ax=ax, hist_kw={'density':True, 'bins':bins}, plot_kw={'color':c, 'label':key})

ax.grid()
ax.legend()
ax.set_xlim(xlim)
ax.set_ylim(bottom=0)
plt.savefig(save_fname)
plt.show()


"""
End plot_ob_hist.py 
"""
