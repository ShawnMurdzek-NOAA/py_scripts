"""
Plot GSI Error Table Profiles

One plot is created for each observation type, with one subplot for each variable. Multiple error
tables can be plotted at once.

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.cm as cm

import pyDA_utils.gsi_fcts as gsi


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

#etable_std_fnames = (['/work2/noaa/wrfruc/murdzek/real_obs/errtables/errtable.rrfs'] + 
#                     ['/work2/noaa/wrfruc/murdzek/real_obs/errtables/1st_iter/errtable.1st_iter.%dday' % i for i in range(1, 8)])
etable_std_fnames = ['/work2/noaa/wrfruc/murdzek/real_obs/errtables/1st_iter_assim_only/errtable.1st_iter.2day',
                     '/work2/noaa/wrfruc/murdzek/real_obs/errtables/1st_iter_assim_only/errtable.1st_iter.7day']
etable_mean_fname = []

#etable_labels = ['original'] + ['%d day' % i for i in range(1, 8)]
#etable_colors = ['k'] + [cm.plasma(i / (len(etable_std_fnames) - 1)) for i in range(1, 8)]
#etable_linestyles = ['--'] + ['-'] * (len(etable_std_fnames) - 1)

etable_labels = ['2day', '7day']
etable_colors = ['r', 'b']
etable_linestyles = ['-', '-']

out_fname = './errtable_profiles.pdf'

# Select observation types to save as PNGs
ob_pngs = []


#---------------------------------------------------------------------------------------------------
# Make Plots
#---------------------------------------------------------------------------------------------------

print('Reading in error tables...')
err_df = {}
#for stat, l in zip(['std', 'mean'], [etable_std_fnames, etable_mean_fnames]):
for stat, l in zip(['std'], [etable_std_fnames]):
    err_df[stat] = {}
    for fname, label in zip(l, etable_labels):
        err_df[stat][label] = gsi.read_errtable(fname)    

print('Making plots...')
errs = ['Terr', 'RHerr', 'UVerr', 'PSerr', 'PWerr']
pdf = PdfPages(out_fname)
for typ in err_df['std'][etable_labels[0]].keys():

    print(typ)

    # Skip this observation if all errors are NaN
    skip = True
    for label in etable_labels:
        for e in errs:
            if not np.all(np.isnan(err_df['std'][label][typ][e])):
                skip = False
                break
    if skip:
        continue

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12, 8), sharey=True)
    plt.subplots_adjust(left=0.07, bottom=0.08, right=0.98, top=0.92)
    for j, e in enumerate(errs):
        ax = axes[int(j/3), j%3]
        for label, c, ls in zip(etable_labels, etable_colors, etable_linestyles):
            for stat in ['std']:
                ax.plot(err_df[stat][label][typ][e], err_df[stat][label][typ]['prs'], c=c, ls=ls,
                        label=label, lw=2)
        ax.legend()
        ax.grid()
        ax.set_yscale('log')
        if typ in (list(range(180, 200)) + list(range(280, 300))):
            ax.set_ylim([1100, 650])
        else:
            ax.set_ylim([1100, 10])
        ax.set_xlabel(e, size=14)
        ax.set_xlim(left=0)
    plt.suptitle('Type = %d' % typ, size=20)
    for j in range(2):
        axes[j, 0].set_ylabel('pressure (mb)', size=14)        

    if typ in ob_pngs:
        plt.savefig('ob%d_errtable_profiles.png' % typ)

    pdf.savefig(fig)
    plt.close(fig)
pdf.close()


"""
End plot_errtable_profiles.py
"""
