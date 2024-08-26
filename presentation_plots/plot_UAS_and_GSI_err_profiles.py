"""
Plot Vertical Error Profiles for UAS and GSI Errtable

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

import pyDA_utils.gsi_fcts as gsi


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# UAS errors
prs = np.array([1100, 100])
uas_errors = {'T':{'err':np.array(2*[0.5]),
                   'prs':prs},
              'RH':{'err':np.array(2*[5]),
                    'prs':prs},
              'UV':{'err':np.array(2*[0.6]),
                    'prs':prs}}
uas_label = 'original'

# GSI errtable file and ob ID to use
gsi_errtable_fname = '/work2/noaa/wrfruc/murdzek/real_obs/errtables/errtable.rrfs'
gsi_thermo_id = 133
gsi_wind_id = 233
gsi_label = 'big errors'

out_fname = 'uas_typ133_errtable.png'


#---------------------------------------------------------------------------------------------------
# Make Plot
#---------------------------------------------------------------------------------------------------

def format_gsi_errtable(errtable_fname, thermo_id, wind_id):
    """
    Read in and format GSI errtable file for plotting
    """

    errtable = gsi.read_errtable(errtable_fname)

    gsi_errors = {'T':{'err':errtable[thermo_id]['Terr'].values,
                       'prs':errtable[thermo_id]['prs'].values},
                  'RH':{'err':10*errtable[thermo_id]['RHerr'].values,
                        'prs':errtable[thermo_id]['prs'].values},
                  'UV':{'err':errtable[wind_id]['UVerr'].values,
                        'prs':errtable[wind_id]['prs'].values}}

    return gsi_errors


def plot_errtables(errors, figsize=(10, 5), labelsize=16, ylim=[1000, 600], subplot_adjust_kw={}):
    """
    Plot error table vertical profiles
    """

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=figsize, sharey=True)
    units = {'T':'K', 'RH':'%', 'UV':'m s$^{-1}$'}

    # Plot data
    for v, ax in zip(units.keys(), axes):
        for name in errors.keys():
            ax.plot(errors[name]['err'][v]['err'], errors[name]['err'][v]['prs'], label=name,
                    **errors[name]['plot_kw'])
        ax.set_xlabel(f"{v} ({units[v]})", size=labelsize)
        ax.set_xlim(left=0)
        ax.set_ylim(ylim)
        ax.set_yscale('log')
        ax.grid()
        ax.tick_params(axis='both', labelsize=(labelsize - 2))

    axes[0].set_ylabel("pressure (hPa)", size=labelsize)
    axes[-1].legend(fontsize=labelsize, loc='upper center')
    plt.subplots_adjust(**subplot_adjust_kw)

    return fig


if __name__ == '__main__':

    # Read in GSI errors
    gsi_errors = format_gsi_errtable(gsi_errtable_fname, gsi_thermo_id, gsi_wind_id)

    # Create plot
    errors = {uas_label:{'err':uas_errors,
                         'plot_kw':{'color':'r',
                                     'lw':4}},
              gsi_label:{'err':gsi_errors,
                         'plot_kw':{'color':'saddlebrown',
                                    'lw':4}}}
    fig = plot_errtables(errors, subplot_adjust_kw={'left':0.1, 'bottom':0.12, 'right':0.98, 'top':0.98, 'wspace':0.1})
    plt.savefig(out_fname)


"""
End plot_UAS_and_GSI_err_profiles.py
"""
