"""
Plot Multiple GSI Errtable Profiles on a Single Plot

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

etable = {'raob':{'fname':'/work2/noaa/wrfruc/murdzek/real_obs/errtables/errtable.rrfs',
                  'thermo':120,
                  'wind':220,
                  'scale':1,
                  'c':'b',
                  'ls':'-'},
          'aircft':{'fname':'/work2/noaa/wrfruc/murdzek/real_obs/errtables/errtable.rrfs',
                    'thermo':133,
                    'wind':233,
                    'scale':1,
                    'c':'r',
                    'ls':'-'},
          'UAS added err':{'fname':'/work2/noaa/wrfruc/murdzek/real_obs/errtables/2nd_iter_assim_only/include_uas/errtable.7day',
                    'thermo':136,
                    'wind':236,
                    'scale':1,
                    'c':'k',
                    'ls':'-'},
          'UAS added err X 2':{'fname':'/work2/noaa/wrfruc/murdzek/real_obs/errtables/2nd_iter_assim_only/include_uas/errtable.7day',
                    'thermo':136,
                    'wind':236,
                    'scale':2,
                    'c':'k',
                    'ls':'--'}}

out_fname = './uas_errtables.png'


#---------------------------------------------------------------------------------------------------
# Make Plots
#---------------------------------------------------------------------------------------------------

print('Reading errtables')
err_df = {}
for key in etable.keys():
    err_df[key] = gsi.read_errtable(etable[key]['fname'])    

print('Making plots...')

errs = {'Terr':'K', 'RHerr':'% /10', 'UVerr':'m s$^{-1}$', 'PSerr':'mb', 'PWerr':'mm'}
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12, 8), sharey=True)
plt.subplots_adjust(left=0.07, bottom=0.08, right=0.98, top=0.92)

for i, e in enumerate(errs.keys()):
    print(e)
    ax = axes[int(i/3), i%3]
    if e == 'UVerr':
        typ = 'wind'
    else:
        typ = 'thermo'
    for key in etable.keys():
        ax.plot(err_df[key][etable[key][typ]][e] * etable[key]['scale'], 
                err_df[key][etable[key][typ]]['prs'], 
                c=etable[key]['c'], ls=etable[key]['ls'],
                label=key, lw=2)

    ax.grid()
    ax.set_yscale('log')
    ax.set_ylim([1100, 600])
    ax.set_xlabel('{e} ({u})'.format(e=e, u=errs[e]), size=14)
    ax.set_xlim(left=0)

axes[1, 0].legend(fontsize=12)
for i in range(2):
    axes[i, 0].set_ylabel('pressure (mb)', size=14)        

plt.savefig(out_fname)


"""
End compare_errtable_profiles.py
"""
