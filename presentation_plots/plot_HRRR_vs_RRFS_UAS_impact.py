"""
Plot UAS Obs Impacts for both the HRRR and RRFS

Used for AMS 2025

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import datetime as dt
import matplotlib.pyplot as plt
import copy

import metplus_OSSE_scripts.plotting.metplus_plots as mp
import metplus_OSSE_scripts.plotting.metplus_tools as mt


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Input simulations
RRFS_sims = {'ctrl':{'color':'k',
                     'dir':'/work2/noaa/wrfruc/murdzek/RRFS_OSSE/metplus_verif_grid_NR/rrfs-workflow_orion/spring/lower_atm/output/GridStat',
                     'ctrl':True},
             'diff_uas150':{'color':'b',
                               'dir':'/work2/noaa/wrfruc/murdzek/RRFS_OSSE/metplus_verif_grid_NR/rrfs-workflow_orion/spring_uas_150km/lower_atm/output/GridStat',
                               'ctrl':False},
             'diff_uas35':{'color':'r',
                              'dir':'/work2/noaa/wrfruc/murdzek/RRFS_OSSE/metplus_verif_grid_NR/rrfs-workflow_orion/spring_uas_35km/lower_atm/output/GridStat',
                               'ctrl':False}}
HRRR_sims = {'ctrl':{'color':'k',
                     'dir':'/work2/noaa/wrfruc/murdzek/HRRR_OSSE/metplus_verif_grid_NR/jet_hrrr/spring/lower_atm/output/GridStat',
                     'ctrl':True},
             'diff_uas150':{'color':'b',
                               'dir':'/work2/noaa/wrfruc/murdzek/HRRR_OSSE/metplus_verif_grid_NR/jet_hrrr/spring_uas_150km/lower_atm/output/GridStat',
                               'ctrl':False},
             'diff_uas35':{'color':'r',
                              'dir':'/work2/noaa/wrfruc/murdzek/HRRR_OSSE/metplus_verif_grid_NR/jet_hrrr/spring_uas_35km/lower_atm/output/GridStat',
                               'ctrl':False}}

# Valid times
valid_times = [dt.datetime(2022, 4, 30, 0)  + dt.timedelta(hours=i) for i in range(0, 156, 1)]
valid_times.remove(dt.datetime(2022, 4, 30, 1))
valid_times.remove(dt.datetime(2022, 5, 6, 8))

# Keyword arguments passed to mp.plot_ua_vprof()
plot_kw = {'fcst_lead':3,
           'file_prefix':'grid_stat_FV3_TMP_vs_NR_TMP',
           'toggle_pts':True,
           'exclude_plvl':[],
           'verbose':False,
           'ylim':[1015, 500],
           'diffs':True,
           'include_ctrl':True,
           'diff_kw':{'match':['FCST_LEAD',
                               'FCST_VAR',
                               'FCST_VALID_BEG',
                               'FCST_LEV',
                               'FCST_UNITS',
                               'VX_MASK',
                               'OBTYPE']},
           'ci':True,
           'ci_lvl':0.95,
           'ci_opt':'t_dist',
           'ci_kw':{'acct_lag_corr':True}}

# X-axis limits
xlim = {'TMP':[-0.4, 1.45],
        'SPFH':[-0.00048, 0.00125],
        'UGRD_VGRD':[-1.2, 4]}

# Output file name
save_fname = 'HRRR_vs_RRFS_UAS_impact_lower_atm.png'


#---------------------------------------------------------------------------------------------------
# Make Plots
#---------------------------------------------------------------------------------------------------

# Create plot
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(13, 10), sharey=True, sharex='col')
plt.subplots_adjust(left=0.11, bottom=0.13, right=0.98, top=0.93, hspace=0.05, wspace=0.05)

# Loop over each variable
for i, (var, lt, stat) in enumerate(zip(['TMP', 'SPFH', 'UGRD_VGRD'], 
                                        ['sl1l2', 'sl1l2', 'vl1l2'],
                                        ['RMSE', 'RMSE', 'VECT_RMSE'])):
    print(f'plotting {var}')
    for j, sims in enumerate([RRFS_sims, HRRR_sims]):
        plot_kw_copy = copy.deepcopy(plot_kw)
        plot_kw_copy['diff_kw']['var'] = [stat]
        _ = mp.plot_ua_vprof(sims, valid_times, ax=axes[j, i], line_type=lt, plot_stat=stat, 
                             plot_param={'FCST_VAR':var, 'OBTYPE':'NR', 'VX_MASK':'shape_mask'}, 
                             **plot_kw_copy)
        axes[j, i].set_xlim(xlim[var])
        axes[j, i].set_title('')
        axes[j, i].tick_params(axis='both', which='both', labelsize=14)

# Add y-axes labels
for i, label in enumerate(['RRFS', 'HRRR-like']):
    axes[i, 0].set_ylabel(f'{label}\npressure (hPa)', size=18)
    for j in range(1, 3):
        axes[i, j].set_ylabel('')

for i, label in enumerate(['temperature RMSE (K)', 'specific humidity RMSE (kg kg$^{-1}$)', 'wind RMSE (m s$^{-1}$)']):
    axes[0, i].set_xlabel('')
    axes[1, i].set_xlabel(label, size=18)

plt.suptitle('3-hr Forecast', size=22)

plt.savefig(save_fname)


"""
End plot_HRRR_vs_RRFS_UAS_impact.py
"""
