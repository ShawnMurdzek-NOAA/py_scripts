"""
Plot GSI O-B Variances and Error Table Profiles

One plot is created for each observation type, with two subplots for each variable (one to show
O-B variances and errtable variances and one to show the number of obs assimilated).

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.cm as cm
import datetime as dt

import pyDA_utils.gsi_fcts as gsi


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# GSI diag files
tmpl_real = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/real_red_data/winter/NCO_dirs/ptmp/prod/rrfs.%s/%s/'
tmpl_osse = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data/winter_perfect/NCO_dirs/ptmp/prod/rrfs.%s/%s/'
dates = [dt.datetime(2022, 2, 1, 9) + dt.timedelta(hours=i) for i in range(159)]
gsi_inputs = {'real': {'dirs':[tmpl_real % (d.strftime('%Y%m%d'), d.strftime('%H')) for d in dates],
                       'c':'b'},
              'OSSE': {'dirs':[tmpl_osse % (d.strftime('%Y%m%d'), d.strftime('%H')) for d in dates],
                       'c':'r'}}

# Error tables (for variances, not means)
etable_inputs = {'7 day': {'fname':'/work2/noaa/wrfruc/murdzek/real_obs/errtables/1st_iter/errtable.1st_iter.7day',
                           'c':'c',
                           'ls':':'}}

# Option to plot the difference in the O-B variances for the first two GSI diag files
plot_diff = True

out_fname = './omb_errtable_profiles_1st_iter_7day.pdf'

# Select observation types to save as PNGs
ob_pngs = []


#---------------------------------------------------------------------------------------------------
# Make Plots
#---------------------------------------------------------------------------------------------------

print('Reading in error table(s)...')
err_df = {}
for key in etable_inputs.keys():
    err_df[key] = gsi.read_errtable(etable_inputs[key]['fname'])    

print()
print('----------------------------')
print('Reading in GSI diag files...')
omb_df = {}
ob_types = []
for key in gsi_inputs.keys():
    print(key)
    omb_df[key] = {}
    for v in ['t', 'q', 'ps', 'pw', 'uv']:
        print(v)
        tmp_fnames = ['%s/diag_conv_%s_%s.%s.nc4' % (path, v, 'ges', d.strftime('%Y%m%d%H'))
                      for path, d in zip(gsi_inputs[key]['dirs'], dates)]
        if v == 'uv':
            tot_df = gsi.read_diag(tmp_fnames)
            for v2 in ['u', 'v']:
                omb_df[key][v2] = tot_df.copy()
                omb_df[key][v2]['Obs_Minus_Forecast_adjusted'] = tot_df['%s_Obs_Minus_Forecast_adjusted' % v2].copy()
        elif v == 'q':
            # Convert q to RH for comparison with RHerr from error table
            omb_df[key][v] = gsi.read_diag(tmp_fnames)
            omb_df[key][v]['q_Obs_Minus_Forecast_adjusted'] = omb_df[key][v]['Obs_Minus_Forecast_adjusted'].copy()
            obs_q = omb_df[key][v]['Observation'].values
            background_q = obs_q - omb_df[key][v]['q_Obs_Minus_Forecast_adjusted'].values
            qs = omb_df[key][v]['Forecast_Saturation_Spec_Hum'].values
            omb_df[key][v]['Obs_Minus_Forecast_adjusted'] = 10 * ((obs_q / qs) - (background_q / qs))
        else:
            omb_df[key][v] = gsi.read_diag(tmp_fnames)

    # Create list of observation types for plotting
    for v in ['t', 'q', 'ps', 'pw', 'u', 'v']:
        for o in omb_df[key][v]['Observation_Type'].values:
            if o not in ob_types:
                ob_types.append(o)

ob_types = np.sort(np.array(ob_types))

print()
print('---------------')
print('Making plots...')
pdf = PdfPages(out_fname)
for typ in ob_types:

    print(typ)

    if typ < 200:
        ncols = 4
        fig, axes = plt.subplots(nrows=2, ncols=ncols, figsize=(12, 8), sharey=True)
        plt.subplots_adjust(left=0.07, bottom=0.08, right=0.98, top=0.92)
        gsi_vars = ['t', 'q', 'pw', 'ps']
        err_vars = ['Terr', 'RHerr', 'PWerr', 'PSerr']
    elif typ > 200:
        ncols = 2
        fig, axes = plt.subplots(nrows=2, ncols=ncols, figsize=(12, 8), sharey=True)
        plt.subplots_adjust(left=0.07, bottom=0.08, right=0.98, top=0.92)
        gsi_vars = ['u', 'v']
        err_vars = ['UVerr', 'UVerr']

    for j, (gv, ev) in enumerate(zip(gsi_vars, err_vars)):
        ax_omb = axes[int((2*j)/ncols), (2*j)%ncols]
        ax_cts = axes[int((2*j)/ncols), (2*j)%ncols+1]

        # Plot error tables variances
        for key in err_df.keys():
            ax_omb.plot(err_df[key][typ][ev]**2, err_df[key][typ]['prs'], c=etable_inputs[key]['c'],
                        ls=etable_inputs[key]['ls'], label='etable %s' % key)

        # Plot O-Bs and ob counts
        prs_ctr = err_df[list(err_df.keys())[0]][typ]['prs'].values
        prs = np.zeros(len(prs_ctr))
        prs[1:] = 0.5 * (prs_ctr[:-1] + prs_ctr[1:])
        prs[0] = prs_ctr[0] + 0.5 * (prs_ctr[1] - prs_ctr[2])
        for key in omb_df.keys():
            subset = omb_df[key][gv].loc[omb_df[key][gv]['Observation_Type'] == typ]
            omb_var = {'all':np.zeros(len(prs)-1), 'assim':np.zeros(len(prs)-1), 
                       'qc<3':np.zeros(len(prs)-1)}
            omb_cts = {'all':np.zeros(len(prs)-1), 'assim':np.zeros(len(prs)-1), 
                       'qc<3':np.zeros(len(prs)-1)}
            for ki in range(len(prs) - 1):
                tmp_subset = subset.loc[(subset['Pressure'] < prs[ki]) & 
                                        (subset['Pressure'] >= prs[ki+1])]
                omb_var['all'][ki] = np.var(tmp_subset['Obs_Minus_Forecast_adjusted'].values)
                omb_cts['all'][ki] = len(tmp_subset)

                tmp_subset = subset.loc[(subset['Pressure'] < prs[ki]) & 
                                        (subset['Pressure'] >= prs[ki+1]) &
                                        (subset['Analysis_Use_Flag'] == 1)]
                omb_var['assim'][ki] = np.var(tmp_subset['Obs_Minus_Forecast_adjusted'].values)
                omb_cts['assim'][ki] = len(tmp_subset)

                tmp_subset = subset.loc[(subset['Pressure'] < prs[ki]) & 
                                        (subset['Pressure'] >= prs[ki+1]) &
                                        (subset['Prep_QC_Mark'] < 3)]
                omb_var['qc<3'][ki] = np.var(tmp_subset['Obs_Minus_Forecast_adjusted'].values)
                omb_cts['qc<3'][ki] = len(tmp_subset)

            if plot_diff:
                if key == list(omb_df.keys())[0]:
                    diff = omb_var['qc<3']
                elif key == list(omb_df.keys())[1]:
                    diff = diff - omb_var['qc<3']
                    ax_omb.plot(diff, prs_ctr[:-1], c='gray', ls='--', 
                                label='diff (qc<3)')

            for label, ls in zip(['all', 'assim', 'qc<3'], ['-', '--', ':']):
                for ax, data in zip([ax_omb, ax_cts], [omb_var[label], omb_cts[label]]):
                    ax.plot(data, prs_ctr[:-1], c=gsi_inputs[key]['c'], ls=ls, 
                            label='%s %s' % (key, label))

        # Add labels and configure axes
        if gv == 'q':
            ax_omb.set_xlabel('RH O$-$B', size=14)
        else:
            ax_omb.set_xlabel('%s O$-$B var' % gv, size=14)
        ax_cts.set_xlabel('ob counts', size=14)
        for ax in [ax_omb, ax_cts]:
            ax.grid()
            ax.set_yscale('log')
            if typ in [153] + list(range(180, 200)) + list(range(280, 300)):
                ax.set_ylim([1100, 600])
            else:
                ax.set_ylim([1100, 10])
        if (gv == 'ps') or (gv == 'v'):
            ax_omb.legend()

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
