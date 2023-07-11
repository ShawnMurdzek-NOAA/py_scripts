"""
Determine Error Variance Needed to Conventional Observation Error Tuning

shawn.s.murdzek@noaa.gov
Date Created: 18 May 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import numpy as np
import pandas as pd
import datetime as dt

import pyDA_utils.gsi_fcts as gsi


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Paths to NCO_dirs folder for real-data and OSSE RRFS runs
path_real = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/real_data/winter/NCO_dirs'
path_osse = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/tune_conv_ob_err/winter_DEBUG/NCO_dirs'

dates = [dt.datetime(2022, 2, 1, 3) + dt.timedelta(hours=i) for i in range(6)]

initial_err_var_fname = '/work2/noaa/wrfruc/murdzek/real_obs/errtable.rrfs'
output_err_var_fname = './errtable.tmp'

initial_err_mean_fname = '/work2/noaa/wrfruc/murdzek/real_obs/errtable_mean.rrfs'
output_err_mean_fname = './errtable_mean.tmp'

# Observation types (set to None for all observation types)
ob_types = [126]


#---------------------------------------------------------------------------------------------------
# Compare O-B Variances For Real and Simulated Obs
#---------------------------------------------------------------------------------------------------

# Extract data
start = dt.datetime.now()
print('start = %s' % start.strftime('%Y%m%d %H:%M:%S'))
print('reading GSI diag files...')
real_omb_fnames = []
osse_omb_fnames = []
for v in ['t', 'q', 'uv', 'pw', 'ps']:
    for d in dates:
        real_omb_fnames.append('%s/ptmp/prod/rrfs.%s/%s_spinup/diag_conv_%s_ges.%s.nc4' % 
                               (path_real, d.strftime('%Y%m%d'), d.strftime('%H'), v, d.strftime('%Y%m%d%H')))
        osse_omb_fnames.append('%s/ptmp/prod/rrfs.%s/%s_spinup/diag_conv_%s_ges.%s.nc4' % 
                               (path_osse, d.strftime('%Y%m%d'), d.strftime('%H'), v, d.strftime('%Y%m%d%H')))
omf_df = {}
omf_df['real'] = gsi.read_diag(real_omb_fnames)
omf_df['osse'] = gsi.read_diag(osse_omb_fnames)

# Read initial errtable
init_var_errtable = gsi.read_errtable(initial_err_var_fname)
new_var_errtable = init_var_errtable.copy()
init_mean_errtable = gsi.read_errtable(initial_err_mean_fname)
new_mean_errtable = init_mean_errtable.copy()

# Compute obs error variance adjustment
if ob_types == None:
    ob_types = np.unique(omf_df['osse']['Observation_Type'])

for typ in ob_types:
    for v, err in zip(['t', 'q', 'uv', 'pw', 'ps'], ['Terr', 'RHerr', 'UVerr', 'PWerr', 'PSerr']):
        print()
        print('computing error variances for type = %d, var = %s' % (typ, v))
        prs = init_var_errtable[typ]['prs'].values
        prs[1:] = 0.5 * (prs[:-1] + prs[1:])
        prs[0] = prs[1] + 0.5 * (prs[1] - prs[2])
        for k in range(len(prs) - 1):

            subset = {}
            for run in ['real', 'osse']:
                subset[run] = omf_df[run].loc[(omf_df[run]['Observation_Type'] == typ) & 
                                              (omf_df[run]['var'] == v) &
                                              (omf_df[run]['Pressure'] < prs[k]) & 
					      (omf_df[run]['Pressure'] >= prs[k+1])]
            if len(subset['osse']) == 0:
                continue 

            ob_err_stat = {}
            for run in ['real', 'osse']:
                ob_err_stat[run] = {}
                for stat, stat_fct in zip(['mean', 'var', 'n'], [np.mean, np.var, len]):
                    if v == 'uv':
                        ob_err_stat[run][stat] = stat_fct(np.concatenate([subset[run]['u_Obs_Minus_Forecast_adjusted'].values,
                                                                          subset[run]['v_Obs_Minus_Forecast_adjusted'].values]))
                    elif v == 'q':
                        obs_q = subset[run]['Observation'].values
                        background_q = obs_q - subset[run]['Obs_Minus_Forecast_adjusted'].values
                        qs = subset[run]['Forecast_Saturation_Spec_Hum'].values
                        omf_q = 10 * ((obs_q / qs) - (background_q / qs))
                        ob_err_stat[run][stat] = stat_fct(omf_q)
                    else:
                        ob_err_stat[run][stat] = stat_fct(subset[run]['Obs_Minus_Forecast_adjusted'])
            new_var_errtable[typ][err][k] = max(0, init_var_errtable[typ][err][k] + 
                                               (ob_err_stat['real']['var'] - ob_err_stat['osse']['var']))
            new_mean_errtable[typ][err][k] = (init_mean_errtable[typ][err][k] + 
                                              (ob_err_stat['real']['mean'] - ob_err_stat['osse']['mean']))
            print(('prs = %6.1f | Real n, mean, var = %6d, %10.3e, %10.3e | ' +
                   'OSSE n, mean, var = %6d, %10.3e, %10.3e | ' +
                   'Diff mean = %10.3e | New OSSE Var = %10.3e') % 
                  (prs[k], ob_err_stat['real']['n'], ob_err_stat['real']['mean'], ob_err_stat['real']['var'], 
                   ob_err_stat['osse']['n'], ob_err_stat['osse']['mean'], ob_err_stat['osse']['var'], 
                   new_mean_errtable[typ][err][k], new_var_errtable[typ][err][k]))

print('saving new errtables')
gsi.write_errtable(output_err_var_fname, new_var_errtable)
gsi.write_errtable(output_err_mean_fname, new_mean_errtable)

print('Done (elapsed time = %s s)' % (dt.datetime.now() - start).total_seconds())

"""
End conv_ob_err_tuning.py 
"""
