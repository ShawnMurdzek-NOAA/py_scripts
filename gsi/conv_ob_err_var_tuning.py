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
import gsi_fcts as gsi
import datetime as dt


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Paths to NCO_dirs folder for real-data and OSSE RRFS runs
path_real = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/real_data/winter/NCO_dirs'
path_osse = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/tune_conv_ob_err/winter1/NCO_dirs'

dates = [dt.datetime(2022, 2, 1, 9) + dt.timedelta(hours=i) for i in range(18)]

initial_err_fname = '/work2/noaa/wrfruc/murdzek/real_obs/errtable.rrfs'
output_err_fname = './errtable.tmp'


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
        real_omb_fnames.append('%s/ptmp/prod/rrfs.%s/%s/diag_conv_%s_ges.%s.nc4' % 
                               (path_real, d.strftime('%Y%m%d'), d.strftime('%H'), v, d.strftime('%Y%m%d%H')))
        osse_omb_fnames.append('%s/ptmp/prod/rrfs.%s/%s/diag_conv_%s_ges.%s.nc4' % 
                               (path_osse, d.strftime('%Y%m%d'), d.strftime('%H'), v, d.strftime('%Y%m%d%H')))
omf_df = {}
omf_df['real'] = gsi.read_diag(real_omb_fnames)
omf_df['osse'] = gsi.read_diag(osse_omb_fnames)

# Read initial errtable
init_errtable = gsi.read_errtable(initial_err_fname)
new_errtable = init_errtable.copy()

# Compute obs error variance adjustment
for typ in np.unique(omf_df['osse']['Observation_Type']):
    for v, err in zip(['t', 'q', 'uv', 'pw', 'ps'], ['Terr', 'RHerr', 'UVerr', 'PWerr', 'PSerr']):
        print()
        print('computing error variances for type = %d, var = %s' % (typ, v))
        prs = init_errtable[typ]['prs'].values
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

            ob_var = {}
            for run in ['real', 'osse']:
                if v == 'uv':
                    ob_var[run] = np.var(np.concatenate([subset[run]['u_Obs_Minus_Forecast_adjusted'].values,
                                                         subset[run]['v_Obs_Minus_Forecast_adjusted'].values]))
                elif v == 'q':
                    obs_q = subset[run]['Observation'].values
                    background_q = obs_q - subset[run]['Obs_Minus_Forecast_adjusted'].values
                    qs = subset[run]['Forecast_Saturation_Spec_Hum'].values
                    omf_q = 10 * ((obs_q / qs) - (background_q / qs))
                    ob_var[run] = np.var(omf_q)
                else:
                    ob_var[run] = np.var(subset[run]['Obs_Minus_Forecast_adjusted'])
            new_errtable[typ][err][k] = max(0, init_errtable[typ][err][k] + (ob_var['real'] - ob_var['osse']))
            print('prs = %.1f, Real Var = %.3f, OSSE Var = %.3f, New OSSE Var = %.3f' % 
                  (prs[k], ob_var['real'], ob_var['osse'], new_errtable[typ][err][k]))

print('saving new errtable')
gsi.write_errtable(output_err_fname, new_errtable)

print('Done (elapsed time = %s s)' % (dt.datetime.now() - start).total_seconds())

"""
End conv_ob_err_var_tuning.py 
"""
