"""
Plot Timeseries of O-B RMSDs For a Single Variable

shawn.s.murdzek@noaa.gov
Date Created: 18 April 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

import pyDA_utils.gsi_fcts as gsi


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# O-B file names. Key is simulation name, value is a list of file names
# O-Bs are found in the "ges" files and O-As are found in the "anl" files
omb_tmpl1 = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data_rrfs-workflow_orion/spring/NCO_dirs/ptmp/prod/rrfs.%s/%s/diag_conv_t_ges.%s.nc4'
omb_tmpl2 = '/work2/noaa/wrfruc/murdzek/RRFS_OSSE/real_red_data_rrfs-workflow_orion/spring/NCO_dirs/ptmp/prod/rrfs.%s/%s/diag_conv_t_ges.%s.nc4'
dates = [dt.datetime(2022, 4, 29, 21) + dt.timedelta(hours=i) for i in range(159)]
omb_fnames = {'OSSE': [omb_tmpl1 % (d.strftime('%Y%m%d'), d.strftime('%H'), d.strftime('%Y%m%d%H')) for d in dates],
              'Real': [omb_tmpl2 % (d.strftime('%Y%m%d'), d.strftime('%H'), d.strftime('%Y%m%d%H')) for d in dates]}

# Observation types
#ob_types = [181, 183, 187, 192, 193]
ob_types = [187]

# Assimilated obs only?
assim_only = True

# Use adjusted values?
use_adjusted = False

# Output directory and string to add to output file names
out_dir = './'
out_str = 'spring_only187'


#---------------------------------------------------------------------------------------------------
# Compute Statistics
#---------------------------------------------------------------------------------------------------

smpl_key = list(omb_fnames.keys())[0]

omb_var = omb_fnames[smpl_key][0].split('/')[-1].split('_')[2]
print('Var = %s' % omb_var)

# Select variable name based on whether adjusted or unadjusted O-Bs are desired
if use_adjusted:
    vname = 'Obs_Minus_Forecast_adjusted'
    adjs_str = 'adjusted'
else:
    vname = 'Obs_Minus_Forecast_unadjusted'
    adjs_str = 'unadjusted'

# Extract data
omb_rmsd = {} 
omb_dates = {}
for sim in omb_fnames:
    omb_rmsd[sim] = np.zeros(len(omb_fnames[sim]))
    omb_dates[sim] = []
    for i, f in enumerate(omb_fnames[sim]):
        df = gsi.read_diag([f])

        # Only keep data for desired observation types
        partial_df = []
        for typ in ob_types:
            partial_df.append(df.loc[df['Observation_Type'] == typ].copy())
        df = pd.concat(partial_df)

        # Only keep data from assimilated obs (if desired)
        assim_str = 'all'
        if assim_only:
            df = df.loc[df['Analysis_Use_Flag'] == 1, :]
            assim_str = 'assim'

        # Compute RMSDs
        omb_rmsd[sim][i] = np.sqrt(np.mean(df[vname].values**2))

        # Extract datetime
        omb_dates[sim].append(dt.datetime.strptime(str(df['date_time'].values[0]), '%Y%m%d%H'))

# Plot data
print('plotting data...')
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
for key, c in zip(omb_rmsd, ['b', 'r']):
    ax.plot(omb_dates[key], omb_rmsd[key], label=f"{key} (avg = {np.mean(omb_rmsd[key]):.2e})", c=c, ls='-') 

ax.grid()
ax.legend(fontsize=12)
ax.set_ylabel(f"{omb_var} O$-$B", size=12)
ax.tick_params(axis='x', rotation=45)

plt.subplots_adjust(left=0.08, bottom=0.15, right=0.99, top=0.88)

# Add additional metrics
ttl = 'TYP ='
for typ in ob_types:
    ttl = f"{ttl} {typ:d},"
ttl = f"{ttl}\nassim = {assim_only}, adjusted = {use_adjusted}"

plt.suptitle(ttl, size=16) 

plt.savefig(f"{out_dir}/omb_rmsd_timeseries_{out_str}_{assim_str}_{adjs_str}.png")
plt.close() 


"""
End diag_timeseries_rmsd_omb.py
"""
