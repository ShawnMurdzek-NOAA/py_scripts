"""
Plot NR, RRFS, and Simulated Observed Soundings at a Single Point

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import datetime as dt

import pyDA_utils.bufr as bufr
import pyDA_utils.plot_model_data as pmd


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Input files (use wrfnat files)
NR_files = ['/work2/noaa/wrfruc/murdzek/nature_run_winter/UPP/20220203/wrfnat_202202031200.grib2']
RRFS_files = ['/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data/winter_updated/NCO_dirs/ptmp/prod/rrfs.20220203/11/rrfs.t11z.natlev.f001.conus_3km.grib2']
sim_obs_files = ['/work2/noaa/wrfruc/murdzek/nature_run_winter/obs/corr_errors_1st_iter_v2/err_csv/202202031200.rap.fake.prepbufr.csv']

# Vertical indices to use for NR/RRFS and max POB (mb) to use for simulated obs
NR_ind = []
RRFS_ind = []
sim_obs_p_range = [1500, 100]

# Stations to plot
station_ids = ["'72520'"]

# Output file name tag
save_tag = ''


#---------------------------------------------------------------------------------------------------
# Make Plots
#---------------------------------------------------------------------------------------------------

for fnameNR, fnameRRFS, fnameOB in zip(NR_files, RRFS_files, sim_obs_files):
    ds_NR = xr.open_dataset(fnameNR, engine='pynio')
    ds_RRFS = xr.open_dataset(fnameRRFS, engine='pynio')
    df_obs = bufr.bufrCSV(fnameOB).df
    for sid in station_ids:

        # Create a DataFrame using only the desired station
        subset_obs = df_obs.loc[(df_obs['SID'] == sid) & (df_obs['subset'] == 'ADPUPA')].copy()
        subset_obs.reset_index(inplace=True)
        lon = subset_obs['XOB'].values[0] - 360.
        lat = subset_obs['YOB'].values[0]
 
        # Plot NR sounding
        fig = plt.figure(figsize=(8, 8))
        NR_out = pmd.PlotOutput([ds_NR], 'upp', fig, 1, 1, 1)
        NR_out.skewt(lon, lat, barbs=False)

        # Plot RRFS sounding
        RRFS_out = pmd.PlotOutput([ds_RRFS], 'upp', fig, 1, 1, 1)
        RRFS_out.skewt(lon, lat, barbs=False, skew=NR_out.skew, hodo_ax=NR_out.h,
                       Tplot_kw={'linewidth':1.5, 'color':'r', 'linestyle':':'},
                       TDplot_kw={'linewidth':1.5, 'color':'b', 'linestyle':':'},
                       Hplot_kw={'linewidth':2, 'linestyle':':'})

        # Plot simulated obs sounding
        Tcond = np.logical_and(subset_obs['TQM'] < 3, subset_obs['TYP'] == 120)
        Wcond = np.logical_and(subset_obs['WQM'] < 3, subset_obs['TYP'] == 220)
        NR_out.skew.plot(subset_obs.loc[Tcond, 'POB'], subset_obs.loc[Tcond, 'TOB'], 'r', linewidth=1.5, linestyle='--')
        NR_out.skew.plot(subset_obs.loc[Tcond, 'POB'], subset_obs.loc[Tcond, 'TDO'], 'b', linewidth=1.5, linestyle='--')

        # Save plot
        fhr = ds_RRFS['TMP_P0_L105_GLC0'].attrs['forecast_time'][0]
        init_time = dt.datetime.strptime(ds_RRFS['TMP_P0_L105_GLC0'].attrs['initial_time'], '%m/%d/%Y (%H:%M)')
        plt.savefig('skewt_%s_%s_%s_f%03d.png' % (save_tag, sid, init_time.strftime('%Y%m%d%H'), fhr))
        plt.close()
        

"""
End plot_NR_RRFS_ob_soundings.py 
"""
