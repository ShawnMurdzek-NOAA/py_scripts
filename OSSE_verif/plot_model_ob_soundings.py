"""
Plot Soundings from Two Modeling Systems and Simulated Obs (Optional) at a Single Point

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
model1_files = ['/work2/noaa/wrfruc/murdzek/nature_run_spring/UPP/20220504/wrfnat_202205041600_er.grib2']
model2_files = ['/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data_rrfs-workflow_orion/spring_uas_35km/NCO_dirs/ptmp/prod/rrfs.20220504/16/rrfs.t16z.natlev.f000.conus_3km.grib2']
sim_obs_files = ['/work2/noaa/wrfruc/murdzek/nature_run_spring/obs/uas_obs_35km/combine_obs/combine_csv/202205041600.rap.fake.prepbufr.csv']

# Vertical indices to use for model1/model2 and max POB (mb) to use for simulated obs
model1_ind = []
model2_ind = []
sim_obs_p_range = [1500, 100]

# Name of model1 and model2
model1_name = 'NR'
model2_name = 'RRFS_spring_uas35'

# Location: RAOB station IDs or (lat, lon) coordinates
use_ob_site = False
ob_subset = 'AIRCAR'
station_ids = ["'UA003694'"]
lon_all = [-95.92877]
lat_all = [41.264866]

# Output file name tag
save_tag = 'uas35'


#---------------------------------------------------------------------------------------------------
# Make Plots
#---------------------------------------------------------------------------------------------------

for fnameM1, fnameM2, fnameOB in zip(model1_files, model2_files, sim_obs_files):
    ds_M1 = xr.open_dataset(fnameM1, engine='pynio')
    ds_M2 = xr.open_dataset(fnameM2, engine='pynio')
    
    if use_ob_site:
        nplot = len(station_ids)
        df_obs = bufr.bufrCSV(fnameOB).df
    else:
        nplot = len(lon_all)
    for j in range(nplot):
  
        print(f'plotting {j+1} of {nplot}')

        # Create a DataFrame using only the desired station
        if use_ob_site:
            subset_obs = df_obs.loc[(df_obs['SID'] == station_ids[j]) & (df_obs['subset'] == ob_subset)].copy()
            subset_obs.reset_index(inplace=True)
            lon = subset_obs['XOB'].values[0] - 360.
            lat = subset_obs['YOB'].values[0]
        else:
            lon = lon_all[j]
            lat = lat_all[j]
 
        # Plot model1 sounding
        fig = plt.figure(figsize=(8, 8))
        M1_out = pmd.PlotOutput([ds_M1], 'upp', fig, 1, 1, 1)
        M1_out.skewt(lon, lat, barbs=False,
                     Tplot_kw={'linewidth':1, 'color':'g', 'linestyle':'-'},
                     TDplot_kw={'linewidth':1, 'color':'g', 'linestyle':'-'},
                     Hplot_kw={'linewidth':1, 'linestyle':'-'})

        # Plot model2 sounding
        M2_out = pmd.PlotOutput([ds_M2], 'upp', fig, 1, 1, 1)
        M2_out.skewt(lon, lat, barbs=False, skew=M1_out.skew, hodo_ax=M1_out.h,
                       Tplot_kw={'linewidth':1.5, 'color':'r', 'linestyle':'--'},
                       TDplot_kw={'linewidth':1.5, 'color':'r', 'linestyle':'--'},
                       Hplot_kw={'linewidth':2, 'linestyle':'--'})

        # Plot simulated obs sounding
        if use_ob_site:
            Tcond = np.logical_and(subset_obs['TQM'] < 3, subset_obs['TYP'] == 120)
            Wcond = np.logical_and(subset_obs['WQM'] < 3, subset_obs['TYP'] == 220)
            M1_out.skew.plot(subset_obs.loc[Tcond, 'POB'], subset_obs.loc[Tcond, 'TOB'], 'b', linewidth=2, linestyle=':')
            M1_out.skew.plot(subset_obs.loc[Tcond, 'POB'], subset_obs.loc[Tcond, 'TDO'], 'b', linewidth=2, linestyle=':')

        # Add title
        if use_ob_site:
            plt.suptitle(f'Green Solid: {model1_name}\nRed Dashed: {model2_name}\nBlue Dotted: Simulated Observation', size=16)
            loc = station_ids[j]
        else:
            plt.suptitle(f'Green Solid: {model1_name}\nRed Dashed: {model2_name}', size=16)
            loc = f'{lat}N_{lon+360}E'

        # Save plot
        fhr = ds_M2['TMP_P0_L105_GLC0'].attrs['forecast_time'][0]
        init_time = dt.datetime.strptime(ds_M2['TMP_P0_L105_GLC0'].attrs['initial_time'], '%m/%d/%Y (%H:%M)')
        plt.savefig('skewt_%s_%s_%s_f%03d.png' % (save_tag, loc, init_time.strftime('%Y%m%d%H'), fhr))
        plt.close()
        

"""
End plot_model_ob_soundings.py 
"""
