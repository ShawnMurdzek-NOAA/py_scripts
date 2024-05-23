"""
Plot Difference Between Two Modeling Systems

Useful for comparing the NR to RRFS forecasts

Odd CartoPy Bug: Sometimes, CartoPy fails to plot all the data from RRFS when using a Lambert
Conformal projection (this has been observed to happen with 700-mb T). To fix this problem, simply
use a different map projection, like Plate Carree. I'm not sure why this bug exists, but using
changing the map projection is a good enough workaround.

Input Arguments
---------------
argv[1] : Input YAML file

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
import xarray as xr
import cartopy.crs as ccrs
import yaml
import sys

import pyDA_utils.plot_model_data as pmd
import pyDA_utils.upp_postprocess as uppp


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Open YAML file
yaml_name = sys.argv[1]
with open(yaml_name, 'r') as fptr:
    param = yaml.safe_load(fptr)

# Extract input parameters
model_info = param['model_info']
field_info = param['field_info']
lat_lim = param['lat_lim']
lon_lim = param['lon_lim']
save_tag = param['save_tag']

# Set indices to None if not defined
for comp in model_info.keys():
    for m in model_info[comp].keys():
        if model_info[comp][m]['indices'] == [None]:
            model_info[comp][m]['indices'] = None

# Set CartoPy projection
proj = ccrs.PlateCarree()


#---------------------------------------------------------------------------------------------------
# Plot Data
#---------------------------------------------------------------------------------------------------

field_list = list(field_info.keys())

# Create plots
for comp_name in model_info.keys():
    mnames = list(model_info[comp_name].keys())
    print()
    print('='*40)
    print(f'Plotting {comp_name}')
    for fname1, fname2 in zip(model_info[comp_name][mnames[0]]['files'], 
                              model_info[comp_name][mnames[1]]['files']):
        print()
        print(f"model 1 = {fname1.split('/')[-1]}")
        print(f"model 2 = {fname2.split('/')[-1]}")
        ds = []
        for f, m in zip([fname1, fname2], mnames):
            ds.append(xr.open_dataset(f, engine='pynio'))

            # Compute derived fields
            if ('WSPD' in field_list) or ('WDIR' in field_list):
                ds[-1] = uppp.compute_wspd_wdir(ds[-1])
            if 'CEIL' in field_list:
                ds[-1] = uppp.compute_ceil_agl(ds[-1], no_ceil=2e4)
                ds[-1]['CEIL'] = ds[-1][model_info[comp_name][m]['ceil_name']]

        for field in field_list:
            print(field)
            prslev = field_info[field]['prslev']
            cmap = field_info[field]['cmap']
            lvls = np.arange(field_info[field]['lvls']['start'], 
                             field_info[field]['lvls']['stop'], 
                             field_info[field]['lvls']['step'])
            diff_lvls = np.arange(field_info[field]['diff_lvls']['start'], 
                                  field_info[field]['diff_lvls']['stop'], 
                                  field_info[field]['diff_lvls']['step'])
            extend = field_info[field]['extend']
            fig = plt.figure(figsize=(8, 11))
            zind = []
            out = []
            fhr = []

            # Plot each field separately
            for ind in range(2):

                if prslev > 0:
                    zind.append(np.where(ds[ind]['lv_ISBL0'] == prslev)[0][0])
                else:
                    zind.append(np.nan)
                out.append(pmd.PlotOutput([ds[ind]], 'upp', fig, 3, 1, ind+1, proj=proj))
                out[-1].contourf(field, cntf_kw={'cmap':cmap, 'extend':extend, 'levels':lvls}, 
                                 ingest_kw={'zind':[zind[-1]]})
                out[-1].config_ax(grid=False)
                if lat_lim[0] != None:
                    out[-1].set_lim(lat_lim[0], lat_lim[1], lon_lim[0], lon_lim[1])
                
                # Extract forecast hour
                try:
                    fhr.append(ds[ind][field].attrs['forecast_time'][0])
                except KeyError:
                    print('Warning: Missing fhr')
                    fhr.append(0)

                # Set title
                if prslev > 0:
                    out[-1].ax_title(txt=f'{prslev/100}-hPa, f{fhr[-1]:03d}h {mnames[ind]}', size=14)
                else:
                    out[-1].ax_title(txt=f'f{fhr[-1]:03d}h {mnames[ind]}', size=14)
        
            # Plot differences
            diff = pmd.PlotOutput(ds, 'upp', fig, 3, 1, 3, proj=proj)
            diff.plot_diff(field, auto=False, cntf_kw={'cmap':'bwr', 'extend':'both', 
                                                       'levels':diff_lvls},
                                              ingest_kw={'zind':zind, 
                                                         'indices':[model_info[comp_name][mnames[0]]['indices'], 
                                                                    model_info[comp_name][mnames[1]]['indices']]})
            diff.config_ax(grid=False)
            if lat_lim[0] != None:
                diff.set_lim(lat_lim[0], lat_lim[1], lon_lim[0], lon_lim[1])
            diff.ax_title(txt=f'difference ({mnames[0]} $-$ {mnames[1]})\n', size=14)
    
            plt.subplots_adjust(left=0.02, bottom=0.04, right=0.98, top=0.97, hspace=0.18)
    
            # Create output file name
            try:
                init_time = dt.datetime.strptime(ds[0][field].attrs['initial_time'], '%m/%d/%Y (%H:%M)')
            except KeyError:
                print('Warning: Missing init_time')
                init_time = dt.datetime(2000, 1, 1, 0, 0)
            if prslev == 0:
                plt.savefig(f"{comp_name}_{save_tag}_{init_time.strftime('%Y%m%d%H')}_f{fhr[0]:03d}_f{fhr[1]:03d}_{field.split('_')[0]}.png")
            else:
                plt.savefig(f"{comp_name}_{save_tag}_{init_time.strftime('%Y%m%d%H')}_f{fhr[0]:03d}_f{fhr[1]:03d}_{prslev/1e2}mb_{field.split('_')[0]}.png")
            plt.close()


"""
End plot_hcrsxn_diffs.py
""" 
