"""
Single Observation Ceilometer Test Using an EnSRF

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import datetime as dt
import metpy.calc as mc
from metpy.units import units
import copy

import probing_rrfs_ensemble as pre
import pyDA_utils.cloud_DA_forward_operator as cfo
from pyDA_utils import enkf
import pyDA_utils.plot_model_data as pmd


#---------------------------------------------------------------------------------------------------
# Functions
#---------------------------------------------------------------------------------------------------

def read_preprocess_ens(yml_fname):
    """
    Read in and preprocess ensemble output

    This step is independent of the observation being assimilated, so it should only need to be done
    once
    """

    # Ceiling field to examine
    ceil_fields = ['HGT_P0_L215_GLC0', 'CEIL_P0_L215_GLC0', 'CEIL_P0_L2_GLC0']
    ceil_names = ['CEIL_LEGACY', 'CEIL_EXP1', 'CEIL_EXP2']
    ceil_miss = [np.nan, np.nan, 20000]

    # Read input data
    param = pre.read_input_yaml(yml_fname)
    ens_obj = pre.read_ensemble_output(param)

    # Preprocess ceiling fields
    ens_obj = pre.preprocess_model_ceil(ens_obj, ceil_fields, ceil_names, ceil_miss)
    ens_obj = pre.preprocess_obs_ceil(ens_obj)

    return ens_obj, param


def run_cld_forward_operator_1ob(ens_obj, ob_sid, ob_idx, ens_name=['mem0001'], hofx_kw={}, verbose=False):
    
    hofx_output = np.zeros(len(ens_name))

    # Select observation for DA
    bufr_df = ens_obj._subset_bufr(['ADPSFC', 'MSONET'])
    dum = bufr_df.loc[bufr_df['SID'] == ob_sid, :]
    cld_ob_df = cfo.remove_missing_cld_ob(dum)
    if verbose: print("Observation:\n", cld_ob_df.loc[:, ['TYP', 'SID', 'XOB', 'YOB', 'CLAM', 'HOCB']])
    
    # Run forward operator
    for i, n in enumerate(ens_name):
        if verbose: print(f'Running forward operator on ensemble member {n}')
        model_ds = ens_obj.subset_ds[n]

        # Create ceilometer forward operator object
        cld_hofx = cfo.ceilometer_hofx_driver(cld_ob_df, model_ds, **hofx_kw)
        if verbose: 
            print("cld_hofx.data['CLAM'] = ", cld_hofx.data['CLAM'])
            print("cld_hofx.data['HOCB'] = ", cld_hofx.data['HOCB'])
            print("cld_hofx.data['ob_cld_amt'] = ", cld_hofx.data['ob_cld_amt'])
            print("cld_hofx.data['hofx'] = ", cld_hofx.data['hofx'])
        hofx_output[i] = cld_hofx.data['hofx'][0][ob_idx]
        cld_amt = cld_hofx.data['ob_cld_amt'][0][ob_idx]
        cld_z = cld_hofx.data['HOCB'][0][ob_idx]
    
    return cld_amt, cld_z, hofx_output, cld_ob_df


def unravel_state_matrix(x, ens_obj, ens_dim=True):
    """
    Unravel state matrix so fields can be plotted
    """

    output = {}
    for v in np.unique(ens_obj.state_matrix['vars']):
        var_cond = ens_obj.state_matrix['vars'] == v
        if ens_dim:
            output[v] = {}
            for i, ens in enumerate(ens_obj.mem_names):
                output[v][ens] = np.reshape(x[var_cond, i], ens_obj.subset_ds[ens][v].shape)
        else:
            output[v] = np.reshape(x[var_cond], ens_obj.subset_ds[ens_obj.mem_names[0]][v].shape)

    return output


def compute_RH(ens_obj, prefix=''):
    """
    Compute RH for the ensemble subset spatial domain
    """

    for mem in ens_obj.mem_names:
        ens_obj.subset_ds[mem][f"{prefix}RH_P0_L105_GLC0"] = ens_obj.subset_ds[mem]['SPFH_P0_L105_GLC0'].copy()
        p = ens_obj.subset_ds[mem]["PRES_P0_L105_GLC0"].values * units.Pa
        T = ens_obj.subset_ds[mem][f"{prefix}TMP_P0_L105_GLC0"].values * units.K
        Q = ens_obj.subset_ds[mem][f"{prefix}SPFH_P0_L105_GLC0"].values * units.kg / units.kg
        ens_obj.subset_ds[mem][f"{prefix}RH_P0_L105_GLC0"].values = mc.relative_humidity_from_specific_humidity(p, T, Q).magnitude * 100
        ens_obj.subset_ds[mem][f"{prefix}RH_P0_L105_GLC0"].attrs['long_name'] = 'relative humidity'
        ens_obj.subset_ds[mem][f"{prefix}RH_P0_L105_GLC0"].attrs['units'] = '%'

    return ens_obj


def compute_ens_stats_3D(ens_obj, field, stat_fct={'mean_':np.mean, 'std_': np.std}):
    """
    Compute various ensemble stats and analysis increments for the desired 3D field
    (i.e., the field is not part of the state vector)

    Save output to the first ensemble member
    """

    mem_names = ens_obj.mem_names
    shape = ens_obj.subset_ds[mem_names[0]][field].shape
    full_array = np.zeros([shape[0], shape[1], shape[2], len(mem_names)])

    # Extract the desired field from each ensemble member
    for i, m in enumerate(mem_names):
        full_array[:, :, :, i] = ens_obj.subset_ds[m][field].values
    
    # Compute stats
    for f in stat_fct.keys():
        ens_obj.subset_ds[mem_names[0]][f"{f}{field}"] = ens_obj.subset_ds[mem_names[0]][field].copy()
        ens_obj.subset_ds[mem_names[0]][f"{f}{field}"].values = stat_fct[f](full_array, axis=3)

    return ens_obj


def compute_ens_incr_3D(ens_obj, field):
    """
    Compute increments for the desired field and save the mean increment to the first ensemble member
    """

    mem_names = ens_obj.mem_names
    incr_sum = np.zeros(ens_obj.subset_ds[mem_names[0]][field].shape)

    # Compute increment for each ensemble member
    for m in mem_names:
        ens_obj.subset_ds[m][f"incr_{field}"] = ens_obj.subset_ds[m][field].copy()
        ens_obj.subset_ds[m][f"incr_{field}"].values = ens_obj.subset_ds[m][f"ana_{field}"].values - ens_obj.subset_ds[m][field].values
        incr_sum = incr_sum + ens_obj.subset_ds[m][f"incr_{field}"].values
    
    # Add average increment to first ensemble member
    ens_obj.subset_ds[mem_names[0]][f"mean_incr_{field}"] = ens_obj.subset_ds[mem_names[0]][field].copy()
    ens_obj.subset_ds[mem_names[0]][f"mean_incr_{field}"].values = incr_sum / len(mem_names)

    return ens_obj


def add_inc_and_analysis_to_ens_obj(ens_obj, enkf_obj):
    """
    Add the analysis increment and analysis fields to the ens_obj for easier plotting
    """

    # Compute increment
    inc_1d = enkf_obj.x_a - enkf_obj.x_b

    # Turn 1D fields into nD fields
    inc_nd = unravel_state_matrix(inc_1d, ens_obj)
    xa_nd = unravel_state_matrix(enkf_obj.x_a, ens_obj)

    # Add fields to ens_obj
    for out, label in zip([inc_nd, xa_nd], ['incr_', 'ana_']):
        for v in out.keys():
            for ens in out[v].keys():
                ens_obj.subset_ds[ens][label+v] = ens_obj.subset_ds[ens][v].copy()
                ens_obj.subset_ds[ens][label+v].values = out[v][ens]

    return ens_obj


def add_ens_mean_std_K_to_ens_obj(ens_obj, enkf_obj):
    """
    Add the ensemble mean, standard deviation, and Kalman gain to the first ensemble member
    """

    # Compute stats
    var_2d = {}
    for x, label1 in zip([enkf_obj.x_b, enkf_obj.x_a], ['', 'ana']):
        for fct, label2 in zip([np.mean, np.std], ['mean', 'std']):
            if label1 == '':
                name = label2
            else:
                name = f'{label2}_{label1}'
            var_2d[name] = unravel_state_matrix(fct(x, axis=1), ens_obj, ens_dim=False)
    var_2d['K'] = unravel_state_matrix(enkf_obj.K, ens_obj, ens_dim=False)

    # Add fields to ens_obj
    ens = ens_obj.mem_names[0]
    for v in var_2d[list(var_2d.keys())[0]]:
        for key in var_2d.keys():
            ens_obj.subset_ds[ens][f"{key}_{v}"] = ens_obj.subset_ds[ens][v].copy()
            ens_obj.subset_ds[ens][f"{key}_{v}"].values = var_2d[key][v]
            if key == 'K':
                units = ens_obj.subset_ds[ens][f"{key}_{v}"].attrs['units']
                ens_obj.subset_ds[ens][f"{key}_{v}"].attrs['units'] = f"{units} / [obs units]"
                ens_obj.subset_ds[ens][f"{key}_{v}"].attrs['long_name'] = "Kalman gain"
        ens_obj.subset_ds[ens][f"mean_incr_{v}"] = ens_obj.subset_ds[ens][v].copy()
        ens_obj.subset_ds[ens][f"mean_incr_{v}"].values = var_2d['mean_ana'][v] - var_2d['mean'][v]

    return ens_obj


def plot_horiz_slices(ds, field, ens_obj, param, nrows=4, ncols=4, 
                      klvls=list(range(16)), figsize=(12, 12), verbose=False, cntf_kw={}, 
                      ob={'plot':False},
                      save_dir=None):
    """
    Plot several horizontal slices at various levels
    """
    
    # Make plot
    fig = plt.figure(figsize=figsize)
    for i, k in enumerate(klvls):
        if verbose:
            print(f"plotting k = {k}")
        plot_obj = pmd.PlotOutput([ds], 'upp', fig, nrows, ncols, i+1)

        # Make filled contour plot
        # Skip plotting if < 2 NaN
        make_plot = np.sum(~np.isnan(ds[field][k, :, :])) > 1
        if make_plot:
            plot_obj.contourf(field, cbar=False, ingest_kw={'zind':[k]}, cntf_kw=cntf_kw)
            cax = plot_obj.cax
            meta = plot_obj.metadata['contourf0']
        else:
            if verbose:
                print(f"skipping plot for k = {k}")
            plot_obj.ax = fig.add_subplot(nrows, ncols, i+1, projection=plot_obj.proj)
        
        # Add location of observation
        if ob['plot']:
            plot_obj.plot(ob['x'], ob['y'], plt_kw=ob['kwargs'])

        plot_obj.config_ax(grid=False)
        plot_obj.set_lim(ens_obj.lat_limits[0], ens_obj.lat_limits[1], 
                         ens_obj.lon_limits[0], ens_obj.lon_limits[1])
        title = 'avg z = {z:.1f} m'.format(z=float(np.mean(ds['HGT_P0_L105_GLC0'][k, :, :] -
                                                            ds['HGT_P0_L1_GLC0'])))
        plot_obj.ax.set_title(title, size=14)
    
    plt.subplots_adjust(left=0.01, right=0.85)

    cb_ax = fig.add_axes([0.865, 0.02, 0.02, 0.9])
    cbar = plt.colorbar(cax, cax=cb_ax, orientation='vertical', aspect=35)
    cbar.set_label(f"{meta['name']} ({meta['units']})", size=14)

    plt.suptitle(field, size=18)

    # Save figure
    if save_dir is None:
        plt.savefig(f"{param['out_dir']}/{field}_{param['save_tag']}.png")
    else:
        plt.savefig(f"{save_dir}/{field}_{param['save_tag']}.png")

    return fig


def plot_horiz_postage_stamp(ens_obj, param, upp_field='TCDC_P0_L105_GLC0', klvl=0, 
                             ob={'plot':False}, save_dir=None):
    """
    Make horizontal cross section postage stamp plots

    Parameters
    ----------
    ens_obj : pyDA_utils.ensemble_utils.ensemble
        Ensemble object
    param : dictionary
        YAML inputs

    Returns
    -------
    fig : plt.figure
        Plot

    """

    # Field name and colorbar limits/colormap. Set klvl to NaN if 2D field
    extend = 'both'
    if save_dir is None:
        save_dir = param['out_dir']
    if (upp_field == 'TCDC_P0_L105_GLC0') or (upp_field == 'ana_TCDC_P0_L105_GLC0'):
        if upp_field == 'TCDC_P0_L105_GLC0':
            save_fname = f"{save_dir}/postage_stamp_cloud_cover_bgd_{param['save_tag']}.png"
        else:
            save_fname = f"{save_dir}/postage_stamp_cloud_cover_ana_{param['save_tag']}.png"
        title = 'Cloud Cover'
        cmap = 'plasma_r'
        clevels = np.arange(0, 100.1, 5)
        extend = 'neither'
    elif upp_field == 'incr_TCDC_P0_L105_GLC0':
        save_fname = f"{save_dir}/postage_stamp_cloud_cover_incr_{param['save_tag']}.png"
        title = 'Cloud Cover'
        cmap = 'bwr'
        clevels = np.arange(-97.5, 97.6, 5)

    elif (upp_field == 'TMP_P0_L105_GLC0') or (upp_field == 'ana_TMP_P0_L105_GLC0'):
        if upp_field == 'TMP_P0_L105_GLC0':
            save_fname = f"{save_dir}/postage_stamp_T_bgd_{param['save_tag']}.png"
        else:
            save_fname = f"{save_dir}/postage_stamp_T_ana_{param['save_tag']}.png"
        title = 'Temperature'
        cmap = 'plasma'
        clevels = np.arange(250, 270, 1)
    elif upp_field == 'incr_TMP_P0_L105_GLC0':
        save_fname = f"{save_dir}/postage_stamp_T_incr_{param['save_tag']}.png"
        title = 'Temperature'
        cmap = 'bwr'
        clevels = np.arange(-5.75, 5.76, 0.5)
    
    elif (upp_field == 'SPFH_P0_L105_GLC0') or (upp_field == 'ana_SPFH_P0_L105_GLC0'):
        if upp_field == 'SPFH_P0_L105_GLC0':
            save_fname = f"{save_dir}/postage_stamp_Q_bgd_{param['save_tag']}.png"
        else:
            save_fname = f"{save_dir}/postage_stamp_Q_ana_{param['save_tag']}.png"
        title = 'Specific Humidity'
        cmap = 'plasma'
        clevels = np.arange(0, 5e-3, 2.4e-4)
        extend = 'max'
    elif upp_field == 'incr_SPFH_P0_L105_GLC0':
        save_fname = f"{save_dir}/postage_stamp_Q_incr_{param['save_tag']}.png"
        title = 'Specific Humidity'
        cmap = 'bwr'
        clevels = np.arange(-5.75, 5.76, 0.5) * 1e-4

    elif (upp_field == 'RH_P0_L105_GLC0') or (upp_field == 'ana_RH_P0_L105_GLC0'):
        if upp_field == 'RH_P0_L105_GLC0':
            save_fname = f"{save_dir}/postage_stamp_RH_bgd_{param['save_tag']}.png"
        else:
            save_fname = f"{save_dir}/postage_stamp_RH_ana_{param['save_tag']}.png"
        title = 'Relative Humidity'
        cmap = 'plasma'
        clevels = np.arange(0, 100.1, 5)
        extend = 'neither'
    elif upp_field == 'incr_RH_P0_L105_GLC0':
        save_fname = f"{save_dir}/postage_stamp_RH_incr_{param['save_tag']}.png"
        title = 'Relative Humidity'
        cmap = 'bwr'
        clevels = np.arange(-5.75, 5.76, 0.5) * 4

    nrows = 5
    ncols = 6
    figsize = (12, 10)

    # Make plot
    fig = ens_obj.postage_stamp_contourf(upp_field, nrows, ncols, klvl=klvl, figsize=figsize, title=title,
                                         plt_kw={'ingest_kw':{'zind':[klvl]}, 
                                                 'cntf_kw':{'cmap':cmap, 'levels':clevels, 'extend':extend}})
    
    # Add location of observation
    if ob['plot']:
        for ax in fig.axes:
            if type(ax) == cartopy.mpl.geoaxes.GeoAxes:
                ax.plot(ob['x'], ob['y'], transform=ccrs.PlateCarree(), **ob['kwargs'])

    plt.savefig(save_fname)

    return fig


if __name__ == '__main__':

    start = dt.datetime.now()

    # Read and preprocess ensemble
    ens_obj, param = read_preprocess_ens(sys.argv[1])
    ens_obj_original = copy.deepcopy(ens_obj)

    # Loop over each observation
    for ob_sid, ob_idx in zip(param['ob_sid'], param['ob_idx']):

        start_loop = dt.datetime.now()
        print()
        print('-----------------------------------------------')
        print(f"Starting single-ob test for {ob_sid} {ob_idx}")

        # Apply cloud DA forward operator
        cld_amt, cld_z, hofx, cld_ob_df = run_cld_forward_operator_1ob(ens_obj, ob_sid, ob_idx, 
                                                                    ens_name=ens_obj.mem_names,
                                                                    hofx_kw={'hgt_lim_kw':{'max_hgt':3500},
                                                                                'verbose':0},
                                                                    verbose=False)
        print('Cloud ceilometer ob hgt =', cld_z)
        print('Cloud ceilometer ob amt =', cld_amt)
        print('Cloud ceilometer H(x) =', hofx)
        print(f"Time to complete forward operator = {(dt.datetime.now() - start_loop).total_seconds()} s")

        # Run EnKF
        enkf_obj = enkf.enkf_1ob(ens_obj.state_matrix['data'], cld_amt, hofx, param['ob_var'])
        enkf_obj.EnSRF()
        print(f"Time to complete forward operator and EnSRF = {(dt.datetime.now() - start_loop).total_seconds()} s")

        # Save output to ens_obj for easier plotting
        ens_obj = add_inc_and_analysis_to_ens_obj(ens_obj, enkf_obj)
        ens_obj = add_ens_mean_std_K_to_ens_obj(ens_obj, enkf_obj)

        # Compute RH as well as ensemble stats for RH
        for p in ['', 'ana_']:
            ens_obj = compute_RH(ens_obj, prefix=p)
            ens_obj = compute_ens_stats_3D(ens_obj, f"{p}RH_P0_L105_GLC0")
        ens_obj = compute_ens_incr_3D(ens_obj, "RH_P0_L105_GLC0")
        param['state_vars'].append('RH_P0_L105_GLC0')

        # Plot ensemble mean and standard deviation
        save_dir = f"{param['out_dir']}/{ob_sid}_{ob_idx}"
        os.system(f"mkdir {save_dir}")
        if param['plot_ens_stats']:
            print('Making ensemble mean and standard deviation plots')
            for meta, cmap in zip(['mean_', 'mean_ana_', 'mean_incr_', 'std_', 'std_ana_', 'K_'], 
                                ['plasma', 'plasma', 'bwr', 'plasma', 'plasma', 'bwr']):
                for v in param['state_vars']:
                    field = meta + v
                    if field not in ens_obj.subset_ds[ens_obj.mem_names[0]]:
                        print(f'field {field} is missing. Skipping.')
                        continue
                    print(f'plotting {field}...')
                    maxval = np.percentile(np.abs(ens_obj.subset_ds[ens_obj.mem_names[0]][field]), 99)
                    minval = np.percentile(np.abs(ens_obj.subset_ds[ens_obj.mem_names[0]][field]), 1)
                    if cmap == 'bwr':
                        clevels = np.linspace(-1, 1, 30) * maxval
                    else:
                        clevels = np.linspace(minval, maxval, 30)
                    if (maxval == 0) and (minval == 0):
                        clevels = np.linspace(-1, 1, 30)
                    fig = plot_horiz_slices(ens_obj.subset_ds[ens_obj.mem_names[0]], 
                                            field,
                                            ens_obj,
                                            param,
                                            nrows=param['plot_nrows'],
                                            ncols=param['plot_ncols'],
                                            klvls=param['plot_klvls'], 
                                            cntf_kw={'cmap':cmap, 'levels':clevels, 'extend':'both'},
                                            ob={'plot':True,
                                                'x':cld_ob_df['XOB'].values[0] - 360, 
                                                'y':cld_ob_df['YOB'].values[0], 
                                                'kwargs':{'marker':'*', 'color':'k'}},
                                            save_dir=save_dir)
                    plt.close(fig)

        # Make postage stamp plots
        if param['plot_postage_stamp']:
            print('Making postage stamp plots')
            for meta in ['', 'incr_', 'ana_']:
                for v in param['state_vars']:
                    print(f'plotting {meta}{v}...')
                    fig = plot_horiz_postage_stamp(ens_obj, param, upp_field=f'{meta}{v}', 
                                                klvl=param['plot_postage_stamp_klvl'],
                                                ob={'plot':True,
                                                    'x':cld_ob_df['XOB'].values[0] - 360, 
                                                    'y':cld_ob_df['YOB'].values[0], 
                                                    'kwargs':{'marker':'*', 'color':'k'}},
                                                save_dir=save_dir)
                    plt.close(fig)
        
        # Clean up
        print(f'total time for {ob_sid} {ob_idx} (forward operator, EnSRF, and plots) = {(dt.datetime.now() - start_loop).total_seconds()} s')
        ens_obj = copy.deepcopy(ens_obj_original)

    print(f'total elapsed time = {(dt.datetime.now() - start).total_seconds()} s')


"""
End single_ceilometer_ob_enkf.py
"""
