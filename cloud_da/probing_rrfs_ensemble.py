"""
Examining Cloud Information in the RRFS Ensemble

This script can either be run stand alone or within the probing_rrfs_ensemble.ipynb jupyter notebook.

Command-Line Inputs:
    yml_fname = Input YAML file name

shawn.s.murdzek@noaa.gov
Environment: pygraf (Jet)
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import yaml
import metpy.calc as mc
from metpy.units import units
import metpy.constants as const
import pandas as pd
import matplotlib
import copy
import sys

import pyDA_utils.plot_model_data as pmd
import pyDA_utils.bufr as bufr
import pyDA_utils.ensemble_utils as eu
import pyDA_utils.upp_postprocess as uppp
import pyDA_utils.cloud_DA_forward_operator as cfo
import pyDA_utils.cloud_DA_forward_operator_viz as cfov


#---------------------------------------------------------------------------------------------------
# Functions
#---------------------------------------------------------------------------------------------------

def read_input_yaml(fname):
    """
    Read YAML inputs 

    Parameters
    ----------
    fname : string
        YAML input file
    
    Returns
    -------
    param : dictionary
        YAML inputs

    """

    with open(fname, 'r') as fptr:
        param = yaml.safe_load(fptr)
    
    return param


def read_ensemble_output(param, verbose=1):
    """
    Reads UPP ensemble output

    Parameters
    ----------
    param : dictionary
        YAML inputs
    verbose : int, optional
        Verbosity level, by default 1

    Returns
    -------
    ens_obj : pyDA_utils.ensemble_utils.ensemble
        Ensemble object

    """

    # UPP output for each ensemble member
    str_format = param['str_format']
    prslev_fnames = {}
    natlev_fnames = {}
    for i in range(1, param['nmem']+1):
        prslev_fnames['mem{num:04d}'.format(num=i)] = str_format.format(num=i, lev='prslev')
        natlev_fnames['mem{num:04d}'.format(num=i)] = str_format.format(num=i, lev='natlev')
    prslev_vars = param['prslev_vars']

    # Read ensemble output
    start = dt.datetime.now()
    ens_obj = eu.ensemble(natlev_fnames, 
                          extra_fnames=prslev_fnames, 
                          extra_fields=prslev_vars, 
                          bufr_csv_fname=param['bufr_fname'], 
                          lat_limits=[param['min_lat'], param['max_lat']], 
                          lon_limits=[param['min_lon'], param['max_lon']],
                          zind=param['z_ind'], 
                          state_fields=param['state_vars'], 
                          bec=param['do_bec'])
    if verbose > 0:
        print('elapsed time = {t:.2f} s'.format(t=(dt.datetime.now() - start).total_seconds()))
        print('Shape of subset =', ens_obj.subset_ds['mem0001']['TMP_P0_L105_GLC0'].shape)

    return ens_obj


def preprocess_model_ceil(ens_obj, ceil_fields, ceil_names, ceil_miss):
    """
    Convert model ceiling heights to AGL and set missing values to NaN

    Parameters
    ----------
    ens_obj : pyDA_utils.ensemble_utils.ensemble
        Ensemble object
    ceil_fields : list of strings
        UPP ceiling fields
    ceil_names : list of strings
        Descriptive ceiling names
    ceil_miss : list of floats
        Missing values for each ceiling field

    Returns
    -------
    ens_obj : pyDA_utils.ensemble_utils.ensemble
        Ensemble object
    
    Notes
    -----
    This should really be added to pyDA_utils.ensemble_utils as a method
        
    """

    for m in ens_obj.mem_names:
        print(f'computing ceiling AGL heights for {m}')
        ens_obj.subset_ds[m] = uppp.compute_ceil_agl(ens_obj.subset_ds[m], no_ceil=np.nan)
        #for old_name, new_name, miss in zip(ceil_fields, ceil_names, ceil_miss):
        #    ens_obj.subset_ds[m][new_name] = ens_obj.subset_ds[m][old_name]
        #    ens_obj.subset_ds[m][new_name].values = (mc.geopotential_to_height(ens_obj.subset_ds[m][old_name].values * units.m * const.g).magnitude -
        #                                             mc.geopotential_to_height(ens_obj.subset_ds[m]['HGT_P0_L1_GLC0'].values * units.m * const.g).magnitude)
        #    if not np.isnan(miss):
        #        ens_obj.subset_ds[m][new_name].values[np.isclose(ens_obj.subset_ds[m][old_name].values, miss)] = np.nan
        #    ens_obj.subset_ds[m][new_name].attrs['long_name'] = 'ceiling height'
        #    ens_obj.subset_ds[m][new_name].attrs['units'] = 'm AGL'

    return ens_obj


def preprocess_obs_ceil(ens_obj, bufr_no_ceil=2e4):
    """
    Remove Obs With Missing Ceilings, Then Set No Ceiling to -1

    Parameters
    ----------
    ens_obj : pyDA_utils.ensemble_utils.ensemble
        Ensemble object
    bufr_no_ceil : float, optional
        Value that initially corresponds to "no ceiling" in BUFR file, by default 2e4

    Returns
    -------
    ens_obj : pyDA_utils.ensemble_utils.ensemble
        Ensemble object

    Notes
    -----
    This should really be added to pyDA_utils.ensemble_utils as a method
        
    """

    ens_obj.verif_obs.df = ens_obj.verif_obs.df.loc[~np.isnan(ens_obj.verif_obs.df['CEILING'])]
    ens_obj.verif_obs.df.reset_index(inplace=True, drop=True)
    ens_obj.verif_obs.df.loc[np.isclose(ens_obj.verif_obs.df['CEILING'], bufr_no_ceil), 'CEILING'] = -1

    return ens_obj


def plot_cloud_cover_horiz_postage_stamp(ens_obj, param):
    """
    Make horizontal cross section postage stamp plots: Cloud Cover at a Single Vertical Level

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
    upp_field = 'TCDC_P0_L105_GLC0'
    save_fname = '{d}/postage_stamp_cloud_cover_{tag}.png'.format(d=param['out_dir'], tag=param['save_tag'])
    title = 'Cloud Cover'
    klvl = 0
    cmap = 'plasma_r'
    clevels = np.arange(0, 100.1, 5)
    extend = 'neither'

    nrows = 5
    ncols = 6
    figsize = (12, 10)

    # Make plot
    fig = ens_obj.postage_stamp_contourf(upp_field, nrows, ncols, klvl=klvl, figsize=figsize, title=title,
                                         plt_kw={'ingest_kw':{'zind':[klvl]}, 
                                                 'cntf_kw':{'cmap':cmap, 'levels':clevels, 'extend':extend}})
    plt.savefig(save_fname)

    return fig


def plot_ceiling_postage_stamp(ens_obj, param, upp_field, title, plot_obs=True, bufr_field='CEILING'):
    """
    Make horizontal cross section postage stamp plots: Ceilings (with BUFR obs)

    Parameters
    ----------
    ens_obj : pyDA_utils.ensemble_utils.ensemble
        Ensemble object
    param : dictionary
        YAML inputs
    upp_field : string
        UPP name of ceiling field
    title : string
        Title for the plot
    plot_obs : boolean
        Option to plot BUFR obs, by default True
    bufr_field : string, optional
        Name of observed ceiling field, by default 'CEILING'

    Returns
    -------
    fig : plt.figure
        Plot

    """

    # Field name and colorbar limits/colormap. Set klvl to NaN if 2D field
    save_fname = '{d}/postage_stamp_ceil_{tag}.png'.format(d=param['out_dir'], tag=param['save_tag'])
    klvl = np.nan
    clevels = np.arange(0, 1000.1, 100)
    extend = 'max'

    bufr_subset = ['ADPSFC']
    bufr_nonan = np.nan     # Only use a row if this obs is not NaN (set to NaN to not use)

    nrows = 5
    ncols = 6
    figsize = (12, 10)

    # Create colormap with "no ceiling" obs (-1) in white
    cmap = copy.copy(matplotlib.cm.plasma)
    cmap.set_under('white')

    # Make plot
    fig = ens_obj.postage_stamp_contourf(upp_field, nrows, ncols, klvl=klvl, figsize=figsize, title=title,
                                         plt_kw={'ingest_kw':{'zind':[klvl]}, 
                                                 'cntf_kw':{'cmap':cmap, 'levels':clevels, 'extend':extend}})
    
    if plot_obs:
        ens_obj.plot_bufr_obs(fig.get_axes()[:-1], bufr_field, bufr_subset, nonan_field=bufr_nonan, 
                              scatter_kw={'cmap':cmap, 'vmin':clevels[0], 'vmax':clevels[-1], 
                                          'edgecolors':'k', 'linewidths':0.75, 's':50})
        
    plt.savefig(save_fname)

    return fig


def plot_ens_mean_std(ens_obj, param, upp_field='TCDC_P0_L105_GLC0', title='Cloud Cover'):
    """
    Plot Ensemble Mean and Standard Deviation

    Parameters
    ----------
    ens_obj : pyDA_utils.ensemble_utils.ensemble
        Ensemble object
    param : dictionary
        YAML inputs
    upp_field : string
        UPP name of ceiling field, by default 'TCDC_P0_L105_GLC0'
    title : string
        Title for the plot, by default 'Cloud Cover'

    Returns
    -------
    fig : plt.figure
        Plot

    """

    # Field name and colorbar limits/colormap. Set klvl to NaN if 2D field
    save_fname = '{d}/ens_mean_std_{tag}.png'.format(d=param['out_dir'], tag=param['save_tag'])
    klvl = 0
    cmap = 'plasma_r'
    vmin = 0
    vmax = 100

    nrows = 1
    ncols = 2
    figsize = (8, 5)

    # Make plot
    fig = plt.figure(figsize=figsize)
    for i, stat in enumerate(['mean', 'std']):
        ax = ens_obj.plot_state_vector(stat, upp_field, fig, nrows, ncols, i+1, zind=klvl, 
                                    pcm_kw={'cmap':cmap, 'vmin':vmin, 'vmax':vmax})


    plt.suptitle(title, size=20)
    plt.savefig(save_fname)

    return fig


def plot_bec_horiz_1var(ens_obj, param, var='TCDC_P0_L105_GLC0'):
    """
    Plot Spatial Distributions of BECs For a Particular Gridpoint and Variable

    Parameters
    ----------
    ens_obj : pyDA_utils.ensemble_utils.ensemble
        Ensemble object
    param : dictionary
        YAML inputs
    var : str, optional
        variable to plot BECs for, by default 'TCDC_P0_L105_GLC0'

    Returns
    -------
    fig : plt.figure
        Plot

    """

    # Coordinate of point to plot BECs in relation to (this is the "target" point)
    plot_klvl = param['bec_klvl']
    save_fname = '{d}/BEC_spatial_dist_{tag}.png'.format(d=param['out_dir'], tag=param['save_tag'])

    # Determine gridpoint closest to the target (lat, lon) coordinate
    tmp_ds = ens_obj.subset_ds[ens_obj.mem_names[0]]
    lat_all = tmp_ds['gridlat_0'].values
    lon_all = tmp_ds['gridlon_0'].values
    target_i, target_j = np.unravel_index(np.argmin((lat_all - param['bec_lat'])**2 + 
                                                    (lon_all - param['bec_lon'])**2), lon_all.shape)

    # Determine indices corresponding with the target variable and location
    target_var_idx = np.where(ens_obj.state_matrix['vars'] == var)[0]
    target_idx = target_var_idx[0] + np.ravel_multi_index([[param['bec_klvl']], [target_i], [target_j]], 
                                                          tmp_ds[var].shape)

    # Make plot
    fig = plt.figure(figsize=(12, 7.5))
    for i, v in enumerate(param['state_vars']):
        ax = ens_obj.plot_state_vector('be_cov', v, fig, 1, len(param['state_vars']), i+1, zind=plot_klvl,
                                       bec_idx=target_idx, ctr_cmap_0=True,
                                       pcm_kw={'cmap':'bwr'},
                                       cbar_kw={'orientation':'horizontal'})
        ax.plot(lon_all[target_i, target_j], lat_all[target_i, target_j], 'k*', ms=10, transform=ccrs.PlateCarree())

    plt.savefig(save_fname)

    return fig


def plot_ens_deviation_hist(ens_obj, param):
    """
    Plot Histograms of Deviations from Ensemble Mean

    ens_obj : pyDA_utils.ensemble_utils.ensemble
        Ensemble object
    param : dictionary
        YAML inputs

    Returns
    -------
    fig : plt.figure
        Plot

    """

    save_fname = '{d}/error_hist_{tag}.png'.format(d=param['out_dir'], tag=param['save_tag'])

    fig, axes = plt.subplots(nrows=1, ncols=len(param['state_vars']), figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.08, right=0.99, top=0.88, wspace=0.2)
    for i, v in enumerate(param['state_vars']):
        axes[i] = ens_obj.plot_ens_dev_hist(v, axes[i], hist_kw={'bins':50})

    axes[0].set_ylabel('counts', size=14)
    plt.suptitle('Deviations from the Ensemble Mean', size=22)
    plt.savefig(save_fname)

    return fig


def plot_ceil_ob_ranks(ens_obj, param, ceil_names, bufr_field='CEILING', verbosity=1):

    klvl = np.nan
    save_fname_spatial = ['{d}/{ceil}_ranks_%s_{tag}.png'.format(d=param['out_dir'], tag=param['save_tag'], ceil=c) for c in ceil_names]
    save_fname_hist = ['{d}/{ceil}_ranks_%s_hist_{tag}.png'.format(d=param['out_dir'], tag=param['save_tag'], ceil=c) for c in ceil_names]
    save_fname_txt = '{d}/ceiling_stats_{tag}.txt'.format(d=param['out_dir'], tag=param['save_tag'])

    bufr_subset = ['ADPSFC']
    bufr_nonan = np.nan    # Only use a row if this ob is not NaN (set to NaN to not use)

    out_fptr = open(save_fname_txt, 'w')

    # Subset obs and remove missing values
    subset_obs = ens_obj._subset_bufr(bufr_subset, nonan_field=bufr_nonan)

    for mfield, fname_spatial, fname_hist in zip(ceil_names, save_fname_spatial, save_fname_hist):

        if verbosity > 0:
            print()
            print('----------------------------')
            print('Field = {f}'.format(f=mfield))
            print()
        out_fptr.write('\n----------------------------\n')
        out_fptr.write('Field = {f}\n'.format(f=mfield))
        start = dt.datetime.now()
        
        # Interpolate ensemble members to ob locations
        ens_interp = ens_obj.interp_model_2d(mfield, subset_obs['YOB'].values, subset_obs['XOB'].values-360., 
                                            zind=klvl, method='nearest', verbose=True)

        # Compute height from geopotential height
        ob_vals = subset_obs[bufr_field].values
        ens_vals = mc.geopotential_to_height(ens_interp.loc[:, ens_obj.mem_names].values * units.m * const.g).to('m').magnitude

        # Set "no ceiling" to an arbitrarily large number
        ob_vals[np.isclose(ob_vals, -1)] = 1e8
        ens_vals[np.isnan(ens_vals)] = 1e9

        # Determine ranks
        ranks = np.zeros(len(subset_obs))
        for i in range(len(ranks)):
            ranks[i] = np.sum(ob_vals[i] > ens_vals[i, :])
        
        # Partitioning: Correct nulls (neither obs nor ensemble has a ceiling)
        correct_nulls = np.sum(np.logical_and(np.isclose(ob_vals, 1e8),
                                            np.all(np.isclose(ens_vals, 1e9), axis=1)))
        if verbosity > 0:
            print('correct null ceilings = {nulls} (out of {tot})'.format(nulls=correct_nulls, tot=len(ob_vals)))
        out_fptr.write('correct null ceilings = {nulls} (out of {tot})\n'.format(nulls=correct_nulls, tot=len(ob_vals)))
        
        # Partitioning: Missed ceilings (obs has ceiling but ensemble doesn't)
        misses = np.sum(np.logical_and(~np.isclose(ob_vals, 1e8),
                                    np.all(np.isclose(ens_vals, 1e9), axis=1)))
        if verbosity > 0:
            print('missed ceilings = {misses} (out of {tot})'.format(misses=misses, tot=len(ob_vals)))
        out_fptr.write('missed ceilings = {misses} (out of {tot})\n'.format(misses=misses, tot=len(ob_vals)))

        # Partitioning: False alarm ceilings (ensemble has ceiling, but obs do not)
        falarm = np.sum(np.logical_and(np.isclose(ob_vals, 1e8),
                                    ~np.all(np.isclose(ens_vals, 1e9), axis=1)))
        if verbosity > 0:
            print('false alarm = {falarm} (out of {tot})'.format(falarm=falarm, tot=len(ob_vals)))
        out_fptr.write('false alarm ceilings = {falarm} (out of {tot})\n'.format(falarm=falarm, tot=len(ob_vals)))
        
        # Partitioning: Hits (both obs and ensemble have a ceiling)
        hit_idx = np.where(np.logical_and(~np.isclose(ob_vals, 1e8),
                                        ~np.all(np.isclose(ens_vals, 1e9), axis=1)))[0]
        hits = len(hit_idx)
        if verbosity > 0:
            print('hit ceilings = {hits} (out of {tot})'.format(hits=hits, tot=len(ob_vals)))
        out_fptr.write('hit ceilings = {hits} (out of {tot})\n'.format(hits=hits, tot=len(ob_vals)))
        
        # Make plots
        for idx, tag in zip([np.where(np.ones(ranks.size))[0], hit_idx], ['all', 'hits']):
            
            rank_subset = ranks[idx]
            
            # Plot ranks: Spatial plot
            if verbosity > 0:
                print('making plots')

            fig = plt.figure(figsize=(8, 8))
            ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
            for i in range(len(rank_subset)):
                if rank_subset[i] == 0:
                    m = 'v'
                elif rank_subset[i] == param['nmem']:
                    m = '^'
                else:
                    m = 'o'
                cax = ax.scatter(subset_obs.loc[i, 'XOB'], subset_obs.loc[i, 'YOB'], c=rank_subset[i], 
                                s=param['msize'], cmap='bwr', vmin=0, vmax=param['nmem'], marker=m, linewidths=1, edgecolors='k',
                                transform=ccrs.PlateCarree(), alpha=1)

            ax.set_title('Ob field: {ob_f}, Model field: {m_f}\n'.format(ob_f=bufr_field, m_f=mfield) +
                        'down triangle: rank 0, up triangle: rank {rmax}'.format(rmax=param['nmem']), size=20)
            ax.set_extent([param['min_lon'], param['max_lon'], param['min_lat'], param['max_lat']])
            ax.coastlines('50m', edgecolor='gray', linewidth=0.25)
            borders = cfeature.NaturalEarthFeature(category='cultural',
                                                    scale='50m',
                                                facecolor='none',
                                                name='admin_1_states_provinces')
            ax.add_feature(borders, linewidth=0.25, edgecolor='gray')
            #lakes = cfeature.NaturalEarthFeature(category='physical',
            #                                     scale='50m',
            #                                     facecolor='none',
            #                                     name='lakes')
            #ax.add_feature(lakes, linewidth=0.25, edgecolor='gray')

            cbar = plt.colorbar(cax, ax=ax, orientation='horizontal', pad=0.02, aspect=25)
            cbar.set_label('Ob rank', size=16)
            plt.savefig(fname_spatial % tag)

            # Plot ranks: Histogram
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
            ax.hist(rank_subset, bins=np.arange(0, param['nmem']+2), edgecolor='k', linewidth=1)
            ax.set_title('Rank Histogram: {v} (n = {n})'.format(v=bufr_field, n=len(rank_subset)), size=18)
            ax.set_xlabel('rank', size=14)
            ax.set_ylabel('count', size=14)
            ax.grid(axis='y')
            ax.set_xlim([0, param['nmem']+1])
            plt.savefig(fname_hist % tag)

        if verbosity > 0:
            print('elapsed time = {t:.2f} s'.format(t=(dt.datetime.now() - start).total_seconds()))
        
    out_fptr.close()

    return None


def run_cld_forward_operator(ens_obj, param, ens_name=['mem0001']):

    for n in ens_name:
        print(f'Running forward operator on ensemble member {n}')
        model_ds = ens_obj.subset_ds[n]
        bufr_obj = ens_obj._subset_bufr(['ADPSFC', 'MSONET'])

        cld_ob_df = cfo.find_bufr_cloud_obs(bufr_obj)

        # Interpolating model columns to obs locations
        cld_hofx = cfo.sfc_cld_forward_operator(cld_ob_df, model_ds, debug=0)
        cld_hofx.interp_model_col_to_ob()

        # Imposing height limits and min cld fraction
        cld_hofx.impose_hgt_limits(hgt_field='model_col_height_agl',
                                fields=['model_col_height_agl', 'model_col_TCDC_P0_L105_GLC0'])
        cld_hofx.impose_min_cld_frac()

        # Set clear obs to have HOCB = 50 m so they don't get removed
        for i in cld_hofx.data['idx']:
            if np.isclose(cld_hofx.data['CLAM'][i][0], 0):
                cld_hofx.data['HOCB'][i][0] = 50
        cld_hofx.impose_hgt_limits(hgt_field='HOCB',
                                fields=['CLAM', 'HOCB'])
        cld_hofx.clean_obs()
        for i in cld_hofx.data['idx']:
            if np.isclose(cld_hofx.data['CLAM'][i][0], 0):
                cld_hofx.data['HOCB'][i][0] = np.nan

        # Adding clear obs and interpolating model clouds in column to ob heights
        cld_hofx.add_clear_obs()
        cld_hofx.decode_ob_clam()
        cld_hofx.interp_model_to_obs()

        # Compute RMSD
        rmsd = cld_hofx.compute_OmB()

        # Make plots
        cld_hofx_viz = cfov.sfc_cld_forward_operator_viz(cld_hofx)

        cld_hofx_viz.scatterplot()
        plt.savefig(f"{param['out_dir']}/cld_amt_scatterplot_{n}_{param['save_tag']}.png")

        cld_hofx_viz.hist(plot_param={'field':'OmB'})
        plt.savefig(f"{param['out_dir']}/OmB_hist_{n}_{param['save_tag']}.png")

        cld_hofx_viz.vert_columns(lon_lim=[-124, -70], lat_lim=[21, 49])
        plt.savefig(f"{param['out_dir']}/hofx_vcols_{n}_{param['save_tag']}.png")

        cld_hofx_viz.composite_cld_cover(lon_lim=[param['min_lon'], param['max_lon']], 
                                         lat_lim=[param['min_lat'], param['max_lon']])
        plt.savefig(f"{param['out_dir']}/composite_cld_cover_{n}_{param['save_tag']}.png")
    
    return rmsd


if __name__ == '__main__':

    # YAML inputs
    yml_fname = sys.argv[1]
    
    # Ceiling field to examine
    # NOTE: the v0.6.2 RRFS output appears to have an older formulation for the cloud base, so we'll omit for now
    #ceil_fields = ['HGT_P0_L215_GLC0', 'CEIL_P0_L215_GLC0', 'CEIL_P0_L2_GLC0', 'HGT_P0_L2_GLC0']
    #ceil_names = ['ceiling', 'ceil_exp1', 'ceil_exp2', 'cld_base']
    #ceil_miss = [np.nan, np.nan, 20000, -5000]
    ceil_fields = ['HGT_P0_L215_GLC0', 'CEIL_P0_L215_GLC0', 'CEIL_P0_L2_GLC0']
    ceil_names = ['CEIL_LEGACY', 'CEIL_EXP1', 'CEIL_EXP2']
    ceil_miss = [np.nan, np.nan, 20000]

    # Read input data
    param = read_input_yaml(yml_fname)
    ens_obj = read_ensemble_output(param)

    # Preprocess ceiling fields
    ens_obj = preprocess_model_ceil(ens_obj, ceil_fields, ceil_names, ceil_miss)
    ens_obj = preprocess_obs_ceil(ens_obj)

    # Create plots
    if param['do_bec']:
        _ = plot_cloud_cover_horiz_postage_stamp(ens_obj, param)
        _ = plot_ceiling_postage_stamp(ens_obj, param, ceil_names[0], ceil_names[0])
        _ = plot_ens_mean_std(ens_obj, param)
        _ = plot_bec_horiz_1var(ens_obj, param)
        _ = plot_ens_deviation_hist(ens_obj, param)

       # Apply cloud DA forward operator
        rmsd = run_cld_forward_operator(ens_obj, param)
        print(f'Forward operator RMSD = {rmsd}')

    _ = plot_ceil_ob_ranks(ens_obj, param, ceil_names)


"""
End probing_rrfs_ensemble.py
"""
