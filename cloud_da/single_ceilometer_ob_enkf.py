"""
Single Observation Ceilometer Test Using an EnSRF

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import sys
import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import datetime as dt

import probing_rrfs_ensemble as pre
import pyDA_utils.cloud_DA_forward_operator as cfo
from pyDA_utils import enkf


#---------------------------------------------------------------------------------------------------
# Functions
#---------------------------------------------------------------------------------------------------

def run_cld_forward_operator_1ob(ens_obj, param, ens_name=['mem0001'], hofx_kw={}, verbose=False):
    
    hofx_output = np.zeros(len(ens_name))

    # Select observation for DA
    bufr_df = ens_obj._subset_bufr(['ADPSFC', 'MSONET'])
    dum = bufr_df.loc[bufr_df['SID'] == param['ob_sid'], :]
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
        hofx_output[i] = cld_hofx.data['hofx'][0][param['ob_idx']]
        cld_amt = cld_hofx.data['ob_cld_amt'][0][param['ob_idx']]
    
    return cld_amt, hofx_output, cld_ob_df


def unravel_state_matrix(x, ens_obj):
    """
    Unravel state matrix so fields can be plotted
    """

    output = {}
    for v in np.unique(ens_obj.state_matrix['vars']):
        output[v] = {}
        var_cond = ens_obj.state_matrix['vars'] == v
        for i, ens in enumerate(ens_obj.mem_names):
            output[v][ens] = np.reshape(x[var_cond, i], ens_obj.subset_ds[ens][v].shape)

    return output


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


def plot_horiz_postage_stamp(ens_obj, param, upp_field='TCDC_P0_L105_GLC0', klvl=0, 
                             ob={'plot':False}):
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
    if (upp_field == 'TCDC_P0_L105_GLC0') or (upp_field == 'ana_TCDC_P0_L105_GLC0'):
        if upp_field == 'TCDC_P0_L105_GLC0':
            save_fname = f"{param['out_dir']}/postage_stamp_cloud_cover_bgd_{param['save_tag']}.png"
        else:
            save_fname = f"{param['out_dir']}/postage_stamp_cloud_cover_ana_{param['save_tag']}.png"
        title = 'Cloud Cover'
        cmap = 'plasma_r'
        clevels = np.arange(0, 100.1, 5)
        extend = 'neither'
    elif upp_field == 'incr_TCDC_P0_L105_GLC0':
        save_fname = f"{param['out_dir']}/postage_stamp_cloud_cover_incr_{param['save_tag']}.png"
        title = 'Cloud Cover'
        cmap = 'bwr'
        clevels = np.arange(-97.5, 97.6, 5)

    elif (upp_field == 'TMP_P0_L105_GLC0') or (upp_field == 'ana_TMP_P0_L105_GLC0'):
        if upp_field == 'TMP_P0_L105_GLC0':
            save_fname = f"{param['out_dir']}/postage_stamp_T_bgd_{param['save_tag']}.png"
        else:
            save_fname = f"{param['out_dir']}/postage_stamp_T_ana_{param['save_tag']}.png"
        title = 'Temperature'
        cmap = 'plasma'
        clevels = np.arange(250, 270, 1)
    elif upp_field == 'incr_TMP_P0_L105_GLC0':
        save_fname = f"{param['out_dir']}/postage_stamp_T_incr_{param['save_tag']}.png"
        title = 'Temperature'
        cmap = 'bwr'
        clevels = np.arange(-5.75, 5.76, 0.5)
    
    elif (upp_field == 'SPFH_P0_L105_GLC0') or (upp_field == 'ana_SPFH_P0_L105_GLC0'):
        if upp_field == 'SPFH_P0_L105_GLC0':
            save_fname = f"{param['out_dir']}/postage_stamp_Q_bgd_{param['save_tag']}.png"
        else:
            save_fname = f"{param['out_dir']}/postage_stamp_Q_ana_{param['save_tag']}.png"
        title = 'Specific Humidity'
        cmap = 'plasma'
        clevels = np.arange(0, 5e-3, 2.4e-4)
        extend = 'max'
    elif upp_field == 'incr_SPFH_P0_L105_GLC0':
        save_fname = f"{param['out_dir']}/postage_stamp_Q_incr_{param['save_tag']}.png"
        title = 'Specific Humidity'
        cmap = 'bwr'
        clevels = np.arange(-5.75, 5.76, 0.5) * 1e-4

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
    param = pre.read_input_yaml(yml_fname)
    ens_obj = pre.read_ensemble_output(param)

    # Preprocess ceiling fields
    ens_obj = pre.preprocess_model_ceil(ens_obj, ceil_fields, ceil_names, ceil_miss)
    ens_obj = pre.preprocess_obs_ceil(ens_obj)

    # Apply cloud DA forward operator
    cld_amt, hofx, cld_ob_df = run_cld_forward_operator_1ob(ens_obj, param, 
                                                            ens_name=ens_obj.mem_names,
                                                            hofx_kw={'hgt_lim_kw':{'max_hgt':3500},
                                                                     'verbose':0},
                                                            verbose=False)
    print('Cloud ceilometer H(x) =', hofx)
    print('Cloud ceilometer ob =', cld_amt)

    # Run EnKF
    enkf_obj = enkf.enkf_1ob(ens_obj.state_matrix['data'], cld_amt, hofx, param['ob_var'])
    enkf_obj.EnSRF()

    # Save output to ens_obj for easier plotting
    ens_obj = add_inc_and_analysis_to_ens_obj(ens_obj, enkf_obj)

    # Make plots
    for meta in ['', 'incr_', 'ana_']:
        for v in param['state_vars']:
            print(f'plotting {meta}{v}...')
            fig = plot_horiz_postage_stamp(ens_obj, param, upp_field=f'{meta}{v}', klvl=param['plot_klvl'],
                                           ob={'plot':True,
                                               'x':cld_ob_df['XOB'].values[0] - 360, 
                                               'y':cld_ob_df['YOB'].values[0], 
                                               'kwargs':{'marker':'*', 'color':'k'}})
            plt.close(fig)
    
    print(f'total elapsed time = {(dt.datetime.now() - start).total_seconds()} s')


"""
End single_ceilometer_ob_enkf.py
"""