"""
Single Observation Ceilometer Test Using an EnSRF

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import sys
import numpy as np

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
    
    return cld_amt, hofx_output


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
    param = pre.read_input_yaml(yml_fname)
    ens_obj = pre.read_ensemble_output(param)

    # Preprocess ceiling fields
    ens_obj = pre.preprocess_model_ceil(ens_obj, ceil_fields, ceil_names, ceil_miss)
    ens_obj = pre.preprocess_obs_ceil(ens_obj)

    # Apply cloud DA forward operator
    cld_amt, hofx = run_cld_forward_operator_1ob(ens_obj, param, 
                                                 ens_name=ens_obj.mem_names,
                                                 hofx_kw={'hgt_lim_kw':{'max_hgt':3500},
                                                          'verbose':0})
    print('Cloud ceilometer H(x) =', hofx)
    print('Cloud ceilometer ob =', cld_amt)

    # Run EnKF
    enkf_obj = enkf.enkf_1ob(ens_obj.state_matrix['data'], cld_amt, hofx, param['ob_var'])
    enkf_obj.EnSRF()

    # Compute analysis increments
    inc_1d = enkf_obj.x_a - enkf_obj.x_b
    inc_nd = unravel_state_matrix(inc_1d, ens_obj)


"""
End single_ceilometer_ob_enkf.py
"""