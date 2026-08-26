"""
Compare Obs Errors Between GSI and JEDI

Using GSI's errtable and JEDI's YAML files. Each JEDI YAML file must only contain a single 
observation type / variable pair and use the naming conventions from RRFSv2X 
(<msg type>_<var><type>.yaml)

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import datetime as dt
import sys
import argparse
import copy
import numpy as np
import pandas as pd

import pyDA_utils.gsi_fcts as gsi


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Input files and directories
gsi_errtable_fname = '/gpfs/f6/wrfruc/scratch/Shawn.S.Murdzek/aircft_retros/errtable.rrfs'
jedi_yaml_dir = '/gpfs/f6/wrfruc/scratch/Shawn.S.Murdzek/aircft_retros/rrfs-workflow/parm/observers'

# Variables and types to check
ob_vars = ['t', 'q', 'uv']
ob_types = [130, 131, 134, 135, 230, 231, 234, 235]


#---------------------------------------------------------------------------------------------------
# Program
#---------------------------------------------------------------------------------------------------

# Read GSI errtable
gsi_errtable = gsi.read_errtable(gsi_errtable_fname)

# Dictionary to go between GSI errtable variables and JEDI YAML variables
var_dict = {'t':'Terr', 'q':'RHerr', 'uv':'UVerr'}

# Loop over each variable and type
for t in ob_types:
    for v in ob_vars:

        # Read JEDI YAML
        yaml_fname = f"{jedi_yaml_dir}/aircft_{v}{t}.yaml"
        try:
            with open(yaml_fname, 'r') as fptr:
                
                # Extract JEDI errors
                for l in fptr:
                    try:
                        if l.strip()[:7] == 'errors:':
                            jedi_err = np.array([float(s) for s in l.strip()[9:-1].split(',')])
                    except IndexError:
                        continue
        except FileNotFoundError:
            print(f"  Missing file. Skipping: {yaml_fname}")
            continue

        # Extract GSI errors
        gsi_err = gsi_errtable[t][var_dict[v]]
    
        # Multipl GSI errors by 0.1 for q (b/c errors are in %/10)
        if v == 'q':
            gsi_err = gsi_err * 0.1

        # Check errors
        if np.allclose(gsi_err, jedi_err):
            print(f"++Errors match for {v} {t}")
        else:
            print(f"--Errors do not match for {v} {t}")
        

"""
End compare_gsi_jedi_errs.py
"""
