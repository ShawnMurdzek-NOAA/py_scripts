"""
Plot CFAD for UPP Output

shawn.s.murdzek@noaa.gov
Date Created: 8 February 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mcm

import pyDA_utils.plot_model_data as pmd


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

fname = '/scratch1/BMC/wrfruc/murdzek/nature_run_tests/nature_run_spring_alphah1/output/202204291300/UPP/wrfprs_202204291400.grib2'
save_fname = './cfad_cref_alphah1.png'
name = 'NR: alphah1'

# CFAD parameters
var = 'REFD_P0_L100_GLC0'
zvar = 'lv_ISBL0'
bin_edges = np.arange(15, 75, 5)


#---------------------------------------------------------------------------------------------------
# Plot CFAD
#---------------------------------------------------------------------------------------------------

# Sample colors from plasma colorbar and define plotting levels
plevels = np.array([1, 100, 500, 1000, 5000, 10000, 25000, 50000, 75000, 100000, 150000])
nlvl = len(plevels)
colors = [mcm.plasma(i / (nlvl-1)) for i in range(nlvl)]

fig = plt.figure(figsize=(8, 8))
ds = xr.open_dataset(fname, engine='pynio')
out = pmd.PlotOutput([ds], 'upp', fig, 1, 1, 1)
out.cfad(var, zvar, bins=bin_edges, prs=True, cntf_kw={'colors':colors, 'levels':plevels, 'extend':'max'})
out.ax_title(txt=name, size=14)
plt.savefig(save_fname)
plt.close()


"""
End plot_cfad.py
"""
