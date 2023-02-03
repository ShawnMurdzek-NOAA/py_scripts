"""
Perform a Quick Sanity Check for Various UPP Fields

shawn.s.murdzek@noaa.gov
Date Created: 24 January 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import cartopy.crs as ccrs


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

upp_nat = '/scratch1/BMC/wrfruc/murdzek/nature_run_spring_mynn2/output/202204291300/UPP/wrfnat_202204291400.grib2'
upp_prs = '/scratch1/BMC/wrfruc/murdzek/nature_run_spring_mynn2/output/202204291300/UPP/wrfprs_202204291400.grib2'
save_fname = 'upp_fields.pdf'

# Hybrid level fields (excluding the "_P0_L105_GLC0" suffix)
hyb_fields = ['TMP', 'SPFH', 'CLWMR', 'ICMR', 'RWMR', 'SNMR', 'GRLE', 'SPNCR', 'SPNCS',
              'SPNCG', 'UGRD', 'VGRD', 'VVEL', 'PRES', 'HGT', 'NCONCD', 'NCCICE', 'FRACCC',
              'REFD', 'TKE', 'DZDT']
hyb_z = [20] * len(hyb_fields)

# Isobaric fields (excluding the "_P0_L100_GLC0" suffix)
iso_fields = ['TMP', 'DPT', 'UGRD', 'VGRD', 'HGT', 'DZDT']
iso_z = [5] * len(iso_fields)

# Other 3D fields
o3d_fields = ['UGRD_P0_L103_GLC0', 'VGRD_P0_L103_GLC0', 'CAPE_P0_2L108_GLC0', 'CIN_P0_2L108_GLC0']
o3d_z = [0] * len(o3d_fields)

# 2D fields
flat_fields = ['PRES_P0_L1_GLC0', 'TMP_P0_L103_GLC0', 'DPT_P0_L103_GLC0', 'HGT_P0_L215_GLC0',
               'VIS_P0_L1_GLC0', 'PWAT_P0_L200_GLC0', 'REFC_P0_L200_GLC0', 'HPBL_P0_L1_GLC0',
               'RETOP_P0_L3_GLC0', 'HGT_P0_L3_GLC0', 'CAPE_P0_L1_GLC0', 'CIN_P0_L1_GLC0',
               'ASNOW_P8_L1_GLC0_acc', 'WEASD_P8_L1_GLC0_acc', 'APCP_P8_L1_GLC0_acc',
               'TCDC_P0_L211_GLC0', 'HCDC_P0_L234_GLC0', 'MCDC_P0_L224_GLC0', 'LCDC_P0_L214_GLC0',
               'CEIL_P0_L2_GLC0', 'CEIL_P0_L215_GLC0']


#---------------------------------------------------------------------------------------------------
# Create Plots
#---------------------------------------------------------------------------------------------------

# Add suffixes to UPP field names
for i in range(len(hyb_fields)):
    hyb_fields[i] = hyb_fields[i] + '_P0_L105_GLC0'
for i in range(len(iso_fields)):
    iso_fields[i] = iso_fields[i] + '_P0_L100_GLC0'

# Create master lists
hyb_fields = hyb_fields + o3d_fields + flat_fields
hyb_z = hyb_z + o3d_z + [np.nan] * len(flat_fields)

# Create figures
pdf = PdfPages(save_fname)
for upp, fields, zind in zip([upp_nat, upp_prs], [hyb_fields, iso_fields], [hyb_z, iso_z]):
    ds = xr.open_dataset(upp, engine='pynio')
    lat = ds['gridlat_0'].values
    lon = ds['gridlon_0'].values
    for f, zi in zip(fields, zind):
        print(f)
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
        if np.isnan(zi):
            cax = ax.contourf(lon, lat, ds[f])
            ax.set_title(f)
        else:
            cax = ax.contourf(lon, lat, ds[f][zi, :, :])
            ax.set_title('%s, %d' % (f, zi))
        ax.coastlines('50m')
        cbar = plt.colorbar(cax, orientation='vertical')
        cbar.set_label('%s (%s)' % (ds[f].attrs['long_name'], ds[f].attrs['units']), size=14)
        pdf.savefig(fig)
        plt.close()

pdf.close()


"""
End check_UPP_fields.py 
"""
