
python diag_omf_vs_val.py \
	/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data_rrfs-workflow_orion/winter/NCO_dirs/ptmp/prod/rrfs.%Y%m%d/%H/diag_conv_q_ges.%Y%m%d%H.nc4 netcdf \
	2022020109 \
	2022020723 \
	OSSE_ges \
	--ob_type 188 \
        --plot_log 0 \
        --pseudo_rh 0

python diag_omf_vs_val.py \
	/work2/noaa/wrfruc/murdzek/RRFS_OSSE/real_red_data_rrfs-workflow_orion/winter/NCO_dirs/ptmp/prod/rrfs.%Y%m%d/%H/diag_conv_q_ges.%Y%m%d%H.nc4 netcdf \
	2022020109 \
	2022020723 \
	real_ges \
	--ob_type 188 \
	--plot_log 0 \
	--pseudo_rh 0
