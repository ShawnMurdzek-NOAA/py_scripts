
# Create and submit jobs to create difference plots between two RRFS runs (or RRFS and NR)
#
# This script uses the plot_hcrsxn_diffs.py script found in py_scripts/OSSE_verif/diff_plots

py_script=/work2/noaa/wrfruc/murdzek/src/py_scripts/OSSE_verif/diff_plots/plot_hcrsxn_diffs.py
yaml_tmpl=hcrsxn_diffs_input.yml
slurm_script=/work2/noaa/wrfruc/murdzek/src/py_scripts/OSSE_verif/diff_plots/run_diff.sh

init_start=('2022020109'
	    '2022020112'
	    '2022020200'
	    '2022020212'
            '2022020300'
	    '2022020312'
            '2022020400'
	    '2022020412'
            '2022020500'
	    '2022020512'
            '2022020600'
	    '2022020612'
            '2022020700'
	    '2022020712')
init_end=('2022020111'
	  '2022020123'
	  '2022020211'
	  '2022020223'
	  '2022020311'
          '2022020323'
	  '2022020411'
          '2022020423'
	  '2022020511'
          '2022020523'
	  '2022020611'
          '2022020623'
	  '2022020711'
          '2022020723')

fname1='/work2/noaa/wrfruc/murdzek/RRFS_OSSE/syn_data_rrfs-workflow_orion/winter/NCO_dirs/ptmp/prod/rrfs.{iYYYY}{iMM}{iDD}/{iHH}/rrfs.t{iHH}z.prslev.f{FFF}.conus_3km.grib2'
fname2='/work2/noaa/wrfruc/murdzek/nature_run_winter/UPP/{vYYYY}{vMM}{vDD}/wrfprs_{vYYYY}{vMM}{vDD}{vHH}00_er.grib2'

#-------------------------------------------------------------------------------

home=`pwd`

for i in "${!init_start[@]}"; do
  subdir=${init_start[i]}
  echo "making inputs for ${subdir}"
  mkdir ${subdir}
  cp ${py_script} ${subdir}/
  cp ${slurm_script} ${subdir}/
  cp ${yaml_tmpl} ${subdir}/
  cd ${subdir}
  sed -i "s={FNAME1}=${fname1}=" ${yaml_tmpl}
  sed -i "s={FNAME2}=${fname2}=" ${yaml_tmpl}
  sed -i "s={ISTART}=${init_start[i]}=" ${yaml_tmpl}
  sed -i "s={IEND}=${init_end[i]}=" ${yaml_tmpl}
  sbatch ${slurm_script}
  cd ${home}
done
