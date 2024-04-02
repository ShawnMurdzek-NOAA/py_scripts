#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 04:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=10GB

date
. ~/.bashrc
pygraf

input_files=( '/lfs4/BMC/wrfruc/murdzek/src/py_scripts/cloud_da/cases/S_TX_2022020109/S_TX_2022020109_input.yml' 
              '/lfs4/BMC/wrfruc/murdzek/src/py_scripts/cloud_da/cases/Dakotas_2022020113/Dakotas_2022020113_input.yml'
	      '/lfs4/BMC/wrfruc/murdzek/src/py_scripts/cloud_da/cases/S_NewEngland_2022020121/S_NewEngland_2022020121_input.yml' 
	      '/lfs4/BMC/wrfruc/murdzek/src/py_scripts/cloud_da/cases/CONUS_2022020113/CONUS_2022020113_input.yml' )

for f in ${input_files[@]}; do
  echo
  echo '===================================================='
  echo ${f}
  python probing_rrfs_ensemble.py ${f}
done

report-mem
date
