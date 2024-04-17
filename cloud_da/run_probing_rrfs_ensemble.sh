#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 04:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=10GB

# Note: Sometimes running this script returns the following error: 
#
# _tkinter.TclError: couldn't connect to display "localhost:11.0"
#
# I'm not entirely sure why this happens sometimes. It doesn't seem to pertain to a certain partition
# (I've had the code run successfully and return an error on sjet in separate instances). Submitting
# during a Jet session using port forwarding (e.g., with "jet_tunnel") with sjet appears to work.

date
. ~/.bashrc
pygraf

input_files=( '/lfs4/BMC/wrfruc/murdzek/src/py_scripts/cloud_da/cases/S_TX_2022020109/S_TX_2022020109_input.yml' 
              '/lfs4/BMC/wrfruc/murdzek/src/py_scripts/cloud_da/cases/Dakotas_2022020113/Dakotas_2022020113_input.yml'
	      '/lfs4/BMC/wrfruc/murdzek/src/py_scripts/cloud_da/cases/S_NewEngland_2022020121/S_NewEngland_2022020121_input.yml' 
	      '/lfs4/BMC/wrfruc/murdzek/src/py_scripts/cloud_da/cases/CONUS_2022020113/CONUS_2022020113_input.yml' )

input_files=( '/lfs4/BMC/wrfruc/murdzek/src/py_scripts/cloud_da/cases/S_TX_2022020109/S_TX_2022020109_input.yml' )

for f in ${input_files[@]}; do
  echo
  echo '===================================================='
  echo ${f}
  python probing_rrfs_ensemble.py ${f}
done

report-mem
date
