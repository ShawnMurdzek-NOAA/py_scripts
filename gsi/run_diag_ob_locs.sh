
# Load Python environment
module purge
module use -a /contrib/miniconda3/modulefiles
module load miniconda3/4.12.0
conda activate pygraf
export PYTHONPATH=$PYTHONPATH:/mnt/lfs4/BMC/wrfruc/murdzek/src
module list
conda list

# Dates and layers to loop over
dates=( 2024032411
        2024032412
	2024032423
	2024032500
	2024032518
	2024032606 )

layers=( 'all'
         'upper_trop'
	 'near_sfc'
	 'pbl' )

# Run Python script
echo
for d in ${dates[@]}; do
  echo ${d}
  for l in ${layers[@]}; do
    echo ${l}
    python diag_ob_locs.py ${d} ${l}
  done
done
