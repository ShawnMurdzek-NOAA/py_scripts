
module use /usw/conda/modulefiles
module load miniforge
conda activate /ncrc/home2/Shawn.S.Murdzek/.conda/envs/my_py
export PYTHONPATH=$PYTHONPATH:/ncrc/home2/Shawn.S.Murdzek/src/

ens=('mem001' 'mem002' 'mem003')
for m in ${ens[@]}; do
  echo
  echo ${m}
  python plot_mpas_hcrsxn.py ../../direct_ceilometer_DA/tests/sample_data/mpas/${m}/mpasout.2024-05-27_04.00.00.TEST.nc \
	                     ../../direct_ceilometer_DA/tests/sample_data/mpas/invariant_TEST.nc \
			     --outtag bgd_${m} \
                             -l 2

  python plot_mpas_hcrsxn.py ../../direct_ceilometer_DA/tests/sample_data/mpas/${m}/mpasout.DA.2024-05-27_04.00.00.TEST.nc \
	                     ../../direct_ceilometer_DA/tests/sample_data/mpas/invariant_TEST.nc \
			     --outtag diff_${m} \
			     -l 2 \
                             --file2 ../../direct_ceilometer_DA/tests/sample_data/mpas/${m}/mpasout.2024-05-27_04.00.00.TEST.nc \
			     -c bwr \
                             --vmin -0.25 \
			     --vmax 0.25

  python plot_mpas_hcrsxn.py ../../direct_ceilometer_DA/tests/sample_data/mpas/${m}/mpasout.DA.2024-05-27_04.00.00.TEST.nc \
	                     ../../direct_ceilometer_DA/tests/sample_data/mpas/invariant_TEST.nc \
			     --outtag diff_${m} \
			     -l 2 \
                             --file2 ../../direct_ceilometer_DA/tests/sample_data/mpas/${m}/mpasout.2024-05-27_04.00.00.TEST.nc \
			     -c bwr \
                             --vmin -1 \
			     --vmax 1 \
			     -f theta

  python plot_mpas_hcrsxn.py ../../direct_ceilometer_DA/tests/sample_data/mpas/${m}/mpasout.DA.2024-05-27_04.00.00.TEST.nc \
	                     ../../direct_ceilometer_DA/tests/sample_data/mpas/invariant_TEST.nc \
			     --outtag diff_${m} \
			     -l 2 \
                             --file2 ../../direct_ceilometer_DA/tests/sample_data/mpas/${m}/mpasout.2024-05-27_04.00.00.TEST.nc \
			     -c bwr \
                             --vmin -0.001 \
			     --vmax 0.001 \
			     -f qv
done
