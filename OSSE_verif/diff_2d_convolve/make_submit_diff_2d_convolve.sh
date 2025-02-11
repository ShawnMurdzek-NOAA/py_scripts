
# Wrapper to create and submit a bunch of jobs to run diff_2d_convolve

py_script='/work2/noaa/wrfruc/murdzek/src/py_scripts/OSSE_verif/diff_2d_convolve/diff_2d_convolve.py'
run_script='/work2/noaa/wrfruc/murdzek/src/py_scripts/OSSE_verif/diff_2d_convolve/run_diff.sh'
yml_template='./input_diff_2d_convolve_TEMPLATE.yml'

# Forecast hours
fcst_hr=( '003' '006' )

# The following should all have the same number of entries
fields=( 'REFC_P0_L200_GLC0' 'CEIL' 'CEIL' )
binary=( 'False' 'True' 'True' )
binary_val=( 0 152.4 914.4 )

#-------------------------------------------------------------------------------

for i in "${!fields[@]}"; do
  for fhr in ${fcst_hr[@]}; do

    # Create working directory
    dir="${fields[i]}_${binary_val[i]}_f${fhr}"
    echo "making and submitting job for ${dir}"
    mkdir ${dir}

    # Copy YAML, Python script, and run script
    yml_name='input_diff_2d_convolve.yml'
    run_script_name='run_diff.sh'
    cp ${py_script} ./${dir}
    cp ${run_script} ./${dir}/${run_script_name}
    cp ${yml_template} ./${dir}/${yml_name}

    # Create YAML (can add other YAML placeholders here if desired)
    cd ${dir}
    sed -i "s={FCST_HR}=${fhr}=" ${yml_name}
    sed -i "s={FIELD_NAME}=${fields[i]}=" ${yml_name}
    sed -i "s={BINARY}=${binary[i]}=" ${yml_name}
    sed -i "s={BINARY_VAL}=${binary_val[i]}=" ${yml_name}

    # Submit job
    sbatch ${run_script_name}
    cd ..

  done
done
