#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 01:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition=orion

date
. ~/.bashrc
my_py
conda activate /work2/noaa/wrfruc/murdzek/conda/pygrib

days=( 32 33 34 35 36 37 38 )
#days=( 119 120 121 122 123 124 125 )
for d in ${days[@]}; do
  python convert_nc_to_grib2.py ${d}
done

date
