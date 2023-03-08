#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 02:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=4GB

date
. ~/.bashrc
adb_graphics

date='20220430 20220501 20220502 20220503 20220504 20220505'
for d in ${date}; do
  python extract_wrfnat_fields.py ${d}
done

report-mem
date
