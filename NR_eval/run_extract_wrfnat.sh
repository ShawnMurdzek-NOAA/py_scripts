#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 06:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=4GB
#SBATCH --partition=orion

date
. ~/.bashrc
my_py

dates='20220202 20220203 20220204 20220205 20220206 20220207 20220208'
for d in ${dates}; do
  python extract_wrfnat_fields.py ${d} 
done

report-mem
date
