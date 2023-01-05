#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 00:10:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=6GB

date
. ~/.bashrc
adb_graphics

python plot_hcrsxn_xarray.py

report-mem
date
