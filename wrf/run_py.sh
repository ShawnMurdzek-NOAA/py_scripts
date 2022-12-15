#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 00:10:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=6GB

date
. ~/.bashrc
py_env

python plot_upper_air_hcrsxn_gif.py

report-mem
date
