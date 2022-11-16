#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 00:15:00
#SBATCH --ntasks=1
#SBATCH --partition=bigmem
#SBATCH --mem=12GB

date
. ~/.bashrc
py_env
python plot_wrf_gif.py
report-mem
date
