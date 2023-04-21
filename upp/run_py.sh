#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 05:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=25GB
#SBATCH --partition=orion

date
. ~/.bashrc
my_py

python plot_upp_ceil.py

report-mem
date
