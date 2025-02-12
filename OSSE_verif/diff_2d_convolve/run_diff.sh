#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 08:00:00
#SBATCH --ntasks=1
#SBATCH --partition=orion

date

. ~/.bashrc
my_py

python -u diff_2d_convolve.py input_diff_2d_convolve.yml

date
