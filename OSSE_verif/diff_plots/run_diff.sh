#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 04:00:00
#SBATCH --ntasks=1
#SBATCH --partition=orion

date

. ~/.bashrc
my_py

python plot_hcrsxn_diffs.py hcrsxn_diffs_input.yml

date
