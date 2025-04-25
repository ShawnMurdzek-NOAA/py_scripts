#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 05:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=10GB
#SBATCH --partition=orion

date
. ~/.bashrc
my_py

python -u plot_hcrsxn_rmse_diff.py input_rmse_diff_f03.yml 

date
