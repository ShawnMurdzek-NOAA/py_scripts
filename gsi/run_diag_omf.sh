#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 08:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=20GB
#SBATCH --partition=orion

date
. ~/.bashrc
my_py

python diag_omf_cts_boxplots.py

date
