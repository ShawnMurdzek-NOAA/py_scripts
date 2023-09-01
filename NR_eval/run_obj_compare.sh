#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 08:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=4GB
#SBATCH --partition=orion

date
. ~/.bashrc
my_py

python obj_based_mrms_cref_compare.py

date
