#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 08:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=30GB
#SBATCH --partition=orion

date
. ~/.bashrc
my_py

which python
python create_conv_obs.py

report-mem
date
