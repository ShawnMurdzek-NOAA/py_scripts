#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 03:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition=xjet
#SBATCH --mem=30GB

date
. ~/.bashrc
adb_graphics

python create_conv_obs_wgts_load.py

report-mem
date
