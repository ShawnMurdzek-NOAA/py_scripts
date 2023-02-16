#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 02:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=25GB

date
. ~/.bashrc
adb_graphics

python create_conv_obs_wgts_load.py

report-mem
date
