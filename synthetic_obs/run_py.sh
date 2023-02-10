#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 00:15:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition=tjet
#SBATCH --mem=4GB

date
. ~/.bashrc
adb_graphics

python create_conv_obs_wgts_load.py

report-mem
date
