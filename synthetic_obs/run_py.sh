#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 02:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=25GB

date
. ~/.bashrc
pygraf

python create_conv_obs_bary.py

report-mem
date
