#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 03:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=25GB

date
. ~/.bashrc
#pygraf
adb_graphics

python create_adpupa_obs.py

report-mem
date
