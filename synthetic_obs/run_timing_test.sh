#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 02:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=25GB

date
. ~/.bashrc
#pygraf
adb_graphics

python timing_tests.py

report-mem
date
