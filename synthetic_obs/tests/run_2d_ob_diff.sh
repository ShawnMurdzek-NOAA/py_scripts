#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 02:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=25GB
#SBATCH --partition=orion

date
. ~/.bashrc
my_py

bufr_tag='rap rap_e rap_p'
subsets='all ADPSFC SFCSHP MSONET'
domain='all easternUS westernUS'
for b in ${bufr_tag}; do
  for s in ${subsets}; do
    for d in ${domain}; do
      python plot_ob_diffs_2d.py ${b} ob_diffs_${s}_%s_%s_%s_%s_%s.png ${s} ${d}
    done
  done  
done

report-mem
date
