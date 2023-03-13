#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 03:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=4GB

date
. ~/.bashrc
adb_graphics

field=( 'cref' 'precip1hr' )
domain=( 'all' 'easternUS' )
zoom=( 'regular' 'zoomed' )
for f in ${field[@]}; do
  for d in ${domain[@]}; do
    for k in ${!zoom[@]}; do
      python frequency_histograms_HRRR.py $f $d $k "./HRRR_${f}_${d}_${zoom[$k]}_freq_hist.png" 
    done
  done
done

report-mem
date
