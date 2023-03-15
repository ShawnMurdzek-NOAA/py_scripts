#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 04:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=4GB

date
. ~/.bashrc
adb_graphics

model=( 'NR' 'HRRR' )
field=( 'cref' 'precip1hr' )
domain=( 'all' 'easternUS' )
zoom=( 'regular' 'zoomed' )
for m in ${model[@]}; do
  for f in ${field[@]}; do
    for d in ${domain[@]}; do
      for k in ${!zoom[@]}; do
        python frequency_histograms.py $m $f $d $k "./${m}_${f}_${d}_${zoom[$k]}_freq_hist.png"
      done 
    done
  done
done

report-mem
date
