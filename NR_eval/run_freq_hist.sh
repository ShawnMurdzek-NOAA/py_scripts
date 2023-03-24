#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 06:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem=4GB
#SBATCH --partition=orion

date
. ~/.bashrc
my_py

#model=( 'NR' 'HRRR' )
#field=( 'cref' 'precip1hr' )
#domain=( 'all' 'easternUS' )
model=( 'HRRR' )
field=( 'precip1hr' )
domain=( 'easternUS' )
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
