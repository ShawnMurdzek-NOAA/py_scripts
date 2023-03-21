#!/bin/sh

date
. ~/.bashrc

# Load environment
case $MACHINE in
"JET")
  source ${PY_ENV}/py_jet.env;;
"HERA")
  source ${PY_ENV}/py_hera.env;;
"ORION")
  source ${PY_ENV}/py_orion.env;;
esac

# Determine range of wrfnat files
first_UPP=`date '+%Y%m%d%H' --date="${cycle::8} ${cycle:8:2}00 -3 hours"`
last_UPP=`date '+%Y%m%d%H' --date="${cycle::8} ${cycle:8:2}00 +1 hours"`

if [ ${first_UPP} -lt ${start_UPP} ]; then
  first_UPP=${start_UPP} 
fi

if [ ${last_UPP} -gt ${end_UPP} ]; then
  last_UPP=${end_UPP} 
fi

echo "${cycle}, ${first_UPP}, ${last_UPP}" 

python create_conv_obs.py ${cycle} ${first_UPP} ${last_UPP}
python create_adpupa_obs.py ${cycle} ${first_UPP} ${last_UPP}
python merge_conv_adpupa.py ${cycle}00

report-mem
date
