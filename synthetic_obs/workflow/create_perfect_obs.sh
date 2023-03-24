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
first_UPP=`date '+%Y%m%d%H' --date="${cycle::8} ${cycle:8:2}00 -4 hours"`
last_UPP=`date '+%Y%m%d%H' --date="${cycle::8} ${cycle:8:2}00 +1 hours"`

if [ ${first_UPP} -lt ${start_UPP} ]; then
  first_UPP=${start_UPP} 
fi

if [ ${last_UPP} -gt ${end_UPP} ]; then
  last_UPP=${end_UPP} 
fi

echo "${cycle}, ${first_UPP}, ${last_UPP}" 

cd ${CODE_DIR}
for tag in 'rap rap_e rap_p'; do
  if [ -e ${BUFR_DIR}/${cycle}00.${tag}.prepbufr.csv ]; then

    # Determine if there are any ADPUPA observations
    nADPUPA=`grep -o ${BUFR_DIR}/${cycle}00.rap.prepbufr.csv | wc -l`
    echo "${tag}, nADPUPA = ${nADPUPA}"
    echo

    # Run simulated observation creation code
    cd ${CODE_DIR}
    python create_conv_obs.py ${WRF_DIR} ${BUFR_DIR} ${OUT_DIR}/conv/ ${cycle} ${first_UPP} ${last_UPP} ${tag}
    if [ ${nADPUPA} -gt 0 ]; then
      python create_adpupa_obs.py ${WRF_DIR} ${BUFR_DIR} ${OUT_DIR}/adpupa/ ${cycle} ${first_UPP} ${last_UPP} ${tag}
      python merge_conv_adpupa.py ${OUT_DIR}/conv/${cycle}00.${tag}.fake.prepbufr.csv ${OUT_DIR}/adpupa/${cycle}00.${tag}.fake.adpupa.csv ${OUT_DIR}/perfect/${cycle}00.${tag}.fake.prepbufr.csv
    else
      cp ${OUT_DIR}/conv/${cycle}00.${tag}.fake.prepbufr.csv ${OUT_DIR}/perfect/${cycle}00.${tag}.fake.prepbufr.csv
    fi

  fi
done

report-mem
date
