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

# Different rap prepbufr tags and the time (in min from DHR=0) to terminate radiosonde drift
# calculations
tags=( 'rap' 'rap_e' 'rap_p' 'sfc' )
end_interp=( 50 25 70 )

# Compute 

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
for i in ${!tags[@]}; do
  if [ -e ${BUFR_DIR}/${cycle}00.${tags[i]}.prepbufr.csv ]; then

    # Determine if there are any ADPUPA observations
    nADPUPA=`grep -o 'ADPUPA' ${BUFR_DIR}/${cycle}00.${tags[i]}.prepbufr.csv | wc -l`
    echo "${BUFR_DIR}/${cycle}00.${tags[i]}.prepbufr.csv"
    echo "nADPUPA = ${nADPUPA}"
    echo

    # Run simulated observation creation code
    cd ${CODE_DIR}

    python create_conv_obs.py ${WRF_DIR} \
                              ${BUFR_DIR} \
                              ${OUT_DIR}/conv/ \
                              ${cycle} \
                              ${first_UPP} \
                              ${last_UPP} \
                              ${tags[i]}

    conv_error=$?
    if [ ${conv_error} -gt 0 ]; then
      exit ${conv_error}
    fi

    if [ ${nADPUPA} -gt 0 ]; then
      python create_adpupa_obs.py ${WRF_DIR} \
                                  ${BUFR_DIR} \
                                  ${OUT_DIR}/adpupa/ \
                                  ${cycle} \
                                  ${first_UPP} \
                                  ${last_UPP} \
                                  ${tags[i]} \
                                  ${end_interp[i]}
      adpupa_error=$?

      if [ ${adpupa_error} -eq 0 ]; then
        python merge_conv_adpupa.py ${OUT_DIR}/conv/${cycle}00.${tags[i]}.fake.prepbufr.csv \
                                    ${OUT_DIR}/adpupa/${cycle}00.${tags[i]}.fake.adpupa.csv \
                                    ${OUT_DIR}/perfect/${cycle}00.${tags[i]}.fake.prepbufr.csv
      
        python merge_conv_adpupa.py ${OUT_DIR}/conv/${cycle}00.${tags[i]}.real_red.prepbufr.csv \
                                    ${OUT_DIR}/adpupa/${cycle}00.${tags[i]}.real_red.adpupa.csv \
                                    ${OUT_DIR}/perfect/${cycle}00.${tags[i]}.real_red.prepbufr.csv

      else
        echo
        echo "error (${adpupa_error}) in create_adpupa_obs.py, using ADPUPA obs from create_conv_obs.py instead"
        echo
        cp ${OUT_DIR}/conv/${cycle}00.${tags[i]}.fake.prepbufr.csv ${OUT_DIR}/perfect/${cycle}00.${tags[i]}.fake.prepbufr.csv
        cp ${OUT_DIR}/conv/${cycle}00.${tags[i]}.real_red.prepbufr.csv ${OUT_DIR}/perfect/${cycle}00.${tags[i]}.real_red.prepbufr.csv
      fi

    else
      cp ${OUT_DIR}/conv/${cycle}00.${tags[i]}.fake.prepbufr.csv ${OUT_DIR}/perfect/${cycle}00.${tags[i]}.fake.prepbufr.csv
      cp ${OUT_DIR}/conv/${cycle}00.${tags[i]}.real_red.prepbufr.csv ${OUT_DIR}/perfect/${cycle}00.${tags[i]}.real_red.prepbufr.csv
    fi

  fi
done

date
