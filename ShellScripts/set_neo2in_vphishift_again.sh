#! /bin/bash

fname_vphiref="vphiref.in"

if [ -e ${fname_vphiref} ];
then
   vphishift=$(cat "${fname_vphiref}" | awk '/^vphiref/ {gsub(/[^0-9.\-]/,""); print}')
   echo -e "vphiref = ${vphishift}"
else
   echo -e "Enter rotation frame frequency shift (e.g., 5000 or 1000): "
   read vphishift
   echo -e "vphiref = ${vphishift}" > ${fname_vphiref}
   sync
   sleep 0.1
   sync
   echo -e "\n"
fi

fname_neo2in="neo2.in"

#dirs="$(find . -type d -name "es_*" -printf "%f\n" | sort)"
dirs="$(ls -d es_*)"

for dir in ${dirs};
do
  cd "${dir}"
  echo "Processing ${fname_neo2in} in ${dir}/ ..."
  if [ -e "${fname_neo2in}.sav" ]; then
    if [ ! -e neo2_multispecies_out.h5 ]; then
      echo -e "Output file neo2_multispecies_out.h5 not found in directory ${dir}"

      vphi=$(cat "${fname_neo2in}.sav" | awk '/^ VPHI=/ {gsub(/[^0-9.\-]/,""); print}')

      vphi_new=$(echo -e "${vphi} + ${vphishift}" | bc)

      cat "${fname_neo2in}.sav" | awk '{\
      gsub(/^ VPHI=[ 0-9.\-]+/," VPHI='${vphi_new}'"); \
      sub("LSW_READ_PRECOM=F","LSW_READ_PRECOM=T"); \
      sub("LSW_WRITE_PRECOM=T","LSW_WRITE_PRECOM=F"); \
      print}' > ${fname_neo2in}
    fi
  else
    echo "${fname_neo2in}.sav does not exist. Save existing neo2.in as neo2.in.err!"
    echo "This directory will be excluded from Condor runs! Check input files and run NEO-2 manually!"
    cat ${fname_neo2in} > "${fname_neo2in}.err"
  fi
  sync
  sleep 0.1
  sync
  cd ..
done
