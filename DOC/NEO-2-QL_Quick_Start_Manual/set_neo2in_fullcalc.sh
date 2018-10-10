#! /bin/bash

fname_neo2in="neo2.in"

dirs="$(find . -type d -name "es_*" -printf "%f\n" | sort)"
#echo "${dirs}"

for dir in ${dirs};
do
    cd "${dir}"
    echo "Processing ${fname_neo2in} in ${dir}/ ..."
    if [ -e "${fname_neo2in}.sav" ];
    then
       echo "${fname_neo2in}.sav exists. Do nothing!"
    else
       if [ -e ${fname_neo2in} ];
       then
          cat ${fname_neo2in} > "${fname_neo2in}.sav"
	  
	  cat ${fname_neo2in} | awk '{\
	  sub("ISW_CALC_MAGDRIFT=          0","ISW_CALC_MAGDRIFT=          1"); \
	  sub("LSW_READ_PRECOM=F","LSW_READ_PRECOM=T"); \
	  sub("LSW_WRITE_PRECOM=T","LSW_WRITE_PRECOM=F"); \
	  sub("ISW_MAG_SHEAR=          0","ISW_MAG_SHEAR=          1"); \
	  print $0}' > ${fname_neo2in}
       fi
    fi
    sleep 0.1
    cd ..
done