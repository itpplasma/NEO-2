#! /bin/bash

dirs=$(find . -type d -name "es_*" -printf "%f\n" | sort)

for dir in ${dirs};
do
    cd ${dir} 
    num_dat_files=$(ls -l *.dat 2>/dev/null | wc -l)
    
    if [ ${num_dat_files} -gt "0" ]; 
    then
       echo "Removing ${num_dat_files} *.dat files ..."
       rm *.dat
    fi

    if [ -e neo2_multispecies_out.h5 ]; then
       echo "Removing neo2_multispecies_out.h5 ..."
       rm neo2_multispecies_out.h5
    fi

    if [ -e collop.h5 ]; then
       echo "Removing collop.h5 ..."
       rm collop.h5
    fi

    cd ..
    sleep 0.1
done