#!/bin/bash

basename_directories="n2_vshift"

#dirs="$(find . -type d -name "run_aug_shift_*" -printf "%f\n" | sort)"
dirs=$(ls -d ${basename_directories}*)
#echo "${dirs}"

fname_NTV_torque_int_all="NTV_torque_int_all.dat"
fname_NTV_torque_int="NTV_torque_int.dat"

if [ -e ${fname_NTV_torque_int_all} ]; then
   rm ${fname_NTV_torque_int_all}
fi

for dir in ${dirs}; do
    if [ -e ${fname_NTV_torque_int_all} ]; then
       cat "${dir}/${fname_NTV_torque_int}" >> ${fname_NTV_torque_int_all}
    else
       cat "${dir}/${fname_NTV_torque_int}" > ${fname_NTV_torque_int_all}
    fi
done
