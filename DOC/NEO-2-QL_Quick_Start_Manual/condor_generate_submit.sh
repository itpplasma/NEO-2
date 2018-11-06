#!/bin/bash

echo -e "\nExecutable:       $0"
echo -e "Creating submit file for Condor ...\n"

if [ -e submit ]
then
    rm submit
fi

if [ \( "$1" = "--help" \) ]
then
    echo "Use following syntax to create a valid submit file for Condor:"
    echo -e "$0\n"
else
    dir_list=$(find . -name "es_*" -type d -printf "%f \n")

    # settings for all subsequent runs
    echo -e "Executable = exec_collop_precomp.sh" > "submit"
    echo -e "Universe   = vanilla" >> "submit"
    echo -e "Error      = foo.err" >> "submit"
    echo -e "Log        = foo.log" >> "submit"
    #echo -e "transfer_input_files = ${TRANSFERED_INPUT_FILES}" >> "submit"
    echo -e "run_as_owner = true" >> "submit"
    echo -e "notify_user  = Andreas.Frank.Martitsch@itp.tugraz.at" >> "submit"
    echo -e "notification = Error" >> "submit"
    echo -e "nice_user    = false" >> "submit"
    #echo -e "nice_user    = true" >> "submit"
    #echo -e "requirements = ( TARGET.UtsnameNodename != \"faepop29\" )" >> "submit"
    echo -e "requirements = ( TARGET.OpSysMajorVer == 8 )" >> "submit"
    echo -e "request_cpus = 4" >> "submit"
    echo -e "request_memory = 6*1024" >> "submit"
    #echo -e "rank = Memory\n" >> "submit"
    echo -e "Getenv = true\n" >> "submit"

    for dir in ${dir_list}; 
    do
	# echo -e ${dir}
	echo -e "Output     = out" >> "submit"
	echo -e "initialdir = ${dir}" >> "submit"
	echo -e "Queue\n" >> "submit"
    done
fi
