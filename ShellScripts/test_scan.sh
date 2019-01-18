#!/bin/bash

########################################################################
### The following variables/constant need to be defined (e.g. by .bashrc):
# NEO2PATH

########################################################################
### Variables and constants.

testcase=${1}
which_code=${2}
testpath=`mktemp -d /temp/$LOGNAME/Neo2/Testing/Runs/Test-XXXXXX`
#echo "Testpath created: ${testpath}"

executablename="neo_2.x"
sourcefoldername="NEO-2-PAR"
buildfoldername="Build-Release"

filetocollect="fulltransp.h5"
collectfilename="fulltransp_collected.h5"
comparefilename="fulltransp.h5"

referencepath="/temp/buchholz/Neo2/Testing/Reference/"

return_value="-1"

precomp_dir="precomp"
comp_dir="comp"
program_set_second_stage="set_neo2in_fullcalc.sh"
program_data_collection="python3 <<EOF
import hdf5tools
hdf5tools.copy_hdf5_from_subfolders_to_single_file('./', '$filetocollect', '$collectfilename')
quit()
EOF"

program_data_reshaping="python3 <<EOF
import hdf5tools
hdf5tools.reshape_hdf5_file('$collectfilename', '$comparefilename')
quit()
EOF"

########################################################################
### Function definitions

function wait_till_condor_jobs_are_finished {
  id_scan_local="$1"

  runs_found="1"
  echo "Waiting for jobs with id $id_scan_local to finish."
  while [ $runs_found -gt 0 ] ; do
    sleep 60
    runs_found=`condor_q --nobatch | grep -c "${id_scan_local}."`
#    echo "$id_scan_local $runs_found"
    if [ $runs_found -gt 0 ] ; then
      echo "still $runs_found job(s) ...continue"
    else
      echo "no more jobs left ..finished"
    fi
  done
}

function compare_data {
  referencepath_local=${1}
  testcase_local=${2}
  return_value_loc=0

  objects_to_compare="/gamma_out"
  exponent_def=-12
  accuracy="1.0e$exponent_def"

  for h5file in $comparefilename ; do
    if ! h5diff --relative=${accuracy} -q $referencepath_local/$testcase_local/$h5file $h5file $objects_to_compare ; then
      echo "Inital comparison failed."
      return_value_loc=1
      exponent=$exponent_def
      another_round=yes
      while [ "x$another_round" == "xyes" -a $exponent -le 1 ]  ; do
        h5diff --relative=1.0e$exponent -q $h5file ./`basename $h5file` $objects_to_compare
        if [ "x$?" = "x0" ] ; then
          another_round=no
        else
          exponent=$[$exponent+1]
        fi
      done
      echo
      echo "Differences are up to order of ~$exponent."
      echo
    else
      echo "objects $objects_to_compare of files '$h5file' are equal."
    fi
  done

  exit $return_value_loc
}

########################################################################
### Actual script

cd $testpath

mkdir $precomp_dir

cd $precomp_dir

if [ "x${which_code}" = "xQL" -o "x${which_code}" = "xql" ] ; then
  sourcefoldername="NEO-2-QL"
fi

# Copy executable
cp $NEO2PATH/$sourcefoldername/$buildfoldername/$executablename ./

# Copy also the template files.
cp -L ../../../Template/${testcase}/* ./

# execute neo-2 to create folders
./$executablename

echo "submitting precomp jobs..."

# submit jobs to condor
condor_out=`condor_submit submit_precomp`
id_scan=`echo "$condor_out" | grep "submitted to cluster" | sed s/"^[0-9]* job(s) submitted to cluster "// | sed s/"\."//`

echo "id $id_scan...done"

# wait till jobs are finished.
wait_till_condor_jobs_are_finished $id_scan

echo "precomp finished"

# copy files for second stage
cd ..
cp -R $precomp_dir $comp_dir/
cd $comp_dir

# run script to change input for second stage.
$program_set_second_stage

echo "submitting comp jobs..."

# submit jobs to condor
condor_out=`condor_submit submit`
id_scan=`echo "$condor_out" | grep "submitted to cluster" | sed s/"^[0-9]* job(s) submitted to cluster "// | sed s/"\."//`

echo "id $id_scan...done"

# wait till jobs are finished
wait_till_condor_jobs_are_finished $id_scan

echo "comp finished"

#~ # collect data
#~ #$program_data_collection
python3 <<EOF
import hdf5tools
hdf5tools.copy_hdf5_from_subfolders_to_single_file('./', '$filetocollect', '$collectfilename')
quit()
EOF

echo "data collected"

# reshape data to make comparison easier? each group would need to be given explicitly in the original form.
#$program_data_reshaping
python3 <<EOF
import hdf5tools
hdf5tools.reshape_hdf5_file('$collectfilename', '$comparefilename')
quit()
EOF

echo "data reshaped"

# compare data or part of it, to reference.
compare_data $referencepath $testcase # function call

return_value=$?

echo "data comparison finished"

cd ../..

echo "test scan finished"

exit $return_value
