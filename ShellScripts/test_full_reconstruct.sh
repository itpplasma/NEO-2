#!/bin/bash

########################################################################
### Variables and constants.

testcase=${1}
number_processors=${2}
which_code=${3}

testpath=`mktemp -d /temp/$LOGNAME/Neo2/Testing/Runs/Test-${testcase}-XXXXXX`
echo "Testpath created: ${testpath}"

referencepath='/temp/buchholz/Neo2/Testing/Reference/'
#~ referencepath='../../../Reference/'

executablename="neo_2.x"

totalnumberofstages=4
numberofstage=0
return_value=0

START=0
END=3

########################################################################
### Function definitions

function check_equality_dat {
  referencepath_local=${1}
  testcase_local=${2}
  numberofstage_local=${3}

  return_value_loc=0

  if md5sum -c $referencepath_local/${testcase_local}/MD5sums-reconstruct_$numberofstage_local ; then
    a='' # just to avoid reports of errors.
  else
    return_value_loc=1
  fi

  return $return_value_loc
}

function check_equality_hdf5 {
  referencepath_local=${1}
  testcase_local=${2}
  return_value_loc=0

  # Specify the accuracy for the comparison.
  exponent_def=-6
  accuracy="1.0e$exponent_def"

  # If the run is in parallel, some modifications to the settings have
  # to be made.
  if [ ${number_processors} -gt 1 ] ; then
    echo "##### parallel mode #####"
  fi

  for h5file in `ls $referencepath_local/${testcase_local}/*.h5` ; do
    testfile=`basename $h5file`

    echo "comparing $h5file and ${testfile}"

    echo "from hdf5tools import compare_hdf5_files; import sys; \
    res = compare_hdf5_files('$h5file', '${testfile}', ${accuracy}, [], '$referencepath_local/${testcase_local}/blacklist.txt', True); \
    sys.exit(0 if res[0] else 1)" | python3
    # Avoid setting the return value to zero, if it was already unequal
    # zero.
    if [ return_value_loc -eq 0 ] ; then
      return_value_loc="$?"
    fi
  done

  return $return_value_loc
}

########################################################################
### Actual script

cp ./$executablename $testpath
cd $testpath

echo "Copying template..."
cp -r -L ../../Template/${testcase}/* ./

echo "Running Test ${testcase}..."

if [ ${number_processors} -eq 0 ] ; then
  echo "Sequential mode"
  runcommand="$testpath/$executablename"
else
  echo "Parallel mode np=${number_processors}"
  runcommand="mpirun --oversubscribe -np ${number_processors} $testpath/$executablename"
fi

if [ "x${which_code}" = "xQL" ] ; then
  END=0
  totalnumberofstages=1
fi

numberofstage=$START
while [ $numberofstage -le $END ] ; do
  switch_reconstruction.sh $numberofstage
  $runcommand >> job.log 2>&1
  # check_equality_dat is a function.
  if check_equality_dat $referencepath ${testcase} $numberofstage ; then
    echo "Test ($numberofstage/$totalnumberofstages) passed. Checksums correct."
  else
    echo "Test ($numberofstage/$totalnumberofstages) failed. Not all files are equal."
    return_value=1
  fi
  numberofstage=$(($numberofstage+1))
done

if check_equality_hdf5 $referencepath ${testcase} ; then
  echo "Test of hdf5 files passed."
else
  echo "Test of hdf5 files failed, there are differences."
  return_value=1
fi

if [[ "x0" == "x$return_value" ]] ; then
  echo "Tests passed. No exit code 1 received."
  #Maybe dangerous...
  rm -r $testpath
else
  echo "At least one test failed."
  # Do not remove the folder, as one might want to check the log.
fi

exit $return_value
