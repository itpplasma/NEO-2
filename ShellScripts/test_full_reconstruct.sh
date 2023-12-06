#!/bin/bash

########################################################################
### Variables and constants.

testcase=${1}
number_processors=${2}
which_code=${3}

# Check if the directory /temp/$LOGNAME/Neo2/Testing/Runs/ exists
if [ ! -d "/temp/$LOGNAME/Neo2/Testing/Runs/" ] && [ -z "$NEO2_TEST_PATH" ]; then
    echo "There is no directory defined for storing test data. Please create the /temp/$LOGNAME/Neo2/Testing/Runs/ directory or define a new location in the environment variable NEO2_TEST_PATH."
    exit 1
fi
# Check if NEO2_TEST_PATH is defined
if [ -n "$NEO2_TEST_PATH" ]; then
    # Check if the directory defined in NEO2_TEST_PATH exists
    if [ ! -d "$NEO2_TEST_PATH" ]; then
        echo "The directory defined in NEO2_TEST_PATH does not exist. Please define an existing directory."
        exit 1
    fi
fi
if [ -d "/temp/$LOGNAME/Neo2/Testing/Runs/" ]; then
    testpath=`mktemp -d /temp/$LOGNAME/Neo2/Testing/Runs/Test-${testcase}-XXXXXX`
else
    testpath=`mktemp -d $NEO2_TEST_PATH/Test-${testcase}-XXXXXX`
fi
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

function check_equality_md5sum {
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

function check_equality_dat {
  referencepath_local=${1}
  testcase_local=${2}

  return_value_loc=0

  abs_accuracy="1.0e-6"
  rel_accuracy="1.0e-5"

  for datfile in `ls $referencepath_local/${testcase_local}/*.dat` ; do
    testfile=`basename $datfile`

    #~ compare_data_files.m referencepath_local testcase_local abs_accuracy rel_accuracy
    compare_data_files.py $referencepath_local/${testcase_local}/$testfile $testfile $abs_accuracy $rel_accuracy
    res="$?"

    if [ ${return_value_loc} -eq 0 ] ; then
      return_value_loc="$res"
    fi
  done

  return $return_value_loc
}

# \brief Check if hdf5 files of testcase are equal.
#
# This function checks if all the files of the reference directory are
# equal (within a certain tolerance) to results of the testrun.
#
# Reference files are considered to be all hdf5 files, by extension .h5,
# that are in the reference directory.
# Over these files is iterated in a loop.
#
# The check itself is currently done via an inline python3 script, that
# calls an appropriate function from module hdf5tools, and exits with
# return value depending on output of the function.
#
# Note that at the moment only differences in data are considered. If
# the files differ in what fields are present, then this is ignored.
#
# input:
# ------
# referencepath: string with the path where the reference folders are
#   located.
# testcase: name of the testcase and thus of the subfolder.
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
    res="$?"
    # Avoid setting the return value to zero, if it was already unequal
    # zero.
    if [ ${return_value_loc} -eq 0 ] ; then
      return_value_loc="$res"
    fi
  done

  return $return_value_loc
}

########################################################################
### Actual script

cp ./$executablename ${testpath}/
cd ${testpath}/

echo "Copying template..."
# From the saved test templates location (hard coded)
cp -r -L /temp/AG-plasma/codes/neo-2_test_templates/${testcase}/* ./

echo "Running Test ${testcase}..."

# Note: absolute paths are used for the location of the executable.
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
  # check_equality_md5sum is a function.
  if check_equality_md5sum $referencepath ${testcase} $numberofstage ; then
    echo "Test ($numberofstage/$totalnumberofstages) passed. Checksums correct."
  else
    echo "Test ($numberofstage/$totalnumberofstages) failed. Not all files are equal."
    return_value=1
  fi
  numberofstage=$(($numberofstage+1))
done

if check_equality_dat ${referencepath} ${testcase} ; then
  echo "Test of dat files passed."
else
  echo "Test of dat files failed, there are differences or files which could not be compared."
  return_value=1
fi

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
