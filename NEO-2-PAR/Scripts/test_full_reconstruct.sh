#!/bin/bash

########################################################################
### Variables and constants.

testcase=${1}
number_processors=${2}

testpath=`mktemp -d /temp/$LOGNAME/Neo2/Testing/Runs/Test-XXXXXX`
echo "Testpath created: ${testpath}"

referencepath='/temp/buchholz/Neo2/Testing/Reference/'
#~ referencepath='../../../Reference/'

totalnumberofstages=4
numberofstage=0
return_value=0

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

  # Specify path that should not be compared, because they are expected
  # to differ.
  # The /metadata object contains for example start and end times of the
  # run, which are expected to differ.
  exclude_paths="--exclude-path=/metadata --exclude-path=/Testcase1/NEO-2/neo2_config/metadata"
  # Specify the accuracy for the comparison. Our default is the default
  # of h5diff.
  accuracy=""

  # If the run is in parallel, some modifications to the settings have
  # to be made.
  if [ ${number_processors} -gt 1 ] ; then
    echo "##### parallel mode #####"
    # the object parallel_storage is indicator if the run was parallel
    # or not, as the reference should be made with a single processor,
    # this is expected to differ.
    exclude_paths="$exclude_paths --exclude-path=/Testcase1/NEO-2/taginfo/parallel_storage"
    # With finite precision order of summation plays a role, thus small
    # differences are to be expected. To account for this, we assume
    # values are equal, if the relative difference is small enough.
    accuracy="--relative=1.0e-13"
  fi

  for h5file in `ls $referencepath_local/${testcase_local}/*.h5` ; do
    if h5diff ${accuracy} ${exclude_paths} -q $h5file ./`basename $h5file` ; then
      a=''
    else
      echo "comparing $h5file and ./`basename $h5file`"
      echo "comparison command is 'h5diff ${accuracy} ${exclude_paths}'"
      h5diff ${exclude_paths} $h5file ./`basename $h5file`
      return_value_loc=1
    fi
  done

  return $return_value_loc
}

########################################################################
### Actual script

cp ./neo_2.x $testpath
cd $testpath

echo "Copying template..."
cp -r ../../Template/${testcase}/ .

cd ${testcase}
echo "Running Test ${testcase}..."

if [ ${number_processors} -eq 0 ]; then
  echo "Sequential mode"
  runcommand="$testpath/neo_2.x"
else
  echo "Parallel mode np=${number_processors}"
  runcommand="mpirun -np ${number_processors} $testpath/neo_2.x"
fi

for numberofstage in 0 1 2 3 ; do
  switch_reconstruction.sh $numberofstage
  $runcommand >> job.log 2>&1
  # check_equality_dat is a function.
  if check_equality_dat $referencepath ${testcase} $numberofstage ; then
    echo "Test ($numberofstage/$totalnumberofstages) passed. Checksums correct."
  else
    echo "Test ($numberofstage/$totalnumberofstages) failed. Not all files are equal."
    return_value=1
  fi
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
