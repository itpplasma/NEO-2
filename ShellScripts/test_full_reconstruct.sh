#!/bin/bash

########################################################################
### Variables and constants.

testcase=${1}
number_processors=${2}
which_code=${3}

testpath=`mktemp -d /temp/$LOGNAME/Neo2/Testing/Runs/Test-XXXXXX`
echo "Testpath created: ${testpath}"

referencepath='/temp/buchholz/Neo2/Testing/Reference/'
#~ referencepath='../../../Reference/'

executablename="neo_2_par.x"

totalnumberofstages=4
numberofstage=0
return_value=0

START=0
END=3

check_absolute_if_relative_fails=no #yes

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
  # Specify the accuracy for the comparison.
  exponent_def=-12
  accuracy="1.0e$exponent_def"

  # If the run is in parallel, some modifications to the settings have
  # to be made.
  if [ ${number_processors} -gt 1 ] ; then
    echo "##### parallel mode #####"
    # The object parallel_storage is indicator if the run was parallel
    # or not, as the reference should be made with a single processor,
    # this is expected to differ.
    exclude_paths="$exclude_paths --exclude-path=/Testcase1/NEO-2/taginfo/parallel_storage"
  fi

  for h5file in `ls $referencepath_local/${testcase_local}/*.h5` ; do
    # Note: return value of h5diff can not be used, as it can be 1, even
    # if no difference is reported.
    # Thus the use of grep -c, this will return 1 if the count is zero,
    # and return 0 if there is at least one occurence.
    if h5diff --relative=${accuracy} ${exclude_paths} $h5file ./`basename $h5file` | grep -c "difference" ; then
      echo "comparing $h5file and ./`basename $h5file`"
      echo "comparison command is 'h5diff --relative=${accuracy} ${exclude_paths}'"
      h5diff --relative=${accuracy} ${exclude_paths} $h5file ./`basename $h5file`

      exponent=$exponent_def
      another_round=yes
      while [ "x$another_round" == "xyes" -a $exponent -le 1 ]  ; do
        # \bug For some reason the line below does not work.
        h5diff --relative=1.0e$exponent ${exclude_paths} -q $h5file ./`basename $h5file` | grep -c "difference"
        if [ "x$?" = "x0" ] ; then
          #~ exponent=`let $exponent + 1`
          exponent=$[$exponent+1]
        else
          another_round=no
        fi
      done
      echo
      echo "Differences are up to order of ~$exponent."
      echo

      if [[ "x$check_absolute_if_relative_fails" == "xyes" ]] ; then
        # Make another check, this time with absolute differences.
        if h5diff --delta=${accuracy} ${exclude_paths} -q $h5file ./`basename $h5file` ; then
          echo "Second try succesfull."
        else
          echo "Pass with absolute difference also failed."
          return_value_loc=1
        fi
      else
        return_value_loc=1
      fi
    fi
  done

  return $return_value_loc
}

########################################################################
### Actual script

if [ "x${which_code}" = "xQL" ] ; then
  END=0
  executablename="neo_2_ql.x"
fi

cp ./$executablename $testpath
cd $testpath

echo "Copying template..."
cp -r -L ../../Template/${testcase}/ .

cd ${testcase}
echo "Running Test ${testcase}..."

if [ ${number_processors} -eq 0 ] ; then
  echo "Sequential mode"
  runcommand="$testpath/$executablename"
else
  echo "Parallel mode np=${number_processors}"
  runcommand="mpirun -np ${number_processors} $testpath/$executablename"
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
