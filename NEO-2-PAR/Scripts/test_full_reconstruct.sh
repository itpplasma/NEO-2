#!/bin/bash

testcase=${1}
number_processors=${2}

testpath=`mktemp -d /temp/$LOGNAME/Neo2/Testing/Runs/Test-XXXXXX`
echo "Testpath created: ${testpath}!"

referencepath='/temp/buchholz/Neo2/Testing/Reference/'
#~ referencepath='../../../Reference/'

totalnumberofstages=4
numberofstage=0

########################################################################
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

########################################################################

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
  #~ if md5sum -c $referencepath/${testcase}/MD5sums-reconstruct_$numberofstage ; then
  # check_equality_dat is a function.
  if check_equality_dat $referencepath ${testcase} $numberofstage ; then
    echo "Test ($numberofstage/$totalnumberofstages) passed. Checksums correct."
  else
    echo "Test ($numberofstage/$totalnumberofstages) failed. Not all files are equal."
    exit 1
  fi
done

echo "Tests passed. No exit code 1 received."
exit 0

#Maybe dangerous...
rm -r $testpath
