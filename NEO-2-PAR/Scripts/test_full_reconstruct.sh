#!/bin/bash

testcase=${1}
number_processors=${2}

testpath=`mktemp -d /temp/$LOGNAME/Neo2/Testing/Runs/Test-XXXXXX`
echo "Testpath created: ${testpath}!"

referencepath='/temp/gernot_k/Neo2/Testing/Reference/'
#~ referencepath='../../../Reference/'

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

switch_reconstruction.sh 0
$runcommand >> job.log 2>&1

if md5sum -c $referencepath/${testcase}/MD5sums-reconstruct_0 ; then
  echo "Test (1/3) passed. Checksums correct."
else
  echo "Test (1/3) failed. Not all files are equal."
  exit 1
fi

switch_reconstruction.sh 1
$runcommand >> job.log 2>&1

if md5sum -c $referencepath/${testcase}/MD5sums-reconstruct_1 ; then
  echo "Test (2/3) passed. Checksums correct."
else
  echo "Test (2/3) failed. Not all files are equal."
  exit 1
fi

switch_reconstruction.sh 2
$runcommand >> job.log 2>&1

if md5sum -c $referencepath/${testcase}/MD5sums-reconstruct_2 ; then
  echo "Test (3/3) passed. Checksums correct."
else
  echo "Test (3/3) failed. Not all files are equal."
  exit 1
fi

echo "Tests passed. No exit code 1 received."
exit 0

#Maybe dangerous...
rm -r $testpath
