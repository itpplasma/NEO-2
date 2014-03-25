#!/bin/bash

testpath=`mktemp -d /temp/gernot_k/Neo2/Testing/Runs/Test-XXXXXXXX`
#testpath='/temp/gernot_k/Neo2/Testing/Runs/000'

cp ./neo_2.x $testpath
cd $testpath

echo "Copying template..."
cp -r ../../Template/${1}/ .

cd ${1}
echo "Running Test ${1}..."

if [ ${2} -eq 0 ]; then
	echo "Sequential mode"
	../neo_2.x >> job.log 2>&1
else
	echo "Parallel mode np=${2}"
	mpirun -np ${2} ../neo_2.x >> job.log 2>&1
fi

if md5sum -c ../../../Reference/${1}/MD5sums ; then
	echo "Test passed. Checksums correct."
	exit 0
else
	echo "Test failed. Not all files are equal."
	exit 1
fi

#Maybe dangerous...
#rm -r $testpath
