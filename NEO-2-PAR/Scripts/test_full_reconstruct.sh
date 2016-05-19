#!/bin/bash

testpath=`mktemp -d /temp/gernot_k/Neo2/Testing/Runs/Test-XXXXXX`
echo "Testpath created: ${testpath}!"
#testpath='/temp/gernot_k/Neo2/Testing/Runs/000'

cp ./neo_2.x $testpath
cd $testpath

echo "Copying template..."
cp -r ../../Template/${1}/ .

cd ${1}
echo "Running Test ${1}..."

if [ ${2} -eq 0 ]; then
	echo "Sequential mode"

	switch_reconstruction.sh 0
	../neo_2.x >> job.log 2>&1

	if md5sum -c ../../../Reference/${1}/MD5sums-reconstruct_0 ; then
        	echo "Test (1/3) passed. Checksums correct."
	        #exit 0
	else
       		echo "Test (1/3) failed. Not all files are equal."
        	exit 1
	fi

        switch_reconstruction.sh 1
        ../neo_2.x >> job.log 2>&1

        if md5sum -c ../../../Reference/${1}/MD5sums-reconstruct_1 ; then
                echo "Test (2/3) passed. Checksums correct."
                #exit 0
        else
                echo "Test (2/3) failed. Not all files are equal."
                exit 1
        fi

        switch_reconstruction.sh 2
        ../neo_2.x >> job.log 2>&1

        if md5sum -c ../../../Reference/${1}/MD5sums-reconstruct_2 ; then
                echo "Test (3/3) passed. Checksums correct."
                #exit 0
        else
                echo "Test (3/3) failed. Not all files are equal."
                exit 1
        fi

else
	echo "Parallel mode np=${2}"

        switch_reconstruction.sh 0
	mpirun -np ${2} ../neo_2.x >> job.log 2>&1

        if md5sum -c ../../../Reference/${1}/MD5sums-reconstruct_0 ; then
                echo "Test (1/3) passed. Checksums correct."
                #exit 0
        else
                echo "Test (1/3) failed. Not all files are equal."
                exit 1
        fi

        switch_reconstruction.sh 1
        ../neo_2.x >> job.log 2>&1

        if md5sum -c ../../../Reference/${1}/MD5sums-reconstruct_1 ; then
                echo "Test (2/3) passed. Checksums correct."
                #exit 0
        else
                echo "Test (2/3) failed. Not all files are equal."
                exit 1
        fi

        switch_reconstruction.sh 2
	mpirun -np ${2} ../neo_2.x >> job.log 2>&1

        if md5sum -c ../../../Reference/${1}/MD5sums-reconstruct_2 ; then
                echo "Test (3/3) passed. Checksums correct."
                #exit 0
        else
                echo "Test (3/3) failed. Not all files are equal."
                exit 1
        fi

fi

echo "Tests passed. No exit code 1 received."
exit 0

#Maybe dangerous...
#rm -r $testpath
