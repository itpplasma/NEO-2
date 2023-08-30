#!/bin/bash

# Quick and dirty shell script to create references for a neo-2-par
# test.
# Runs all four stages and writes the md5sum to a file - a different one
# for each stage.

executable=neo_2.x
#~ runcommand="./$executable"
runcommand="mpirun -np 2 ./$executable"

rm *.dat
rm *.h5

#~ cp ../pcas1.dat ./

for numberofstage in 0 1 2 3 ; do
  echo "########################################################################"
  echo "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"
  switch_reconstruction.sh $numberofstage
  $runcommand
  md5sum *.dat | grep -v "MD5sums-reconstruct_" > MD5sums-reconstruct_$numberofstage
  echo "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"
  echo "########################################################################"
done
