#!/bin/bash

# Manual definition via Condor submit file
CPUs=$1

echo "Job requested $CPUs CPUs, starting mpirun..."
sed -i "s#prop_reconstruct.*=.*#prop_reconstruct = 0#" neo2.in
echo "Reconstruct 0 - Calculate Greensfunctions of ripples"
OMP_NUM_THREADS=2 mpirun.openmpi -mca orte_tmpdir_base "/tmp/" -np $CPUs ./neo_2.x
sed -i "s#prop_reconstruct.*=.*#prop_reconstruct = 1#" neo2.in
echo "Reconstruct 1 - Solve for boundary condition on fieldline start"
OMP_NUM_THREADS=2 ./neo_2.x
sed -i "s#prop_reconstruct.*=.*#prop_reconstruct = 2#" neo2.in
echo "Reconstruct 2 - Calculate the solution g(eta,phi) at the interpolation grid knots."
OMP_NUM_THREADS=2 mpirun.openmpi -mca orte_tmpdir_base "/tmp/" -np $CPUs ./neo_2.x
sed -i "s#prop_reconstruct.*=.*#prop_reconstruct = 3#" neo2.in
echo "Reconstruct 3 - Save results to .h5 and remove individual files"
OMP_NUM_THREADS=2 ./neo_2.x

echo "Job done."
