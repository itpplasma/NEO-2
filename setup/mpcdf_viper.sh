#!/bin/bash

module load git cmake ninja 
module load gcc openmpi
module load mkl
module load hdf5-serial netcdf-serial fftw-serial
module load gsl/2.7

export CC=gcc
export CXX=g++
export FC=gfortran
