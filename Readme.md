# Installation

## Prerequisites
git, cmake, make, ninja, gcc/g++/gfortran, MPI (OpenMPI or MPICH), BLAS/LAPACK (OpenBLAS or MKL), SuiteSparse, FFTW, GSL, HDF5/NetCDF with Fortran libraries. If SuiteSparse is not available on your system, it will be built automatically by NEO-2. On Debian or Ubuntu run

    sudo apt install git cmake make ninja-build gcc g++ gfortran
    sudo apt install openmpi-bin openmpi-common libopenmpi-dev
    sudo apt install libopenblas-dev libsuitesparse-dev
    sudo apt install libfftw3-dev libgsl-dev libhdf5-dev libnetcdf-dev libnetcdff-dev

## Build
Run

    make

inside the NEO-2 directory. This will handle the whole build process, including `cmake` calls.

You obtain a `build` directory with subdirectories

* `NEO-2-QL` with binary `neo_2_ql.x` for the tokamak version.
* `NEO-2-PAR` with binary `neo_2_par.x` for the stellarator version.

If you have additional compilers and MPI implementations (e.g. Intel) installed, be sure to set the environment like

    export CC=gcc
    export CXX=g++
    export FC=gfortran
    export MPI_HOME=/usr

before the build. You may also need to `deactivate` Python environments such as Anaconda.

# Run

NEO-2 requires Boozer files for MHD equilibria in text format (`.bc` ending). These can be generated via executables / functions in https://github.com/itpplasma/libneo/

* `efit_to_boozer` for tokamaks
* `vmec_to_boozer` for stellarators: https://github.com/itpplasma/libneo/blob/main/python/libneo/boozer.py
