# NEO-2

[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/71e55b695f624fd9b5e18f97acdf95fd)](https://app.codacy.com/gh/itpplasma/NEO-2/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_coverage)

# Installation

## Prerequisites
git, cmake, make, ninja, gcc/g++/gfortran, MPI (OpenMPI or MPICH), BLAS/LAPACK (OpenBLAS or MKL), SuiteSparse, FFTW, GSL, HDF5/NetCDF with Fortran libraries. If SuiteSparse is not available on your system, it will be built automatically by NEO-2. For code coverage, lcov is also required. On Debian or Ubuntu run

    sudo apt install git cmake make ninja-build gcc g++ gfortran
    sudo apt install openmpi-bin openmpi-common libopenmpi-dev
    sudo apt install libopenblas-dev libsuitesparse-dev
    sudo apt install libfftw3-dev libgsl-dev libhdf5-dev libnetcdf-dev libnetcdff-dev
    sudo apt install lcov  # Optional, for code coverage

## Build

Run

    make

inside the NEO-2 directory. This will handle the whole build process, including `cmake` calls.

You obtain a `build` directory with subdirectories

* `NEO-2-QL` with binary `neo_2_ql.x` for the tokamak version.
* `NEO-2-PAR` with binary `neo_2_par.x` for the stellarator version.

## Testing

To run tests:

```bash
    make test
```

To generate code coverage report:

```bash
    make coverage
```

This will build with coverage instrumentation, run all tests, and display a coverage summary. Coverage data files are generated in the `build` directory.

To view detailed coverage report:

    # After running make coverage, view the summary
    lcov --summary build/coverage_filtered.info
    
    # Generate HTML report for detailed viewing
    cd build && genhtml coverage_filtered.info --output-directory coverage_html
    # Open build/coverage_html/index.html in a web browser

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
