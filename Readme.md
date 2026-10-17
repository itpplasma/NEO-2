# NEO-2

NEO-2-QL supports named Coulomb collision modes and experimental physical
Krook damping. See [collision model inputs and normalization](DOC/collision_models.md).

[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/71e55b695f624fd9b5e18f97acdf95fd)](https://app.codacy.com/gh/itpplasma/NEO-2/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_coverage)

## Installation

### Prerequisites

git, cmake, make, ninja, gcc/g++/gfortran, MPI (OpenMPI or MPICH), BLAS/LAPACK (OpenBLAS or MKL), SuiteSparse, FFTW, and GSL. HDF5/NetCDF compatibility is provided through libneo and Fortio. If SuiteSparse is not available on your system, it will be built automatically by NEO-2. For code coverage, lcov is also required. On Debian or Ubuntu run

    sudo apt install git cmake make ninja-build gcc g++ gfortran
    sudo apt install openmpi-bin openmpi-common libopenmpi-dev
    sudo apt install libopenblas-dev libsuitesparse-dev
    sudo apt install libfftw3-dev libgsl-dev
    sudo apt install lcov  # Optional, for code coverage

### Build

Run

    make

inside the NEO-2 directory. This will handle the whole build process, including `cmake` calls.

To build against a specific libneo branch, tag, or commit SHA, pass `LIBNEO_REF` on the make command line:

    make LIBNEO_REF=<branch|tag|sha>

The equivalent cmake option is `-DLIBNEO_REF=<branch|tag|sha>`. When unset, the build uses the tested libneo revision pinned in `COMMON/ProjectConfig.cmake.in`. An ambient `LIBNEO_REF` in the shell does not affect the build.

You obtain a `build` directory with subdirectories

* `NEO-2-QL` with binary `neo_2_ql.x` for the tokamak version.
* `NEO-2-PAR` with binary `neo_2_par.x` for the stellarator version.

## Testing

To run tests:

    make test

To generate code coverage report:

    make coverage

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

## Merging

- Pull requests, including every layer of a stacked PR, run the fast CI tier.
  Label a PR `full-ci` to run the full tier on it before queueing.
- Review a stack layer by layer. When it is ready, collapse it into one PR into
  `main` (retarget the top layer to `main`); the lower PRs then show as merged.
- Click "Merge when ready" (or `gh pr merge --merge --auto`). The PR enters the
  merge queue, the full CI tier runs on the exact merge result, and on success
  it lands on `main` as a merge commit. `main` itself is not re-tested on push.
- A branch need not be up to date with `main`. Merge `main` into a PR only
  when it conflicts.

## Run

NEO-2 requires Boozer files for MHD equilibria in text format (`.bc` ending). These can be generated via executables / functions in https://github.com/itpplasma/libneo/

* `efit_to_boozer` for tokamaks
* `vmec_to_boozer` for stellarators: https://github.com/itpplasma/libneo/blob/main/python/libneo/boozer.py
