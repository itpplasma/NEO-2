# NEO-2

[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/71e55b695f624fd9b5e18f97acdf95fd)](https://app.codacy.com/gh/itpplasma/NEO-2/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_coverage)

## Installation

### Prerequisites

git, cmake, make, ninja, gcc/g++/gfortran, MPI (OpenMPI or MPICH), BLAS/LAPACK (OpenBLAS or MKL), FFTW, and GSL. HDF5/NetCDF compatibility is provided through libneo and Fortio. Native CPU FortSparse is the default sparse solver; its tested revision is fetched during configuration. Optional AMD/METIS ordering libraries are permitted, while SuperLU, MUMPS and the GPL UMFPACK helper are disabled in this configuration. SuiteSparse is required only for the explicit legacy configuration below. For code coverage, lcov is also required. On Debian or Ubuntu run

    sudo apt install git cmake make ninja-build gcc g++ gfortran
    sudo apt install openmpi-bin openmpi-common libopenmpi-dev
    sudo apt install libopenblas-dev
    sudo apt install libsuitesparse-dev  # Optional ordering libraries / legacy provider
    sudo apt install libfftw3-dev libgsl-dev
    sudo apt install lcov  # Optional, for code coverage

### Build

Run

    make

inside the NEO-2 directory. This will handle the whole build process, including `cmake` calls.

To build against a specific libneo branch, tag, or commit SHA, pass `LIBNEO_REF` on the make command line:

    make LIBNEO_REF=<branch|tag|sha>

The equivalent cmake option is `-DLIBNEO_REF=<branch|tag|sha>`. When unset, the build uses the tested libneo revision pinned in `COMMON/ProjectConfig.cmake.in`. An ambient `LIBNEO_REF` in the shell does not affect the build.

The sparse provider is selected at configure time:

```sh
# Default: native CPU factorization, with no external sparse-solver linkage.
cmake --preset default -DNEO2_NATIVE_SPARSE=ON
# Explicit historical UMFPACK provider, for compatibility comparisons.
cmake --preset default -DNEO2_NATIVE_SPARSE=OFF
```

`FORTSPARSE_REF` selects the tested dependency revision; an explicit
`-Dfortsparse_SOURCE_DIR=/path/to/fortsparse` uses a local source checkout.
The native provider preserves `sparse_solve` and its real/complex, dense/COO/CSC
and multiple-RHS overloads. Operations 0/1/2/3 mean factor-solve-free,
factor-only, solve with retained factors (auto-factor if absent), and free.
Method 2 requests FortSparse's original-matrix iterative refinement; method 3
omits it. Changing method on a retained factor refactors the supplied matrix
to update the retained refinement policy. Numeric solver errors stop with a
diagnostic instead of returning an unchecked result.

Native sparse memory controls are in `/settings/` for both QL and PAR:

```fortran
&settings
  sparse_host_memory_limit_bytes = 2147483648 ! 2 GiB per MPI process; 0: unlimited
  sparse_low_memory = .true.                 ! RAM-only LOW_MEMORY schedule
/
```

Defaults are zero/unlimited and `.false.`, preserving the existing schedule.
The limit is a 64-bit count of **solver-owned host array payload**, including
factors, refinement, scratch and retained-factor overlap. It excludes NEO-2's
matrix assembly and dense/COO-to-CSC conversion, caller RHS/results, fixed
runtime descriptors, allocator overhead and process RSS. Leave separate node
memory for those costs and for every MPI process; this is not a whole-job limit.
NEO-2's operators use general real/complex LU, whose pivot-dependent fill has no
exact pre-factor memory forecast. The exact real-SPD planning API is not used or
claimed for these operators. A refused cap reports FortSparse RESOURCE_LIMIT
with required/available payload bytes; operating-system allocation failure is
a distinct status. No fallback silently ignores the requested cap.

A changed cap or LOW_MEMORY setting on a retained factor triggers refactorization
of the supplied matrix before the next solve, just as a method change does.
Freeing a factor remains available after changing settings. Both controls are
recorded in HDF5 configuration; the byte cap is a one-element int64 dataset,
so values above 2 GiB retain their exact integer value. Negative caps are rejected;
the explicit legacy provider accepts only the default memory controls.
A terminal native failure aborts the existing MPI job when MPI is active, so
peers cannot hang in scheduler communication. No collective is added to sparse
operations, and no distributed sparse solver is introduced.

All audited PAR, QL axi and Arnoldi calls use normal solves. Arnoldi's real
branch still splits complex RHS into real/imaginary solves against one real
factor; its complex branch uses a complex factor. There is one active shared
factor per process, as in the historical module. This adapter does not add
thread-safe concurrent handles, distributed factors or new transpose callers.
The spline LAPACK band path remains available. Native capabilities do not
establish large-application performance or full scientific-run admission.

Focused checks after building:

```sh
cmake --build --preset default --target test_native_sparse_adapter test_spline_sparse_assembly test_arnoldi_breakdown
ctest --test-dir build -R 'native_sparse_|spline_sparse_assembly_test|arnoldi_breakdown_test' --output-on-failure
```

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

## Run

NEO-2 requires Boozer files for MHD equilibria in text format (`.bc` ending). These can be generated via executables / functions in https://github.com/itpplasma/libneo/

* `efit_to_boozer` for tokamaks
* `vmec_to_boozer` for stellarators: https://github.com/itpplasma/libneo/blob/main/python/libneo/boozer.py
