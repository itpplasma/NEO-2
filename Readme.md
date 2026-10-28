# NEO-2

NEO-2 solves the linearized drift kinetic equation to compute neoclassical
transport coefficients in tokamaks and stellarators, including the toroidal
torque (neoclassical toroidal viscosity) caused by non-axisymmetric magnetic
perturbations. It is developed at the Institute of Theoretical and
Computational Physics, TU Graz.

[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/71e55b695f624fd9b5e18f97acdf95fd)](https://app.codacy.com/gh/itpplasma/NEO-2/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_coverage)

The repository contains two variants of the code:

* NEO-2-QL (`neo_2_ql.x`): tokamaks with small non-axisymmetric perturbations,
  multispecies plasmas, toroidal torque.
* NEO-2-PAR (`neo_2_par.x`): stellarators and unperturbed tokamaks, with
  relativistic options (e.g. for current drive calculations), MPI-parallel.

Python tools for preparing runs and reading output are in `python/`.

## Installation

### Prerequisites

git, cmake, make, ninja, gcc/g++/gfortran, MPI (OpenMPI or MPICH), BLAS/LAPACK (OpenBLAS or MKL), SuiteSparse, FFTW, and GSL. HDF5/NetCDF compatibility is provided through libneo and Fortio. If SuiteSparse is not available on your system, it will be built automatically by NEO-2. For code coverage, lcov is also required. On Debian or Ubuntu run

    sudo apt install git cmake make ninja-build gcc g++ gfortran
    sudo apt install openmpi-bin openmpi-common libopenmpi-dev
    sudo apt install libopenblas-dev libsuitesparse-dev
    sudo apt install libfftw3-dev libgsl-dev
    sudo apt install lcov  # Optional, for code coverage

If you have additional compilers and MPI implementations (e.g. Intel) installed, be sure to set the environment like

    export CC=gcc
    export CXX=g++
    export FC=gfortran
    export MPI_HOME=/usr

before the build. You may also need to `deactivate` Python environments such as Anaconda.

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

More details on dependencies and cluster builds are in [Build.md](Build.md).

## Testing

    make test

runs the CTest suite. `make coverage` rebuilds with coverage instrumentation,
runs the tests and prints a summary; see [CONTRIBUTING.md](CONTRIBUTING.md) for
the HTML report.

## Run

NEO-2 requires Boozer files for MHD equilibria in text format (`.bc` ending). These can be generated via executables / functions in https://github.com/itpplasma/libneo/

* `efit_to_boozer` for tokamaks
* `vmec_to_boozer` for stellarators: https://github.com/itpplasma/libneo/blob/main/python/libneo/boozer.py

Annotated example input files are in `DOC/`: `neo2.in.ql-full` and
`neo2.in.par-full` (NEO-2 settings), `neo.in.ql-full` and `neo.in.par-full`
(magnetic field settings), and the example Boozer files `example_bg_field.bc`
and `example_pert_field.bc`.

## Documentation

* `DOC/UserDoc/` — user manual (LaTeX): building, input, output, running,
  theory. Build it with `make doc`.
* `DOC/ExtraDocuments/` — worked run descriptions for NEO-2-PAR (ECCD) and
  NEO-2-QL (NTV).
* `DOC/NEO-2-QL_Quick_Start_Manual/` — quick start manual for NEO-2-QL.
* `DOC/ripple_solver_normalizations_and_output.tex` — normalizations and
  output quantities of the ripple solver.

## Citing

If you use NEO-2 in published work, please cite

* W. Kernbichler, S. V. Kasilov, G. Kapper, A. F. Martitsch, V. V. Nemov,
  C. G. Albert and M. F. Heyn, Plasma Phys. Control. Fusion 58, 104001 (2016),
  https://doi.org/10.1088/0741-3335/58/10/104001

and, for the toroidal torque computed by NEO-2-QL, S. V. Kasilov et al.,
Phys. Plasmas 21, 092506 (2014).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for the review and merge workflow.

## License

MIT, see [LICENSE](LICENSE).
