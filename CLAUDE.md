# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NEO-2 is a Fortran-based scientific computing project for neoclassical transport calculations in fusion plasma devices. It has two main variants:
- **NEO-2-QL**: For tokamak (axisymmetric) configurations
- **NEO-2-PAR**: For stellarator (3D) configurations with parallel computation support

## Build Commands

```bash
# Full build (creates both neo_2_ql.x and neo_2_par.x)
make

# Run tests
make test

# Generate documentation
make doc

# Clean build directory
make clean
```

The build system uses CMake with Ninja as the generator. Executables are created in:
- `build/NEO-2-QL/neo_2_ql.x`
- `build/NEO-2-PAR/neo_2_par.x`

## Build System Notes

- We have an outer Makefile on the project top level that handles cmake magic
- Just run make to build and make test to test and make coverage to generate coverage
- Never run make inside build directory, just use ninja but prefer outer make

## Testing

Golden record tests are run via GitHub Actions on PRs. To run tests locally:
```bash
make test
```

For PR submissions, the CI runs:
- Lorentz and QL tests on all pushes
- PAR tests only on pull requests (due to longer runtime)

## Code Architecture

### Core Components
- `/COMMON/`: Shared Fortran modules for both variants (magnetic field calculations, collision operators, splines)
- `/NEO-2-QL/`: Tokamak-specific implementation
- `/NEO-2-PAR/`: Stellarator-specific implementation with MPI parallelization

### Python Package (`/python/`)
Structured package with modules for:
- `neo2_mars/`: MARS integration tools
- `neo2_ql/`: QL-specific utilities
- `neo2_util/`: General utilities
- `neo2par/`: PAR-specific utilities

Install with: `cd python && pip install -e .`

### Key Dependencies
- **Fortran**: Primary language (F90/95/2008)
- **MPI**: For parallel computation
- **Libraries**: BLAS/LAPACK, SuiteSparse (UMFPACK), GSL/FGSL, HDF5, NetCDF, FFTW3
- **Python**: numpy, h5py, omfit_classes, libneo

### Input/Output
- Input: Boozer coordinate files (`.bc`) for MHD equilibria
- Output: HDF5 files with transport coefficients and physics quantities

## Development Guidelines

- Follow Fortran 90+ standards, avoid Fortran 77 style
- Use existing modules in `/COMMON/` for shared functionality
- Maintain separation between QL and PAR implementations
- Output data in HDF5 format for efficiency
- Keep MPI code isolated to PAR variant

## Testing Guidelines

- Tests must be non-shallow, non-tautological, and as small and fast as possible
- Avoid any file or network I/O wherever possible