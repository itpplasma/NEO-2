# QL moment reduction: validation record

Date: 2026-11-01. Base: `8c8abdae67c3864a77c5f28ebcb07f95bb1c2c6f`.

## Environment and reproducibility

GNU Fortran 14.2.0 (Debian 14.2.0-19), CMake 3.31.6, Ninja, system BLAS/LAPACK,
`OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`. Release uses CMake's `-O3`; Debug uses `-g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`.
Both use `-ffree-line-length-none` for the pre-existing native routine.

The standalone CMake commands are in `ql_moment_reduction.md`. The local
validation workspace contains the new sources and the exact native
`integral_part` fragment read through the GitHub connector at the pinned base;
it is not a full NEO-2 checkout/build. The repository CMake extracts the same
routine from the full source without changing it.

SHA256 of the extracted routine, including its final newline:

    28c9b7f86e4b6aa8896633e95fabf0a889ec2b158d03be306f5ba51245cf11e9

## Results actually run

Both configurations: **2/2 CTest tests passed**, 204 collision/native-integration
assertions and 53 generic Schur assertions. No test tolerance or golden record
was changed to achieve these results.

| Quantity | Debug | Release |
|---|---:|---:|
| Native integral vs UCV, worst scaled max error over 42 comparisons | 7.152448e-18 | 6.938894e-18 |
| Native-kernel manufactured solve, reconstructed field max error | 2.482534e-16 | 2.224779e-16 |
| Native-kernel manufactured solve, original residual 2-norm | 3.158580e-15 | 3.051826e-15 |
| Native-kernel manufactured solve, GMRES steps | 5 | 5 |
| Independent full LAPACK vs generic reduced solve, max field error | 4.002966e-16 | 4.002966e-16 |
| Full/reduced eigenvalue matching, max absolute error | 4.343363e-17 | 2.564591e-17 |

The native-kernel solve has 270 full coordinates and 180 moment coordinates,
three species, and a manufactured diagonal sparse-inverse surrogate. It is not
a high-collisionality performance benchmark. The small eigenvalue test includes
a singular collision factor and verifies the extra full-space zeros separately.

Negative tests verify rejection of omitted regularizer terms, inaccurate
inverses, incompatible sources, wrong dimensions, nonfinite/failed callbacks,
exhausted iteration budgets, and separate-species energy overconstraints.
The energy test also solves the independently assembled full bordered matrix.

Four deliberately corrupted source variants were built separately in Debug;
all four were rejected by the owning tests, with the original sources untouched:

| Mutation | Observed failing assertion |
|---|---|
| Remove U's factor 1/2 | Native UCV equality |
| Remove the forward bounce contribution | Native UCV equality |
| Remove the original-equation acceptance check | Bordered success is not original-equation success |
| Replace the conjugate inner product by an unconjugated one | Complex reduced solve |

The two new modules also compile cleanly with `-std=f2008 -Wall -Wextra
-Wimplicit-interface -Wsurprising -fcheck=all` and no warnings. This statement
does not apply to untouched legacy sources.

## Not performed / no claim

- Full NEO-2 build or full repository `make test`.
- Production streaming matrices or equilibrium input decks.
- MPI launch, distributed collective validation or production solver selection.
- Existing golden-record runs and transport-observable comparisons.
- The n1e21 and lag-14 cluster cases mentioned in #247/#254.
- A wall-clock or memory speedup over current Arnoldi-Richardson.
- Improved conditioning, a cure for the sparse-inverse noise floor, or a
  Stafford-specific reduction theorem.

The implementation is kept off the production selection path until these gates
are addressed. Reported floating-point agreement verifies these fixtures; it
is not a general analytic conditioning or convergence bound.
