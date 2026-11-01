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

## Opt-in performance microbenchmark (synthetic, 2026-11-01)

Run `cmake --build build/moment-release --target bench_ql_moment_vs_full`, then
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 build/moment-release/bench_ql_moment_vs_full`.
The benchmark compares the same complex restarted GMRES implementation in two
coordinates: full-space `W=I, D=UCV` and moment-space `W=V, D=UC`. Both use
**the same** manufactured full original operator `L=P-UCV`, forward-substitution
`P^-1` from a bidiagonal streaming surrogate, collision Legendre deposition
patterns (including parity and bounce), two species and manufactured coefficients.
The original equation and complete reconstructed vector are checked after every
solve. The required original residual is `2e-10*||q||_2`; restart is 20 and the
maximum is 120 Krylov iterations. Eight timed samples follow a warmup; only
solve time is measured (no model or coefficient setup or sparse factorization).

Measured on the available container AMD EPYC 9V74 virtual CPU, GNU Fortran 14.2,
Release `-O3`, one BLAS/OpenMP thread, process affined to CPU 2:

| Full dimensions | Moment dimensions | Full median CPU | Moment median CPU | Speedup | Full/moment Krylov steps | Original residual full/moment |
|---:|---:|---:|---:|---:|---:|---:|
| 9,792 | 1,152 | 22.832 ms | 9.932 ms | 2.30x | 17 / 17 | 3.22e-9 / 2.66e-9 |
| 52,224 | 2,048 | 424.876 ms | 111.986 ms | 3.79x | 53 / 56 | 5.59e-8 / 1.47e-8 |
| 145,600 | 3,200 | 1,084.468 ms | 265.901 ms | 4.08x | 50 / 54 | 6.19e-8 / 1.18e-8 |

Maximum full-versus-moment reconstructed vector differences were 9.3e-11,
3.4e-9 and 5.6e-9 respectively. The full-space and moment-space algorithms
satisfy the same original-equation tolerance; the reduced solver used 0, 3
and 4 more Krylov steps respectively, but fewer and shorter global-vector
orthogonalizations. Timings vary with CPU scheduling. Earlier repetitions
showed approximately 2.3–2.9x, 3.5–3.9x, and 3.9–4.2x respectively.

**Limitations.** This is not the production Arnoldi-Richardson baseline, nor
production NEO-2's streaming matrix, actual collision coefficients, MPI,
physical inputs, n1e21 or lag-14 results. In an additional uncommitted stress
setting (gain 1.7 for the largest case), *both* full and moment GMRES exhausted
120 iterations and had O(10) original residual norms. Moment reduction alone
does not cure a highly unfavorable spectrum or inverse noise floor.
No production speedup, wall-time improvement over native NEO-2, or Stafford-
specific acceleration is established.
