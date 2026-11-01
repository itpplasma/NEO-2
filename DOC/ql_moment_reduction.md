# Exact QL collision-moment reduction

## Status and scope

This is a **tested serial reference implementation and native-kernel regression**.
It is not yet a selectable production QL solver. No `neo2.in` setting, existing
`integral_part` implementation, Arnoldi routine, conservation filter, sparse
provider, or golden record is changed. The new modules are built with the QL
library and the tests are registered with CTest.

The production collision routine is exercised verbatim in a test adapter. Its
host data are manufactured and the species allgather is emulated. This tests the
retained collision arithmetic, not MPI communication, a production equilibrium,
or the full sparse streaming operator. Full native QL qualification, distributed
integration, and high-collisionality timing remain open.

This is classical Schur elimination implemented for NEO-2's actual moment
factorization, **not a Stafford-specific theorem or demonstrated speedup**.

## 1. The exact retained collision action

`NEO-2-QL/collision_moment_map.f90` exposes three maps:

* `map%project`: V, the parity-aware Legendre projection (`pleg_bra`).
* `mix_collision_moments`: C, the retained energy/species matrix (`ailmm_aa`).
* `map%lift`: U, Legendre reconstruction (`pleg_ket`) **and spatial deposition**.

Thus the current discrete integral collision action is `C_I = U C V`.
U includes the factor 1/2, forward and backward deposition, unequal neighboring
pitch supports, endpoint treatment, and `addboucol` bounce-boundary terms.
It is not just a pointwise Legendre reconstruction. No new truncation is made.
C may be singular; neither C nor VU is inverted or assumed to be the identity.

The map takes one species' distribution and moment vector. The mixing operation
takes the gathered moments of all species. Arrays are one-based externally:
`coeff(mout+1,min+1,l+1,sout+1,sin+1)`. Moments are flattened in Fortran order,
with energy fastest, then Legendre index, then field-line position. The
constructor builds the same contiguous block offsets as `ind_start` and returns
an error for malformed dimensions or nonfinite coefficients.

For a common field-line mesh, each species has

    n = 2 (lag+1) sum_phi (npassing(phi)+1)
    r = (lag+1) (leg+1) N_phi.

There is no unconditional claim that r<n. Both dimensions and the cost of each
inverse application must be measured for an actual input. Different species may
use different pitch maps, but the current mixing interface requires a common
retained energy/Legendre/field-line moment layout.

## 2. Reduction, reconstruction, and exactness proof

Let the actual original equation be `L f = q`, and suppose an invertible
preconditioning matrix P and maps D,W satisfy

    L = P - D W.

The simplest nonaxisymmetric case takes `D=UC`, `W=V`. Set

    f0 = P^-1 q
    m  = W f
    S  = I - W P^-1 D.

Solve `S m = W f0`, then reconstruct `f = f0 + P^-1 D m`.
Indeed, the reduced equation gives Wf=m, and substitution gives Lf=q.
Conversely, every solution of Lf=q yields this same reduced equation and
reconstruction. This is a bijection even if L or S is singular; uniqueness
requires additional hypotheses. An iterative solver need not converge on a
singular or inconsistent system, and must report that fact.

The implementation applies S matrix-free and never forms an n-by-n matrix or
stores the dense response P^-1 D. GMRES basis storage is O(r * restart), in
addition to O(n) reconstruction work and the externally retained factorization.
The reference implementation has global vectors and ordinary complex inner
products; it must NOT be passed rank-local vectors as though it were an MPI
solver.

### The axisymmetric regularizer cannot be dropped

For the AX branch, the factorized `amat_sp` is generally P=M+R, whereas the
original equation is L=M-UCV. Therefore

    L = P - (UCV + R),

not P-UCV. If R=EF, use

    D = [UC, E],       W = [V; F].

A matrix supported on s rows always has such a factorization: E injects those
rows and F contains their full row values. This is an exact row-support bound,
not a claim that the particular regularizer is small or that this factorization
is rank-minimal. The additional regularizer coordinates count toward r.
The tests deliberately omit R from D W while keeping the original L callback;
the independent residual then rejects the result.

### Spectrum: preserved, not improved

For T=P^-1 D W and Tm=W P^-1 D, the nonzero eigenvalues, with algebraic
multiplicity, agree. For n>=r,

    det(z I_n - T) = z^(n-r) det(z I_r - Tm).

This follows from the block determinant identity det(I-AB)=det(I-BA), and
extends from nonzero z to a polynomial identity. For n<r the roles reverse.
Moment-space reduction removes extra zero modes; it does not remove the
troublesome nonzero spectrum. It can reduce Krylov storage, but does not imply
better conditioning, normality, iteration counts, or a more accurate P^-1.

### Small reduced residual is insufficient

With exact reconstruction and `s = m-Wf`,

    Lf-q = D s.

A small s can be amplified by D; numerical inverse errors add another defect.
The implementation accepts only a freshly evaluated **original** equation
residual, supplied independently of the reduction. An intentionally inaccurate
inverse in the tests converges in reduced space but is correctly rejected after
reconstruction. This is relevant to the inverse-noise issue reported in #247;
this reduction does not claim to cure that issue.

## 3. Conservation and bordered systems

Optional caller-supplied constraints are handled by exact elimination of

    [L K; B 0] [f;lambda] = [q;g].

Here B is the constraint action, including any desired transpose/conjugation;
the solver does not implicitly transpose it. Define X=P^-1 D and G=P^-1 K.
The reduced system is

    [I-WX   WG ] [m     ] = [W f0    ]
    [ BX   -BG ] [lambda]   [g-B f0  ],

with reconstruction `f=f0+X m-G lambda`. The code applies this matrix-free.
It never divides by a scalar `1-Ritz_value`, so a compatible unit feedback mode
can be retained and constrained rather than numerically divided out.

The solver checks all three quantities separately:

    ||Lf-q||_2,  ||Lf+K lambda-q||_2,  ||Bf-g||_2.

Success requires both physical and bordered equation residuals below
`atol+rtol*||q||_2`, and the constraint residual below `constraint_atol`.
Constraint rows must be physically normalized and that absolute tolerance must
be chosen accordingly. A nonzero K lambda is not silently declared an original
kinetic solution. `schur_original_residual` identifies the case where the
bordered equations pass but the original equation does not.

The caller must provide the *actual* conserved combinations and source
compatibility. This implementation does not infer them from a collision-model
name. The two-species energy-exchange test retains the difference mode while
constraining total energy; separately constraining both species produces
nonzero multipliers and is rejected as a solution of the original equation.
This is a manufactured algebraic test, not a diagnosis of a particular NEO-2
input deck. Near-null vectors are not automatically exact invariants.

## 4. API and error contract

`moment_schur_operator_t` contains dimensions and callbacks:

    inverse_p(x,y,info)    P^-1 x (retained factorization, not a changing tolerance)
    feedback(m,y,info)     D m
    moments(f,m,info)      W f
    original(f,y,info)     L f, before inverse/projection/filtering
    force(lambda,y,info)  K lambda (only when constraints exist)
    constraints(f,y,info) B f      (only when constraints exist).

Callbacks act on global concatenated complex(real64) vectors, set every output,
and return zero on success. The driver invokes

    call solve_moment_schur(op, q, g, f, lambda, options, report)

with zero-length g and lambda for the unconstrained problem. All output fields
are initialized; never use a returned iterate as a converged result unless
`report%status == schur_ok`. Finite nonconverged iterates are returned for
inspection. Counters distinguish Krylov steps, reduced operator applications
(including explicit residual checks), and inverse calls.

The restarted complex GMRES uses two-pass modified Gram-Schmidt, LAPACK ZLARTG,
true reduced residuals at every restart, explicit breakdown detection, and a
finite iteration budget. It is not the legacy Arnoldi/Ritz deflation algorithm.
No dimension-dependent finite-termination guarantee is claimed in floating
point arithmetic. Each inverse remains a full kinetic solve; this formulation
does not remove the sparse factorization cost.

Statuses are: success (0), invalid input (1), callback failure (2), nonfinite
arithmetic (3), Krylov breakdown (4), exhausted budget (5), modified original
equation (6), constraint residual failure (7), allocation failure (8), and
reconstruction/bordered residual failure (9). The no-feedback and zero-source
cases are handled explicitly without normalizing a zero seed.

## 5. Reproduce the bounded validation

From the repository root:

```sh
cmake -S TEST/moment_reduction -B build/moment-debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build/moment-debug --parallel
ctest --test-dir build/moment-debug -V
cmake -S TEST/moment_reduction -B build/moment-release -DCMAKE_BUILD_TYPE=Release
cmake --build build/moment-release --parallel
ctest --test-dir build/moment-release -V
```

Only a Fortran compiler, CMake, BLAS, and LAPACK are required for this isolated
suite. GNU Debug adds bounds/runtime checks and invalid/zero/overflow traps.
The full repository test path remains `make test`.

At configure time, CMake extracts `integral_part` directly from the current
`ripple_solver_ArnoldiOrder2_test.f90`, writes a build-local include, and records
its SHA256. It does not rewrite production sources. The mock collective first
checks the independently executed native local projection, then supplies the
other species' projected data. All species are exercised. This prevents a
shared error in the new projection from being hidden by the gather adapter.

The tests cover 42 native-kernel comparisons with complex inputs, odd/even
Legendre parity, three species, ragged pitch grids including zero passing bins,
both bounce settings, endpoint impulses and zero sources. Negative controls
prove that omitted bounce and cross-species terms are detected. An end-to-end
manufactured problem uses this native collision kernel and a diagonal P,
checks the raw native equation, and reconstructs all 270 field coordinates from
180 moment coordinates. That P is deliberately a surrogate, not a captured
NEO-2 streaming matrix.

Independent dense LAPACK tests cover full/reduced solution and spectrum,
singular collision factors, unstable Richardson eigenvalues, actual GMRES
restarts, unit feedback/nullspace modes, incompatible forcing, exact borders,
regularizer retention, invalid inputs, tiny nonzero sources, nonfinite or failed
callbacks, exhausted budgets, and inaccurate inverses.

See [the validation record](ql_moment_reduction_validation.md) for measured
errors, compiler flags and the explicit untested scope.

## 6. Integration gates and coordination

The new branch is independent of the live solver branches; it does not edit
`arnoldi_mod` or `ripple_solver_ArnoldiOrder2_test`.

- #213: the complete Arnoldi projection remains a separate correction.
- #223: retain original-equation acceptance checks during integration.
- #247: compare at equal *observable* accuracy, not just GMRES tolerance.
- #248: collision-model-dependent nullspaces must be retained explicitly.
- #254: determine correct physical constraint combinations and report all
  multipliers; do not assume a bordered solution solves the unmodified equation.
- #49: this implementation still uses P^-1; it does not implement the proposed
  alternative full-operator solver without inner solves and does not close #49.

Before enabling a production selection:

1. Connect the real retained sparse inverse and raw M action, extract the actual
   R=P-M support, and verify the augmented factorization on every admitted branch.
2. Implement the species-distributed vector layout, collective inner products,
   callback failure propagation, and residual norms; run multi-rank tests.
3. Compare all source/species channels and transport observables on the same
   golden inputs without changing reference data or tolerances.
4. Run the reported n1e21 and lag-14 decks, including field reconstruction,
   conservation, original residual, setup time, solve time and memory.
5. Report whether any reduction beyond classical moment/Schur algebra exists.
   No Stafford-dependent gain has been established by these tests.

The mathematical reduction and the serial executable reference are complete in
their stated class. These production integration/qualification gates are not.
