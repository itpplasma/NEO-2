# Order2 integral-part solve convergence

The `settings` namelist accepts `epserr_iter`, `niter` and `n_arnoldi`.
Order2 defaults are `1e-7`, `100` and `500`; explicit values are honored.
The other ripple solver retains its historical `epserr_iter=1e-5` default.
Consumed controls appear in `neo2_config.h5` under `settings/epserr_iter`,
`settings/niter` and `settings/n_arnoldi`.

For every source solve, rank zero prints one `NEO2_SOLVE_JSON` record.
The record identifies `system`, `source` (1 through 3), and
`driving_species` (zero based), and gives `tolerance`, `iteration_limit`,
`arnoldi_limit`, `iterations`, `normalized_equation_residual`, `converged`
and `reason`. A failed solve prints its failure record and terminates on
all ranks before transport coefficients can be published.

The normalized equation residual is the maximum species normwise backward error:

    ||q + C_integral f - M f||_1 /
      (||q||_1 + ||C_integral f||_1 + ||M f||_1).

Each species has its own numerator and denominator; their ratios appear in
`normalized_equation_residual_by_species` and the scalar is their maximum.
`residual_norm` is `max_species_normwise_backward_error_l1`. Here `M` is the
assembled kinetic operator (`amat_symm` for the axisymmetric system,
`amat_sp` for the nonaxisymmetric system); `C_integral` is applied with the
same coupled collision action used by the equation. Evaluation occurs
before the sparse inverse, conservation-nullspace projection, and Arnoldi
stabilization. Thus the reported error is neither the preconditioned
Richardson update nor a proof of a transport discretization error bound.
The separate `density_constraint_residual` and `energy_constraint_residual`
diagnostics report normalized moment constraints and never replace the
equation residual.

## Axisymmetric conservation constraints

The axisymmetric iteration solves the bordered system

    [L  K] [f     ]   [q]
    [B' 0] [lambda] = [0],     L = M - C_integral,

where the columns of `K` are the species' Maxwellian density and (with the
integral part) energy kets and `B` holds the matching moment bras. Its
Richardson step with preconditioner `[P K; B' 0]`, `P` the factorized
regularized matrix, reduces to the plain step `y = f + P^-1 (q - L f)`
followed by

    f_new = y - G (B' G)^-1 B' y,     G = P^-1 K.

A fixed point satisfies `B' f = 0` and `L f + K lambda = q`. The kets are only
approximate null vectors of the coupled operator, so `K lambda` is the part of
the equation residual that the constraints cost; it is included in the
reported residual. The earlier filter subtracted half of the ket component,
`y - (1 - tol) K (B' y)/(B' K)` with `tol = 0.5`, whose fixed point carries
the residual `(1 - tol) P K c` instead.

Zero numerator and zero denominator give zero residual. Invalid controls,
nonfinite residuals, Arnoldi/subtraction failures, and iteration exhaustion
are failures. Existing iteration progress lines remain available.
