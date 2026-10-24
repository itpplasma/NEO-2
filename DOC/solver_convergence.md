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
Density and energy projection remains in the iteration. The separate
`density_constraint_residual` and `energy_constraint_residual` diagnostics
report normalized moment constraints and never replace the equation residual.
For coupled species, per-species energy projection need not be a true
nullspace gauge; an exposed residual floor requires investigation.

Zero numerator and zero denominator give zero residual. Invalid controls,
nonfinite residuals, Arnoldi/subtraction failures, and iteration exhaustion
are failures. Existing iteration progress lines remain available.

`NEO2_AX_ENERGY_FILTER_POLICY=full` is an opt-in numerical iteration strategy for
resolved `full_fp`. It removes the full centered energy component using the native
bra, ket and denominator, while retaining the historical density damping.
Unset or `native` preserves the historical arithmetic. Invalid values, other
collision models and combination with `NEO2_AX_FILTER_POLICY=none` are rejected.
The physical operator, source, NA arithmetic and original residual plus enabled
moment acceptance are unchanged. This changes the iteration map and requires a
fresh Arnoldi construction. It is not a normalization repair or a claim that the
removed direction is a physical right null vector.

`NEO2_AX_ONLY_FEASIBILITY=1` is a default-off, diagnostic-only boundary for the
three-species full-energy strategy with original acceptance and unchanged
1e-5/100/500 controls. It rejects conflicting probes and invalid invocations.
Only after all nine distinct AX source/driver solves have passed the existing
native equation and enabled-moment validator does it stop before NA, with
`AX_STRICT_FEASIBILITY NO_TRANSPORT`. Missing, duplicate or failed solves cannot
produce this marker. Unset or `0` retains normal execution. The marker attests
AX feasibility only and never admits transport.
