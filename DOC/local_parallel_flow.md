# Local parallel flow of the standard QL drives

The default is off. Setting `lsw_local_parflow=.true.` in `&ntv_input` writes
`local_parallel_flow.h5` next to `neo2_multispecies_out.h5`. It requires
`isw_ripple_solver=3` and `lsw_multispecies=.true.`. All other outputs are
unchanged.

## Where the average is taken

`solve_eqs` in `ripple_solver_ArnoldiOrder2_test.f90` forms
`qflux = 0.5*Re(conjg(flux_vector) . f)`. Row 2 of `flux_vector` holds the RK
step weight times the parallel-flow moment `flux_vector_plot(2,:)`
(`convol_curr` times `weightlag(2,m)`, with `sign_of_bphi`) at each step of the
field line. The axisymmetric (AX) rows are then doubled. `compute_Dijab`
divides by `y(6)` and converts `qflux(2,:)` into `D31`, `D32` and `D33`.
`compute_ParFlow` contracts them with the forces into `ParFlow_*_spec`, which
is `<V_par B>`. The new output keeps the per-step moments before the sum over
steps.

## Grid and fields

Arrays are written in Fortran order. h5py shows them reversed.

- `phi_mfl`, `theta`: Boozer angles of the RK steps on the field line.
  `theta = boozer_theta_beg + iota*(phi_mfl - boozer_phi_beg)` is not reduced
  modulo 2*pi. The line covers one poloidal turn.
- `bhat` is B/Bref and `Bref` is in G.
- `surface_weight` is `dphi/(bhat*h^phi)` with the RK weights.
  `surface_volume` is `y(6)`, the normalisation of D3j and ParFlow.
  The flux-surface average used by NEO-2 is
  `<x> = sum(surface_weight*x)/surface_volume`.
- `D3j_AX_real/imag`, `D3j_NA_real/imag` have the shape
  (step, j, drive species b, flux species a), where j=1,2,3 multiplies A1, A2
  and A3. Their units are those of `D3j_*` in `neo2_multispecies_out.h5`.
- When `isw_calc_Er` is nonzero, `A_forces` (A1, A2, A3; species) uses the
  formulas of `compute_ParFlow`. A3 is nonzero only on its multispecies
  `isw_calc_Er=1` branch. `VparB_*` (step, species) equals
  `-sum_b sum_j D3j*A_j^b`, the local V_par*B in cm*G/s, positive along the
  equilibrium field. `j_par_*` is `sum_a e*Z_a*n_a*VparB_a/(Bref*bhat)` in
  statA/cm^2, also positive along B.

## Averages

For the AX problem `Re <D3j_AX> = D31_AX`, `D32_AX` or `D33_AX`, and
`Re <VparB_AX> = ParFlow_AX_spec`. The imaginary AX part is roundoff.

The non-axisymmetric (NA) solution is linear in the perturbation and carries
the toroidal phase of `bnoverb0`. The physical field at a sampled step is
`Re[value]`, and elsewhere on the surface it is
`Re[value*exp(i*n*(phi - phi_mfl))]`, with `n` the perturbation's `m_phi`.
The existing NA rows are `(1/2) Re <.>` of a single field-line sum:
`D3j_NA = 0.5*Re <D3j_NA local>` and `ParFlow_NA_spec = 0.5*Re <VparB_NA>`.
This sum does not average over the toroidal phase. The full flux-surface average of
the n/=0 helical flow vanishes.

`ctest -R local_parallel_flow` checks both identities to roundoff on the
synthetic `phi_ma` fixture with an added n=1 |B| perturbation. It also checks
that switching the output on leaves every legacy dataset bitwise unchanged.
Set `NEO2_LOCAL_PARFLOW_BASELINE` to a binary without this change to compare
the switched-off run with it bitwise.
