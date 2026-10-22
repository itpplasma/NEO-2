# Prescribed misaligned potential in the local QL solver

The default is off. In `&ntv_input`, `phi_m` is a complex potential in
statvolt, `delta_m` is a complex surface-label corrugation in cm, and
`phi_ma_m` is the signed poloidal mode number. Both complex amplitudes default
to zero. The toroidal mode is the perturbation file's existing `m_phi`,
including its field-period multiplier. The real drive is
`Re[(phi_m + Er*delta_m)*exp(i*(phi_ma_m*theta + m_phi*phi))]`.
For a pure potential drive, set `delta_m=(0,0)` and supply a zero-|B|
perturbation file with the desired toroidal mode. Delta is a label corrugation,
not the geometric displacement of the surface, whose sign is opposite.
Delta must belong to the same perturbed equilibrium and radial label, and
cannot be inferred from a scalar |B| perturbation.

`Phi_0'=-Er` follows from the force in `ntv_mod.f90:3725`–3727. The solver
recovers Er from its native rotation using the inverse relation at
`ntv_mod.f90:1580`, whose forward definition is at `ntv_mod.f90:2379`.
The signed Jacobian follows the same geometry conversion as that routine,
including the left-handed Boozer branch. A zero misaligned amplitude adds no
quadrature or solve. A configured aligned pair writes a zero extra current.

The supported path is the local tokamak Arnoldi O2 NA solve:
`magnetic_device=0`, `mag_coordinates=1`, `isw_axisymm=1`,
`isw_ripple_solver=3`, `isw_qflux_NA=1`, `prop_reconstruct=0`, and explicit
species profiles (`lsw_multispecies=T`, one MPI rank per species). The new
moment requires an initialized nonrelativistic collision basis. One additional
collective RHS contracts the A1 and A2 forces before solving and shares the
existing sparse factorization. The legacy force columns are preserved.

`phi_ma_current.h5` follows the existing HDF5 output style. `j_par_real` and
`j_par_imag` are the total extra local current in statA/cm², positive along
the equilibrium magnetic field, sampled at `theta` and `phi_mfl`. They retain
the toroidal phase. `surface_weight` is the existing RK volume quadrature.
The ordinary line average and Bhat-weighted line average are separate. The
latter closes to `current_row_real/imag`; the transport current row measures
⟨j_parallel B⟩/Bref. No Fourier half-product enters a complex linear current.
The full-surface average `j_par_fsa_real/imag` is exactly zero for n≠0 on this
axisymmetric background. A single field-line average need not vanish.

`normalize_phi_ma_j_hat` is the normalization boundary specified by
neo2-phi1 `derivation/L1/L1.md`. For one species with A2=0 and nonzero A1,
Er and Phi_MA, the file also contains sampled `J_hat_real/imag`, using signed
charge and `v_T,K=sqrt(2*T/m)/sqrt(2)`. SI conversion constants are stored
beside the data. This normalization alone does not reconcile NEO-2's native
rotation frame with the laboratory ExB frequency used by a KILCA comparison.

Run the focused tests with `ctest -R phi_ma --output-on-failure`. Set
`NEO2_PHI_MA_BASELINE` to an origin/main executable to enable the bitwise
legacy-dataset comparison; that test reports skipped when the binary is absent.
The small committed fixture is synthetic and checks linear behavior, not
collisionless susceptibility or resolution convergence.
