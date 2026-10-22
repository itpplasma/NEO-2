# Phi_1 basis density responses (L2)

**Hypothesis — scope.** Opt in with `isw_phi1=1`, `phi1_num_m=K`, and
`phi1_m=m1,...,mK` in `&ntv_input`. Each unit-statvolt harmonic uses the
single nonzero toroidal mode of the magnetic perturbation, with phase
`exp(i*(m*theta+n*phi))`. This layer returns responses, without QN closure.
The direct nonrelativistic Boozer tokamak O2 route, explicit kinetic species,
a neutral equilibrium, one MPI rank per species and fixed Om_tE are required.
Use either this basis mode or L1's prescribed Phi^MA input in one run.

**Hypothesis — normalization.** D7–D10 in
[neo2-phi1](https://github.com/itpplasma/neo2-phi1) define the source and
moments. L1's x^-1 projection and RK source stencil are reused. Separate
A1/A2 responses are contracted with the physical driving beta and forces
once. The independent density-shift collision source is retained for unequal
T and vanishes at common T. These are discrete responses, not physical
convergence or benchmark acceptance.

**Hypothesis — output contract.** The enabled-only `phi1` HDF5 group contains
`m`, `n`, phase origins, species tags, RK surface weights, the Gram matrix,
`density_dB_re/im` [cm^-3] and `density_response_re/im` [cm^-3/statvolt].
The response has Fortran axes (test harmonic, drive harmonic, species),
reversed by HDF5. The density bra includes both streaming branches and
turning cells, without a Fourier half factor. `collision_source_norm`
diagnoses D10. The default-off path keeps L1 outputs unchanged.

**Hypothesis — checks.** `phi1_density_test` uses the independently probed
D8 fixture, whose header gives provenance. `phi1_response_test` reuses L1's
synthetic circular geometry, reflects its Z sine coefficient to use
right-handed angles, and supplies two synthetic species. It compares
the disabled executable with `NEO2_PHI1_BASELINE` (an L1 executable), tests
basis independence, the aligned m=0 and zero magnetic-drive nulls, and D10.
