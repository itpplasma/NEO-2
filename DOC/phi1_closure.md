# Galerkin Phi_1 closure (L3)

**Hypothesis — discrete contract.** `isw_phi1=2` extends the L2 basis-response
mode. The master assembles `M=sum(e_a R_a)-chi G`,
`rhs=-sum(e_a density_dB_a)` and solves the complex system with pivoted LU.
The coefficient vector and diagnostics are broadcast to the kinetic ranks.
The backward residual must be <=1e-8 and reciprocal condition >=1e-12.
`isw_phi1=1` remains the independent density-response mode, and zero is off.
The route restrictions and units are those in [phi1_response.md](phi1_response.md).

**Hypothesis — observables.** The enabled `phi1` HDF5 group adds complex
coefficients, matrix/RHS, sampled physical density, residual/conditioning,
magnetic/ExB/total particle and energy flux, and species/total torque. Fluxes
use `ftilde` with the real single-n quadratic factor 1/2. The ExB contraction
uses the sampled density/energy moments and exactly the RK weights used by QN.
`qn_flux_check=sum(e_a Gamma_ExB_a)-chi_Phi1_Ur` is the D9 physical-density
charge-flux check. The Boltzmann subtraction is retained explicitly.

**Hypothesis — physics boundary.** D9 in
[neo2-phi1](https://github.com/itpplasma/neo2-phi1) does not establish the
additional toroidal ordering that drops second-order magnetic-drift terms.
The reported torque is the requested `ftilde*v_g1` product with the existing
torque conversion. Its interpretation as complete physical NTV remains open,
as do kinetic/basis convergence and comparison with SFINCS `includePhi1`.
The unequal-temperature source follows the retained D10 collision action.

**Hypothesis — exchange contract.** The master writes `phi1.bc`, one surface
in the ten-column Boozer perturbation layout. The mandatory quantity comment
identifies statvolt, and `bmnc=Re(Phi_j)`, `bmns=-Im(Phi_j)`. Geometry columns
are zero. Signed m rows, n and nper follow the magnetic perturbation, while
the global geometry header and interpolated surface metadata come from the
Boozer background. This first writer requires inp_swi=9 and lab_swi=9 or 10.
The single surface is not a radial profile. NEO-RT's constant-in-s handling
and its physical comparison are the next layer.

**Hypothesis — checks.** `phi1_qn_test` uses the complex nonidentity-Gram
fixture with provenance, zero forcing and singular/ill-conditioned rejection.
`phi1_closure_test` checks a driven, zero and unequal-T synthetic case, the
independent unscaled matrix residual, D9 and the production `neo_read_pert`
round trip, including the coefficient sign and surface header.
