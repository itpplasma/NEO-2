# Native work-only response diagnostic

Set both `lsw_response_witness=.true.` and
`lsw_work_only_response=.true.` in `&ntv_input`. Both default to false.
The additional option performs one combined physical work-driver solve using
the nonaxisymmetric factorization and Arnoldi cache before their normal cleanup.
It retains the actual MPI-coupled integral collision operator. Original source
channels, solved responses and legacy qflux outputs are preserved.

This diagnostic currently admits multispecies Arnoldi O2 NTV, prescribed
`isw_calc_Er=2`, an axisymmetric background (`isw_axisymm=1,n0b=0`),
Boozer lab 8/9 with the positive native chart, no magnetic shear, flat density
and temperature profiles, nonzero prescribed rotation and nonzero toroidal
perturbation mode. All ranks must use the same input deck, as required by the
native receiving-species MPI solver. Every receiving rank participates even
when its local work RHS vanishes.

The complete complex relative perturbation must be constant on every consumed
field-line node, including boundary extensions. Its variation must be at most
`128*epsilon*abs(b0)`, and its derivative times `max(1,abs(phi_end-phi_begin))`
must satisfy the same bound. Its complex phase must agree between ranks.
Periodic seam closure and the actual source-geometry split are checked.
The axisymmetric physical response's imaginary L2 norm must be at most
`128*epsilon` times its real L2 norm. These are scope checks, not convergence
targets. Unsupported controls or failed checks stop the opt-in run.

Physical source weights are obtained from the existing
`compute_electromagnetic_work` API with the exact prescribed-Er reporting
operation order. All driving species, including off-diagonal responses, enter
the reconstruction. Actual A3/A2 weights must vanish. The admitted radial
source split is `Bwork=Bforward+b0*Baxis`; channel 2/Ware is excluded.

The new datasets are appended to each uniquely named raw witness file.
Complex vectors use separate `_real` and `_imag` datasets, following the
original native equation-row ordering, including constraint and seam rows.

| Dataset prefix `work_only_` | Meaning |
|---|---|
| `Bwork,Bcomplete,Y1physical,M,Ywork` | Physical combined work RHS, complete RHS, original response, `b0*real(Y0physical)`, additional response |
| `LnYwork,residual` | True unpreconditioned `L_N*Ywork` and `L_N*Ywork-Bwork` |
| `LnM,RN,Z,LnZ` | `L_N*M`, `Bcomplete+L_N*M-Bwork`, `Y1physical+M-Ywork`, `L_N*Z` |
| `identity_residual,forward_residual` | `L_N*Z-RN`, `L_N*Y1physical-Bcomplete` |
| `global_l2_norms` | Across all receiving MPI ranks: work residual, Bwork, identity residual, original forward residual |
| `receiving_raw_moments` | Native `Re(bra^H*value)/2` for Y1, Ywork, M, Z, then `W(Z)-W(Y1)+W(Ywork)` |
| `shift_work_scale,shift_work_gate_pass` | Absolute cancellation scale and roundoff check for actual native shift work |
| `raw_to_torque` | Multiply each raw receiving moment by `-2*n_a*T_a/(scale*dl_over_b)` |
| `cpu_seconds,wall_seconds,cost_scope` | Per-rank measured cost of the additional solve, four operator applications and vector I/O; excludes preparation/guard overhead |

Source weights, actual early Er, geometry-split defect, perturbation variation,
derivative, axisymmetric imaginary contamination and seam phase are also saved.
The true operator is `L_N*x=A_n*x-integral_part(x)`, including all native rows.
It does not substitute the regularized axisymmetric preconditioner or a sparse
backsolve for that operator. Global residual norms include every receiving rank.

The extra payload is approximately `208*N` bytes per receiving rank for
`N=n_2d_size`, plus scalar metadata: 10.4 MB at N=50000. There are eight
additional complex state vectors live at the solve hook, about `128*N` bytes,
plus geometry, MPI/solver temporaries and HDF buffers. These are size estimates;
measure actual files, peak RSS and wall time in native cluster qualification.
Output streams through the existing witness handle; no matrix is duplicated.

`work_only_diagnostic_complete=1` means the diagnostic returned. Both
`work_only_native_solve_qualified` and `work_only_nr_continuum_identification`
remain zero. Original collision-action and residual-response qualification
flags also remain zero. The saved `RN` is an assembled discrete-coordinate
remainder. It has no established continuum NEO-RT identification, and raw
collision assembly is not relabeled a physical `C*f10`.

The owning `work_only_response_test` uses an analytic complex 3x3 coupled
system with a seam entry, an independent elimination solve, four-phase work
invariance, a deliberately inaccurate solve and unsupported-control checks.
It verifies the receiving-rank sum callback contract, not MPI execution.
Native qualification must additionally establish source/binary closure,
parent/default/off and on/off legacy-output parity, actual source-weight
reporting parity, all-rank residuals and phase/shift checks, work accuracy,
refinement and measured cost under the project cluster policy. No full native
compile, actual coupled solve or collision-model equivalence follows from the
helper test.
