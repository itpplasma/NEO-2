# Native response witness for the multispecies Arnoldi solver

Set `lsw_response_witness=.true.` in `&ntv_input` to preserve native response
coefficients and source assembly before their storage is reused. This diagnostic
requires `lsw_multispecies=.true.`, `isw_ripple_solver=3`, `isw_qflux_NA=1` and
`isw_calc_Er>=1`. The default is false and creates no witness files.

Each receiving MPI rank writes one
`response_witness_p<TAG>_a<SPECIES>_c<INVOCATION>.h5` per solver call.
Species indices start at zero. The invocation counter prevents an overwrite
when the same propagator is solved again. An existing filename is rejected.
HDF streaming applies only to the new witness handle; the previous global
creation setting is restored before the solver resumes. Native outputs,
`collop.h5`, sources and equations are unchanged.

## Saved coefficients

Complex arrays are represented by separate `_real` and `_imag` datasets.
Their stored axes are described in Fortran order below; HDF readers may present
the axes in reverse order. The first axis uses native equation-row numbering.
The last axis of the response/RHS coefficient arrays denotes the driving
species, including coupled responses driven by a different species than the
receiving rank. Its entries follow `driving_species_index=0:Ns-1`; the complex
writer stores ordinary one-based array bounds. Complex field-line arrays follow
the same positional order as the real geometry arrays, beginning at `ibeg`.

| Dataset | Fortran axes and meaning |
|---|---|
| `axisymmetric_rhs` | row, channel, driving species; RHS before the first solve |
| `axisymmetric_y0` | row, channel, driving species; complex solved axisymmetric coefficients |
| `axisymmetric_auxiliary4` | row, driving species; unsolved incompressibility source |
| `collision_differential_b<B>` | row, channel; native sparse differential action on `real(Y0)` |
| `collision_integral_b<B>_k<K>` | row; native complex integral action before real-part subtraction |
| `collision_combined_native` | row, channel, driving species; actual combined native array |
| `mirror_native` | row, channel, driving species; actual native pitch-derivative product |
| `forward_rhs_b<B>` | row, channel; forward RHS before the axisymmetric-response addition |
| `complete_rhs_b<B>` | row, channel; consumed RHS immediately after that addition |
| `nonaxisymmetric_y1` | row, channel, driving species; solved complex nonaxisymmetric coefficients |
| `work_flux_vector` | receiving bra channel, row; actual final native work bra |

Solved channels are `(1,2,3)=(A1,A3,A2)`. Channel 4 is never labeled a
solved distribution. The witness preserves the imaginary part of `Y0` even
though the native collision/mirror construction uses its real part.

The file retains `ind_start,npl,eta,alambd`, the consumed field-line geometry,
complex perturbation and derivative, signed geometry scales, energy source and
bra moments, all four RK factors, species profiles and consumed rotation
parameters. Only defined endpoint ranges are written for `delt_pos`
(`ibeg+1:iend`) and `delt_neg` (`ibeg:iend-1`); their dataset bounds retain these
indices. Native bin widths include each position-dependent final partial bin.
For spatial point `j` and energy index `m=0:lag`, its row block begins at
`ind_start(j)+2*(npl(j)+1)*m`. Positive pitch has full bins then the final
partial bin; negative pitch reverses this order. These coefficients are bin
integrals, not pointwise distribution samples.

Keep the same run's `collop.h5`, source revision, input deck and
`neo2_multispecies_out.h5`. The latter supplies the actual `Er`, inductive
drive and `TphiEM_source_weight`; reconstruct with driving-species weights.
A pure electric drive requires its actual `(A3,A2)` weights to vanish for all
species. `collop.h5` supplies the basis definition and pair operators; its
999-point velocity table is a sampling grid, not the energy quadrature.

For `N=n_2d_size` and `Ns=num_species`, the main coefficient payload with
integral collisions is approximately `376*N*Ns+48*N` bytes per receiving rank,
plus geometry, metadata and HDF overhead. D/T/e at `N=50000` is approximately
59 MB per rank. This estimate is not a measured pilot size. Output is linear
in the state size; no equation matrix is copied. Budget actual bytes, memory
and write time on the cluster before a larger diagnostic campaign.

## Qualification boundary

`complete=1` is written only after the native nonaxisymmetric solve returns
and all requested arrays are saved. It certifies a complete capture, not
solver convergence. `physical_collision_action_qualified` and
`residual_response_qualified` remain zero.

In particular, the differential collision action has local RK factors removed,
while the integral action includes neighboring-node and endpoint weights.
`collision_combined_native` must not be relabeled a pointwise physical `C*f10`
without an independently checked reconstruction. The saved separate actions,
weights and complete RHS permit that check without another lossy projection.
For relative perturbation constant along the field line, the complete-source
identity can then be tested under its stated coordinate assumptions. This
observer does not perform a residual-driver or adjoint solve and does not
establish equivalence to a NEO-RT response at finite collisionality.

The focused I/O tests verify exact complex reconstruction, untouched inputs,
default/disabled behavior, writer-setting restoration and overwrite rejection.
Native qualification additionally requires the project's full registered checks,
parent/default/off output parity, on/off numerical parity, all ranks/propagators'
complete files, actual input/drive readback and independent reconstruction of
the source and work contractions. A later same-operator residual response and
numerical projection/convergence bounds are separate gates.
