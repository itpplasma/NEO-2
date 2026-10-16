# Explicit Boozer field-provider lifecycle

`neo_field_provider` exposes `neo_field_load(path, input_switch)`,
`neo_field_prepare(s)` and `neo_field_free()` for library consumers that do not
own `neo.in`. Input switches 8 and 9 describe the existing .bc column formats.
The file-header minor/major radii are in metres, as in `neo_read`; the loader
publishes those values and the tabulated radial interval. Prepare rejects an
out-of-range or nonfinite surface. Free is safe before load and after free.
Loading a different file releases the old spectra and angular caches.

The implementation calls the existing `neo_init` and `magfie` routines, with
laboratory 10's geometry-derived Jacobian/current orientation. It does not
introduce alternative field equations. Angular interpolation uses 1025 poloidal
and 9 toroidal nodes by default, cubic radial interpolation and magnetic shear.
Optional `theta_points` and `phi_points` select other grids before cache
allocation; each requires at least 5 points. Those are
reference-provider settings qualified on the committed fixtures, not a general
accuracy bound for arbitrary spectra. No error enclosure is provided. Consumers
must qualify other inputs at their required observable accuracy.

The historic provider has one global instance. Load/prepare/evaluate/free must
be serialized; calls cannot interleave with an independent NEO-2 calculation.
This API does not initialize MPI or implement distributed field ownership.
`magfie` retains its historic conventions: `bmod` is in tesla, while the metric,
unit-vector components and exported hat quantities use the documented CGS
coordinate conventions. Setting `partpa_mod:bmod0=1e-4` makes the exported field
and current hats directly usable in gauss and gauss-centimetres, respectively.

The perturbation module exposes `neo_load_pert(path, input_switch)` and
`neo_free_pert()` independently of the background input format.
`neo_magfie_pert_amplitude(x, amplitude)` evaluates the existing Fourier field at
phi=0, divides by the existing `bmod0`, and clamps the radial point to its file
interval. It requires one toroidal harmonic. The existing full field remains
`neo_magfie_pert(x, field)` and retains the toroidal phase. Type 8 uses the
provider's retained negative-m phase; type 9 uses `(bmnc-i*bmns)*exp(i*m*theta)`.
Consumers that need a cosine envelope must explicitly take the real part for
type 8; that envelope is distinct from the provider's full complex coefficient.

`field_provider_conventions` checks both handedness signs and iota signs,
reloads different field strengths, returns to the original profile, and checks
signed-m perturbations, sine/cosine phases, normalization and repeated free.
The failure cases reject unsupported formats, out-of-range surfaces, angular
grids with fewer than five nodes and a multi-toroidal-harmonic amplitude request.

The CGS `magfie_result=0` curl now includes the angular derivatives of its
nonzero radial covariant field. The retained [FortSym definition and independent
oracles](complete-covariant-curl.md) specify its signed-Jacobian convention and
qualification limits. NEO-RT separately keeps its zero-radial-field
approximation and marks its toroidal curl component unavailable; this upstream
repair does not change that consumer choice. The unused `magfie_result=1`
legacy output and the missing radial derivative in `neo_magfie_b` remain outside
this qualification. No rigorous field certificate is supplied.

## Laboratory-10 curl orientation

In the stored coordinate order `(s, phi, theta)`, the contravariant curl is
`(d_phi h_theta-d_theta h_phi, d_theta h_s-d_s h_theta,
  d_s h_phi-d_phi h_s)/J`. All three numerators therefore change by the same
factor when the denominator's orientation changes. The historical reduced curl
expressions use `J0=-psi_pr*(curr_pol+iota*curr_tor)/B**2*1e6`. Laboratory 10
replaces that sign by the geometry sign: `J=abs(J0)*sign(J_geometry)`.
The previous reduced formulas required a matching orientation factor
`sign(J0)*sign(J_geometry)`. The complete generated curl instead divides by the
final signed J directly, so that additional factor is no longer applied.
The laboratory-specific Jacobian/current convention itself is unchanged.

The independent regression differentiates the returned covariant field at two
step sizes and Richardson-extrapolates. It checks the determined axisymmetric
radial and poloidal components at five angles, both chart orientations and both
iota signs, with the existing `1e-8` relative gate. A laboratory-8 case checks the
unchanged Jacobian/current orientation path. The additional complete-curl test
checks all three components against derivatives of the actual returned
covariant field on the separately documented manufactured 3D charts. This
checks interpolation consistency, not arbitrary-spectrum physical accuracy.
