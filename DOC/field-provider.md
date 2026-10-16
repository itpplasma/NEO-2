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
and 9 toroidal nodes, cubic radial interpolation and magnetic shear. Those are
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
The failure cases reject unsupported formats, out-of-range surfaces and a
multi-toroidal-harmonic amplitude request.

The original general `magfie` curl implementation omits derivatives of its
nonzero radial covariant field from two components. This lifecycle API does not
claim to repair that scientific limitation. NEO-RT explicitly uses its existing
zero-radial-field approximation and marks the toroidal curl component unavailable;
its separately retained curl identity tests qualify the determined axisymmetric
components. No complete general-field curl or rigorous certificate is claimed.
