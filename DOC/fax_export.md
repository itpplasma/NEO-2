# Axisymmetric distribution export

Set `lsw_export_fax = .true.` in `ntv_input` with the order-2 Arnoldi
ripple solver (`isw_ripple_solver=3`). The default is off. An axisymmetric
solve writes `fax_background_spec<RANK>.h5` before its solution is reused
or discarded. Use one MPI rank per species. A background-only calculation
can set `isw_qflux_NA=0`; the export also works when the NTV solve follows.
It does not change either equation or any computed coefficient.

The files are a set: **every field-species file is required**. The MPI rank
is the zero-based field species; `species_tag` maps this order to profile tags.
The three source columns use NEO-2's internal order: density, parallel electric
field, temperature. `source_vector_all`'s final axis enumerates the **driving**
species, including off-diagonal collision responses. No axis is summed away.

## HDF5 version 1

`format_version=1`. Array dimensions below are in Fortran order. HDF5/Python
readers see reversed dimensions; the supplied Python reader reverses them once.
HDF5 bounds attributes retain the Fortran index origins.

- `g_real`, `g_imag`: `(solution_index, source=1:3, driving_species=0:N-1)`.
  These are the native axisymmetric solution's **pitch-band integrals**, not
  point values. Each position contains blocks `m=0:lag`; each block contains
  `2*(npassing+1)` entries. Positive parallel velocity occupies the first
  half in increasing band order; negative velocity occupies the second half
  in reversed band order. The incomplete last band ends at `eta=1/bhat`.
- `offset`, `npassing`, `ibeg`, `eta`: the native index and pitch support.
  At position `t`, positive-band `i=1:npassing(t)+1`, block `m`, the Fortran
  solution index is `offset(t)+2*(npassing(t)+1)*m+i`. The corresponding
  negative-band index is `offset(t)+2*(npassing(t)+1)*(m+1)-i+1`.
- `phi`, `bhat`, `dlogb_dphi`, `h_contravariant_phi_abs`: field-line position/branch geometry;
  `bhat=B/B_ref`. `phi_direction=sign(b^phi)` is saved before the native
  solver replaces the contravariant field component by its absolute value.
  Native sigma labels increasing/decreasing phi; physical signed pitch is
  `lambda=phi_direction*lambda_native`. The position-chart branch is likewise
  `sign(phi_direction*dlogb_dphi)`. These are not labels for multiple wells.
- `basis_coeff`: `(power=0:degree, block=0:lag, segment)`. At
  `x=v/v_T=sqrt(E/T_field)`, evaluate the polynomial in
  `x-basis_centers(segment)`, using `basis_edges` to choose a segment.
  `basis_extend` selects polynomial continuation beyond the final edge;
  otherwise the basis is zero there. Supported bases are the native
  normalized associated Laguerre basis (0), quadratic B-splines (11, order 3),
  and the monoenergetic Lorentz block (unity). Unsupported bases and
  relativistic models fail explicitly when export is requested.
- `collop_base_exp`, `collop_base_prj`, `bspline_order`, `scalprod_alpha`,
  `scalprod_beta`, `lag`, `isw_lorentz`: basis and operator metadata. The
  historical variable `lag` is a highest **basis index**; basis 11 is not
  Laguerre. B-spline coefficients include the actual knot spacing and tail.
- `mass_g`, `charge_number`, `density_cm3`, `temperature_erg`, `species_tag`:
  all species, including field and driving species. Gradients, when allocated,
  are `dn_ds_cm3` and `dT_ds_erg`, with respect to `s_tor`.
- `B_ref_G`, `iota`, `psi_pr_hat`, `scalefac_kG`, `bcovar_theta_hat`,
  `bcovar_phi_hat`, `y_ntv`, `Om_tE`, `isw_calc_Er`: native normalization.
  In particular `<|grad s|> = y_ntv(7)/y_ntv(6)` in Fortran indexing.
- `flux_real`, `flux_imag`, `qflux`: the native axisymmetric flux functional
  and its public (factor-two corrected) coefficient matrix for this field
  species, `(flux, source, driving_species)`. `flux_geometry` removes the
  energy weight and repeats the native position/pitch quadrature in every
  block. This allows an independent energy integral to recover `qflux`.
- `sample_x`, `sample_basis`, `sample_step`, `sample_band`, `sample_sign`,
  `sample_row`, `sample_g`: native `phi_exp` evaluations and eight distribution
  band moments at six energies. `sample_g` has axes `(source, driving_species,
  energy_sample, position_sample)`. Divide by the corresponding band width
  to obtain a band mean. Samples are before profile-force multiplication.

## Physical normalization

The definitions are those in `ripple_solver_normalizations_and_output.tex`,
Eqs. `multispec:fa1norm` and `multispec:fa1basisfun`:

`f_AX^a/f_M^a = sum_b,k,m beta_k^b A_k^b g_m,kb^a phi_m(v/v_Ta)`.

Here `beta_1=beta_2=rho_b/<|grad s|>`, `beta_3=B_ref`,
`rho_b=m_b v_Tb c/(z_b e B_ref)` (signed charge), and the physical source
order is density, temperature, parallel electric field. The file uses internal
order **1,3,2**. `A1=<|grad s|> n'/n - 1.5 A2 - z e Er/T`,
`A2=<|grad s|> T'/T`, `A3=z e <E_parallel B>/(T <B^2>)`.
All quantities are CGS. An electric field must be supplied explicitly or read
from the completed native multispecies output; the export occurs before the
final force-balance output and its `Om_tE` is not a substitute for that field.

## Reader, oracles and consumer table

With numpy, scipy and h5py installed:

```sh
python python/src/neo2_ql/fax.py fax_background_spec*.h5 --out oracle.json
python python/src/neo2_ql/fax.py fax_background_spec*.h5 --out table_oracle.json \
  --table fax_energy.dat --neo2-output neo2_multispecies_out.h5
```

The oracle compares native `phi_exp` and distribution samples at tolerance
`2e-10*max(1,abs(reference))`. It independently integrates the Maxwellian energy
weights (powers 4,3,6 in internal source order), then contracts the complete
pitch/spatial distribution to recover the native qflux matrix, with a relative
matrix-scale tolerance `2e-8` per driving species, also checked entrywise
with a zero floor of `1e-12` of that matrix scale. The monoenergetic Lorentz
model has no Maxwellian energy integral. These are reconstruction checks;
they do not establish solver or physical transport convergence.

For native CTest registration, configure `NEO2_FAX_ORACLE_FILES` with a
semicolon-separated list of native exports and `Python3_EXECUTABLE` with the
Python interpreter containing these packages. `fax_basis_test` is always
registered and compares the exported Laguerre/B-spline polynomials, including
the tail, with the original basis routines.

`--table` writes GORILLA_APPLETS `FAX_ENERGY_V1` (documented by its consumer).
It retains every block, source, driving species and field species, and
interpolates band means and their pitch derivatives onto a common position/
pitch grid. This single-well position chart rejects multiple branch reversals;
it does not qualify radial coverage, interpolation convergence, collision
operators or physical device transport. Supply `--consumer` with the native
APPLETS energy-test executable to compare its reconstructed distribution with
the native band moments (matrix-scale tolerance 2%, including interpolation).
