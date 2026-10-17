# Collision models in NEO-2-QL

`&collision collision_model` selects an operator. The default, `'native'`,
preserves the existing collision switches and their numerical implementation.
The named Coulomb modes resolve these switches as follows:

| `collision_model` | `isw_lorentz` | `isw_energy` | `isw_integral` | Operator |
|---|---:|---:|---:|---|
| `native` | unchanged | unchanged | unchanged | Existing input contract |
| `pitch_angle` | 0 | 0 | 0 | Pitch scattering at the requested energy resolution |
| `test_particle` | 0 | 1 | 0 | Coulomb pitch and energy scattering |
| `full_fp` | 0 | 1 | 1 | Linearized multispecies Coulomb operator with field-particle terms |
| `monoenergetic_lorentz` | 1 | 0 | 0 | Existing monoenergetic Lorentz operator |
| `krook_constant` | 0 | 0 | 0 | Species-diagonal constant physical damping |
| `krook_inverse_cube` | 0 | 0 | 0 | Species-diagonal damping proportional to `(v/vth)^-3` |

The monoenergetic Lorentz mode resets the energy expansion to its native
single-speed representation; it is distinct from `pitch_angle`. The native
monoenergetic mode does not support a multispecies solve. Named Coulomb modes
do not change collision coefficients, discretizations, or the existing species
pair weighting. `full_fp` couples the perturbed species distributions;
test-particle coefficients depend on the background mixture.

## Physical Krook rates and species

Set `krook_nu(1:num_spec)` in `&collision`, ordered like
`species_tag_vec` in `&multi_spec`. Every active entry must be finite and
strictly positive, in **s^-1**. These are damping frequencies, with no factor
of `2*pi` or toroidal harmonic number. For example:

```fortran
&collision
  collision_model = 'krook_constant'
  krook_nu = 100.0d0, 80.0d0, 4000.0d0
  isw_lorentz = 0
  isw_energy = 0
  isw_integral = 0
  isw_momentum = 0
  isw_relativistic = 0
/
```

For a D/T/e mixture, supply three native species entries with charges
`(+1,+1,-1)`, masses in g, densities in cm^-3, temperatures in erg, and the
same physical electric rotation. The rate order follows those entries; the
implementation has no D/T/e-specific species limit beyond the native input
capacity. The thermal speed is `vth_a=sqrt(2*T_a/m_a)` in cm/s. At fixed
normalized speed `u=v/vth_a`, the operators are

\[
C_a[g_a]=-\nu_a G(\eta)g_a,\qquad
C_a[g_a]=-\nu_a G(\eta)u^{-3}g_a.
\]

`krook_geometry_trapped` and `krook_geometry_passing` in `&collision` are
finite, strictly positive, dimensionless factors, both defaulting to 1.
`G(eta)` is the trapped factor for `eta > 1/B_max`, and the passing factor
otherwise. Here `B_max = fieldline%b_abs_max`, the same global maximum used
by the native pitch grid. To match PENTRC's CONSTANT model, set
`collision_model='krook_constant'`, `krook_nu=nu_a` in s^-1,
`krook_geometry_trapped=1/(2*eps)` (supply its numerical value), and
`krook_geometry_passing=1`. The factors also apply to `krook_inverse_cube`.

The pitch unknowns are band-integrated fluxes with a local cubic Taylor
reconstruction in `lambda=abs(v_parallel)/v`, not Legendre coefficients.
The geometry step is integrated exactly in each band's weak mass, splitting
at `lambda_t=sqrt(max(0,1-bhat/B_max))`. For a band `[l,u]`, with
`t=clip(lambda_t,l,u)-l`, its Taylor moment of order `n-1` is
`[G_p*(u-l)^n+(G_t-G_p)*t^n]/n!`, `n=1,...,4`; reconstruction follows this
integration. A crossing band therefore mixes the reconstructed unknowns.
The Arnoldi solver keeps a separate collision mass so rotation and magnetic
drift continue to use the unweighted pitch mass. Unit factors retain the
original arithmetic bit for bit.

Krook damping is nonconserving and has no response coupling between species.
The supplied physical rates can include any declared background collision
channels. They do not turn this scalar relaxation operator into the full
Fokker-Planck operator. The full FP mode retains native multispecies coupling.

Krook requires the nonrelativistic energy-basis representation, and explicitly
rejects Lorentz-only, energy-diffusion, integral, NBI, or precomputed FP-cache
switches. No conservation correction is silently applied. The zero-rate limit
is not a supported Krook solve: the steady axisymmetric operator would regain
nullspaces. Existing native models retain their original input behavior.

`conl_over_mfp_vec` still supplies the native positive inverse-length
normalizers `kappa_a` (negative inputs specify positive `kappa_a` directly).
These cancel from the physical Krook reaction. They remain involved in the
native pitch-grid construction, so pitch convergence is still required. An
FP comparison must independently match the native Coulomb normalization and
background mixture; equating `kappa*vth` to a selected scalar width changes
the collision model.

## Native weak reaction and normalization

The streaming energy mass uses `u^4` and the scalar reaction uses `u^3`.
For the native projection and expansion functions, define

\[
K^{(p)}_{mm'}=\pi^{-3/2}\int_0^\infty
u^{3-p+\alpha}e^{-(1+\beta)u^2}
\phi_m^{\rm prj}(u)\phi_{m'}^{\rm exp}(u)\,du,
\qquad p=0\text{ or }3.
\]

The existing `M_transform_inv` orthogonalization is applied. Constant width
has the same energy projection as native `x1mm`; inverse-cube width has a
finite `u^0` integrand for the standard `alpha=0` basis. No low-speed floor
is introduced. Both modes require integrable scalar-product weights.

The actual native pitch carrier acts on a constant reconstructed distribution
as `-kappa/(2*h_phi) * integral(deta/abs(lambda))`. Its factor of one half
follows from the consuming reconstruction coefficients. Consequently the
reaction matrix stored in `denmm_a` is

\[
\mathrm{denmm}_{K,a}=-\frac{2\nu_a}{\kappa_a v_{{\rm th},a}}K^{(p)}.
\]

The product gives positive damping on the kinetic left-hand side. `denmm_a`
is a scalar weak-mass carrier in these modes; no energy derivative is
computed. Coulomb pitch and field-particle matrices are zero. Both QL ripple
solver paths retain the native energy, pitch, and spatial discretizations.
Krook removes the density and energy nullspaces physically, so the native
periodic regularizer and Arnoldi density projection are disabled for Krook.

The resolved model is written to `neo2_config.h5` under `collision` and to
`collop.h5` under `meta`; active physical rates are `krook_nu_s_inv`.
Both groups also record `krook_geometry_trapped` and `krook_geometry_passing`.
`collop.h5/denmm_a` contains the energy reaction matrices; the geometry-weighted
pitch mass is assembled locally in each ripple solver.
The relaxed variable is labeled `native_thermodynamic_response`; its
Hamiltonian nonadiabatic identification is not yet qualified and
`collision_relaxed_variable_qualified` is zero. The native `TphiNA_spec`
observable is computed from the radial particle flux by the flux-force
relation in `ntv_mod::compute_TphiNA_a`. Equality to a direct electromagnetic
perturbation-work observable under nonconserving Krook damping remains an
additional qualification condition.

## Comparisons and verification

Constant and inverse-cube local damping can represent the corresponding
scalar MARS/PENTRC/NEO-RT rates. A quantitative agreement claim additionally
requires matching the relaxed variable, torque observable, geometry, and
orbit assumptions. PENTRC's harmonic factor `1+h^2/4` acts on orbit Fourier
harmonics (`pentrc/energy.f90`, harmonic branch); `h` is not a pitch Legendre
index. In bounce angle it would require `1-(1/4)*d^2/dtheta_b^2`, with an
orbit-dependent coordinate and normalization. NEO-2's scalar weak pitch mass
(`collision_model_mod::krook_pitch_mass`) contains no orbit-harmonic index or
bounce-angle derivative. No additional scalar local factor reproduces this
harmonic dependence, and it is not implemented here.
Unknown model names are rejected.

The registered `collision_krook_test` checks geometry moments against direct
analytic integrals of global pitch polynomials, including bands crossing the
boundary, full trapped/passing bands, field reversal, exact common-factor
scaling and bitwise identity with the legacy unit-factor pitch mass, for both
Krook laws. It also checks actual native D/T/e reaction
matrices against analytic Laguerre/Gamma moments, constant-to-cubic pitch
weak forms including a turning-point band, and manufactured two-energy-mode
complex resolvents. Input rejection tests cover incompatible and invalid
controls. Native full suites, golden comparisons, and scientific scans remain
separate qualification steps.

Mathematical leaves are generated by `gen/gen_krook.f90` using FortSym;
`COMMON/generated/krook_leaf.inc` records the generator revision and command.
The geometry moments are generated separately by `gen/gen_krook_geometry.f90`
in `COMMON/generated/krook_geometry_leaf.inc`, preserving the existing leaves.
