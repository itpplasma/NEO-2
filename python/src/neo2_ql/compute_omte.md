# Computing $\Omega_{tE}$ from Radial Force Balance

## Overview

The $E \times B$ rotation frequency $\Omega_{tE}$ characterizes the toroidal
precession of plasma due to the radial electric field $E_r$.  This module
provides a hierarchy of models for computing $\Omega_{tE}$ from experimental
profiles without running the full neoclassical drift-kinetic solver in NEO-2.

The models are implemented in
[`compute_omte.py`](compute_omte.py)
with validation tests in
[`test_compute_omte.py`](../../test/test_compute_omte.py)
and a reference fixture extracted from an AUG #30835 NEO-2 run in
[`omte_reference_aug30835.npz`](../../test/data/omte_reference_aug30835.npz).

GitHub tracking: umbrella issue
[#75](https://github.com/itpplasma/NEO-2/issues/75).


## Physics

### Radial force balance

The radial force balance for species $a$ in a magnetised plasma, neglecting
inertia and viscosity, reads (Gaussian CGS units)

$$
Z_a e\, n_a \!\left( \mathbf{E} + \frac{1}{c}\, \mathbf{v}_a \times \mathbf{B} \right)
= \nabla p_a
\tag{1}
$$

where $Z_a$ is the charge number, $e$ the elementary charge
($4.803 \times 10^{-10}\,\text{statC}$), $n_a$ the density, $p_a = n_a T_a$
the scalar pressure, $\mathbf{v}_a$ the mean flow velocity, and $c$ the speed
of light.  This is equation (12.1) in Helander & Sigmar [1].

Taking the component along $\nabla r$ (where $r$ is the effective minor radius)
and decomposing the flow into toroidal and poloidal components:

$$
E_r = \frac{1}{Z_a e\, n_a}\, \frac{\mathrm{d}p_a}{\mathrm{d}r}
    + v_{\varphi,a}\, \frac{B_\vartheta}{c}
    - v_{\vartheta,a}\, \frac{B_\varphi}{c}
\tag{2}
$$

Since $E_r$ is species-independent, any single species suffices.  In practice
the main ion species is used because its profiles (from CXRS) are best
characterised.

This is the standard form given in Ida & Fujita [2] eq. (1), Kim, Diamond &
Groebner [3] eq. (1), and Viezzer et al. [4] eq. (1).

### $E \times B$ rotation frequency

In Boozer flux coordinates with toroidal flux label $s = \Phi_\text{tor} /
\Phi_\text{tor,edge}$, the $E \times B$ toroidal precession frequency is

$$
\Omega_{tE} = \frac{c\, E_r}{\iota\, \sqrt{g}\, B^\varphi}
\tag{3}
$$

where $\iota = 1/q$ is the rotational transform, $\sqrt{g}\, B^\varphi$ is the
Jacobian times the contravariant toroidal magnetic field component (in
$\text{G}\,\text{cm}$), and $E_r$ is the radial electric field with respect
to the effective radius.

This definition is the one used internally in NEO-2
([`ntv_mod.f90`](../../../NEO-2-QL/ntv_mod.f90) line 2315).


## Model hierarchy

### Level 0: Diamagnetic (pressure gradient only)

Set $v_{\varphi,a} = v_{\vartheta,a} = 0$ in eq. (2):

$$
E_r^{(\text{dia})} = \frac{1}{Z_i e\, n_i}\, \frac{\mathrm{d}p_i}{\mathrm{d}r}
\tag{4}
$$

with the pressure gradient

$$
\frac{\mathrm{d}p_i}{\mathrm{d}r}
= \left( T_i \frac{\mathrm{d}n_i}{\mathrm{d}s}
       + n_i \frac{\mathrm{d}T_i}{\mathrm{d}s} \right)
  \langle |\nabla s| \rangle
\tag{5}
$$

where $\langle |\nabla s| \rangle$ (`av_nabla_stor`) converts from the flux
coordinate derivative $\mathrm{d}/\mathrm{d}s$ to the effective-radius
derivative $\mathrm{d}/\mathrm{d}r$.

Substituting into eq. (3):

$$
\Omega_{tE}^{(\text{dia})}
= \frac{c}{\iota\, \sqrt{g}\, B^\varphi}
  \cdot \frac{\langle |\nabla s| \rangle}{Z_i e\, n_i}
  \left( T_i \frac{\mathrm{d}n_i}{\mathrm{d}s}
       + n_i \frac{\mathrm{d}T_i}{\mathrm{d}s} \right)
\tag{6}
$$

**Implementation:**
[`compute_omte_diamagnetic()`](compute_omte.py)

**GitHub issue:**
[#72](https://github.com/itpplasma/NEO-2/issues/72).

**Inputs needed:** $n_i(s)$, $T_i(s)$, their derivatives $\mathrm{d}n_i/\mathrm{d}s$,
$\mathrm{d}T_i/\mathrm{d}s$, charge number $Z_i$, and the geometry quantities
$\iota$, $\sqrt{g} B^\varphi$, $\langle|\nabla s|\rangle$ from a Boozer-coordinate
equilibrium or a previous NEO-2 run.

### Level 1: Measured toroidal rotation

Retain the $v_\varphi B_\vartheta$ term in eq. (2), using the experimentally
measured toroidal rotation velocity from charge-exchange recombination
spectroscopy (CXRS).  The poloidal velocity is still neglected
($v_\vartheta = 0$):

$$
E_r^{(1)} = \frac{1}{Z_i e\, n_i}\, \frac{\mathrm{d}p_i}{\mathrm{d}r}
           + v_{\varphi,i}\, \frac{B_\vartheta}{c}
\tag{7}
$$

This is the standard approach used in experimental $E_r$ determination; see
Viezzer et al. [4] section 2.

The Python implementation is
[`compute_omte_toroidal_rotation()`](compute_omte.py),
with the unified entry point
[`compute_omte_force_balance()`](compute_omte.py)
falling back to Level 0 when `v_phi` and `b_theta` are omitted.

**GitHub issue:**
[#73](https://github.com/itpplasma/NEO-2/issues/73).

### Level 2: Neoclassical poloidal velocity estimate

Include all three terms in eq. (2), estimating $v_{\vartheta,i}$ from
neoclassical theory.  In the banana regime the poloidal velocity is (Kim,
Diamond & Groebner [3] eq. (6)):

$$
v_{\vartheta,i} \approx K_i(\nu^*)\, \frac{1}{Z_i e B_\varphi}\,
                        \frac{\mathrm{d}T_i}{\mathrm{d}r}
\tag{8}
$$

where $K_i$ depends on the collisionality regime:

| Regime | $\nu^*$ range | $K_i$ |
|--------|--------------|-------|
| Banana | $\nu^* \ll 1$ | $-1.17$ |
| Plateau | $\nu^* \sim 1$ | $\approx -0.5$ |
| Pfirsch--Schluter | $\nu^* \gg 1$ | $\approx +0.5$ |

The Python implementation is
[`compute_omte_neoclassical_poloidal()`](compute_omte.py),
which uses
[`compute_poloidal_rotation_neoclassical()`](compute_omte.py)
with a manually selected coefficient `K_i`.
For convenience,
[`compute_omte_neoclassical_poloidal_auto_k()`](compute_omte.py)
selects `K_i` from a simple collisionality map:

| Regime | $\nu^*$ range | $K_i$ |
|--------|--------------|-------|
| Banana | $\nu^* < 0.1$ | $-1.17$ |
| Plateau | $0.1 \le \nu^* < 10$ | $-0.5$ |
| Pfirsch--Schluter | $\nu^* \ge 10$ | $+0.5$ |

**GitHub issue:**
[#74](https://github.com/itpplasma/NEO-2/issues/74).

### Auxiliary model: Strict NEO-2 `Vphi` convention

The solver also supports a repo-native interpretation of the toroidal
rotation input through the reduced `isw_Vphi_loc=0` algebra in
[`compute_Er()`](../../../NEO-2-QL/ntv_mod.f90).  If the neoclassical
transport terms are suppressed but the `Vphi` term is kept exactly as written
in the Fortran, the resulting model is

$$
E_r^{(\text{exact }V_\varphi)} =
\frac{
  V_\varphi (\iota B_\vartheta + B_\varphi)
  + \frac{c\, T_i B_\vartheta}{Z_i e\, \sqrt{g} B^\varphi}
    \frac{1}{p_i}\frac{\mathrm{d}p_i}{\mathrm{d}r}
}{
  \frac{c\, B_\vartheta}{\sqrt{g} B^\varphi}
}
\tag{9}
$$

This is implemented by
[`compute_omte_toroidal_rotation_neo2_convention()`](compute_omte.py).
It is useful for checking consistency with the code path in NEO-2, but it is
not the same thing as the Level 1 component-pair model above and it should not be
interpreted as the next point in the experimental-model hierarchy.

### Level 3: Full neoclassical (NEO-2 internal)

NEO-2 computes $E_r$ self-consistently from the multi-species drift-kinetic
equation by solving for the parallel flows via the axisymmetric transport
coefficients $D_{31}$, $D_{32}$, $D_{33}$ (Kernbichler et al. [5]).

The full formula (for `isw_Vphi_loc=0`, see
[`ntv_mod.f90`](../../../NEO-2-QL/ntv_mod.f90) lines 2239--2315)
has the structure

$$
E_r = \frac{
  V_\varphi (\iota\, B_\vartheta + B_\varphi)
  + \frac{c\, T_i B_\vartheta}{Z_i e\, \sqrt{g} B^\varphi}
    \frac{1}{p_i}\frac{\mathrm{d}p_i}{\mathrm{d}r}
  + \text{(neoclassical terms)}
}{
  \frac{c\, B_\vartheta}{\sqrt{g} B^\varphi}
  + \text{(neoclassical terms)}
}
\tag{10}
$$

where the neoclassical terms involve sums over species of
$D_{31}$, $D_{32}$, $D_{33}$ transport coefficients and the inductive electric
field.  This is activated with `isw_calc_Er=1` in the NEO-2 input.

The mode `isw_calc_Er=2` (added in the `om-te-profile-support` branch,
[#69](https://github.com/itpplasma/NEO-2/issues/69)) bypasses the internal
calculation and reads an externally provided $\Omega_{tE}(s)$ profile, which is
the intended consumer of the output from this module.

### Output-backed exact mode

For debugging and validation there are now two output-backed Python modes in
[`neo2_output_omte.py`](neo2_output_omte.py):

- `compute_omte_from_neo2_output(..., mode="stored")`
  reads the converged `Er` written by NEO-2 and converts it to $\Omega_{tE}$.
- `compute_omte_from_neo2_output(..., mode="transport")`
  reconstructs `Er` from the stored `D31_AX`, `D32_AX`, `D33_AX`,
  `avEparB_ov_avb2`, geometry, and species state for the tested
  `isw_Vphi_loc=0` branch.

The second path mirrors the Fortran `compute_Er()` algebra directly and is the
right bridge between the reduced Python models and the full NEO-2 solve.
It is not a reduced model.  It is an exact reconstruction path from NEO-2
output data.


## Coordinate system and units

All computations use **Gaussian CGS** units, matching NEO-2 internals.
In Gaussian CGS the electric and magnetic fields have the same dimensions:
$[\text{statV/cm}] = [\text{G}] = [\text{g}^{1/2}\,\text{cm}^{-1/2}\,\text{s}^{-1}]$.

| Quantity | Symbol | Unit |
|----------|--------|------|
| Density | $n$ | $\text{cm}^{-3}$ |
| Temperature | $T$ | $\text{erg}$ |
| Elementary charge | $e$ | $4.803 \times 10^{-10}\,\text{statC}$ |
| Speed of light | $c$ | $2.998 \times 10^{10}\,\text{cm/s}$ |
| Magnetic field | $B$ | $\text{G}$ (Gauss) |
| Radial electric field | $E_r$ | $\text{statV/cm}$ |
| $\Omega_{tE}$ | | $\text{rad/s}$ |
| $\sqrt{g}\, B^\varphi$ | | $\text{G}\,\text{cm}$ |
| $\langle|\nabla s|\rangle$ | `av_nabla_stor` | $\text{cm}^{-1}$ |
| $\iota$ | `aiota` | dimensionless |

The flux surface label is $s = s_\text{tor} = \Phi_\text{tor} /
\Phi_\text{tor,edge}$ (normalised toroidal flux).  Profile derivatives
$\mathrm{d}/\mathrm{d}s$ are with respect to this label.

**Unit conversions** commonly needed:

- $T\,[\text{eV}] \to T\,[\text{erg}]$: multiply by $1.602 \times 10^{-12}$
- $n\,[\text{m}^{-3}] \to n\,[\text{cm}^{-3}]$: multiply by $10^{-6}$
- $B\,[\text{T}] \to B\,[\text{G}]$: multiply by $10^4$

These conversions are handled by
[`load_profile_data.py`](load_profile_data.py):
`convert_units_from_si_to_cgs()`.


## Dimensional analysis

The unit chain for $\Omega_{tE} = c\, E_r / (\iota\, \sqrt{g}\, B^\varphi)$:

$$
\frac{\mathrm{d}p}{\mathrm{d}r}
= \underbrace{\frac{\mathrm{d}p}{\mathrm{d}s}}_{\text{erg/cm}^3}
  \times \underbrace{\langle|\nabla s|\rangle}_{\text{cm}^{-1}}
\quad [\text{erg/cm}^4]
$$

$$
E_r = \frac{\mathrm{d}p/\mathrm{d}r}{n\, Z\, e}
\quad \frac{[\text{erg/cm}^4]}{[\text{cm}^{-3}][\text{statC}]}
= [\text{statV/cm}] = [\text{G}]
$$

$$
\Omega_{tE} = \frac{c\, E_r}{\iota\, \sqrt{g}\, B^\varphi}
\quad \frac{[\text{cm/s}][\text{G}]}{[\text{G}\,\text{cm}]}
= [\text{s}^{-1}] \checkmark
$$

The key point is that $\sqrt{g}\, B^\varphi$ has units $[\text{G}\,\text{cm}]$,
not $[\text{G}\,\text{cm}^2]$.  This follows from the flux relation
$\iota\, \sqrt{g}\, B^\varphi = \langle|\nabla s|\rangle \cdot \psi_t'$
where $\psi_t' = \mathrm{d}(\Phi_\text{tor}/2\pi)/\mathrm{d}s$
$[\text{G}\,\text{cm}^2]$.


## Input file conventions

The profile generator
[`generate_multispec_input.py`](generate_multispec_input.py)
writes the radial input datasets
`T_prof`, `n_prof`, `dT_ov_ds_prof`, and `dn_ov_ds_prof`
to the multispecies HDF5 file. These arrays are indexed by
radial point and species and are the authoritative source for any
post-processing that wants to reproduce the run inputs.

In `NEO-2-QL/neo2.f90`, `prepare_multispecies_scan()` reads those HDF5 arrays,
selects one radial point `ind_boozer_s`, and fills the namelist vectors
`T_vec`, `n_vec`, `dT_vec_ov_ds`, and `dn_vec_ov_ds` for that single run.
Later in `neo2.f90` those vectors are copied into the module arrays
`T_spec`, `n_spec`, `dT_spec_ov_ds`, and `dn_spec_ov_ds`, which represent the
state of the currently active surface inside the solver.

Finally, `write_multispec_output_a()` in
[`ntv_mod.f90`](../../../NEO-2-QL/ntv_mod.f90)
writes `T_spec` and `n_spec` back to `neo2_multispecies_out.h5`.
Those fields are therefore solver state for a completed single-surface run,
not replacements for the original radial profile arrays.

The same output writer now also stores the active reconstruction inputs
`dn_spec_ov_ds`, `dT_spec_ov_ds`, `species_tag_Vphi`, `isw_Vphi_loc`, and
`Vphi` in `neo2_multispecies_out.h5`. That makes the output file
self-contained for replaying the exact `compute_Er()` algebra in Python
without reopening the original input deck.

The rotation input `Vphi` written by
[`generate_multispec_input.py`](generate_multispec_input.py)
is stored with the HDF5 unit attribute `rad / s`. In the Kasilov 2014
notation and in the NEO-2 force-balance implementation this is the
contravariant toroidal angular frequency, not a cylindrical velocity.
Accordingly, Level 1 and Level 2 in this Python module accept the
NEO-2-native Boozer component pairs directly:

$$
\Delta E_r^{(1)} = \frac{V^\varphi B_\vartheta^\text{cov}}{c},
\qquad
\Delta E_r^{(2)} = -\frac{V^\vartheta B_\varphi^\text{cov}}{c}.
$$

The helper functions are written so that either

- physical cylindrical velocity/field pairs, or
- Boozer contravariant/covariant pairs

can be supplied, because only the products entering force balance matter.

The reduced repo-native alternative is
[`compute_omte_toroidal_rotation_neo2_convention()`](compute_omte.py),
which preserves the reduced `compute_Er()` algebra without the transport
terms. On the AUG fixture below, that strict form is
much farther from the full NEO-2 result because it still omits the transport
terms in both numerator and denominator.


## Full-output reconstruction

[`compute_omte_from_neo2_output()`](neo2_output_omte.py) has two distinct modes:

- `mode="stored"` reads the `Er` written by NEO-2 and converts it directly to
  $\Omega_{tE}$.
- `mode="transport"` replays the full `compute_Er()` force-balance algebra from
  the stored `D31/D32/D33` coefficients, `avEparB_ov_avb2`, the active species
  gradients, and the measured rotation selector written into
  `neo2_multispecies_out.h5`.

These are exact replay paths from NEO-2 output, not additional reduced-model
levels.

For a current-output file from the patched writer, those two modes must agree
to numerical roundoff. The circular-tokamak regression fixture used by the
tests is regenerated with
[`regenerate_neo2_ql_axisymmetric_fixture.sh`](../../test/data/regenerate_neo2_ql_axisymmetric_fixture.sh),
which reruns a single-surface multispecies axisymmetric case and refreshes
[`neo2_ql_axisymmetric_multispecies_out.h5`](../../test/data/neo2_ql_axisymmetric_multispecies_out.h5).


## Validation against NEO-2

The Level 0 model is validated against full neoclassical NEO-2 output for
ASDEX Upgrade shot #30835 (2 flux surfaces with `isw_calc_Er=1`).
The comparison must use the actual per-surface profile input consumed by the
run (`n_prof`, `T_prof`, `dn_ov_ds_prof`, `dT_ov_ds_prof`). The scalar
species state written to `n_spec` and `T_spec` in the multispecies output is
not the radial profile and produces the wrong Level 0 estimate.

| $s_\text{tor}$ | NEO-2 $\Omega_{tE}$ | Level 0 $\Omega_{tE}$ | Ratio |
|:-:|:-:|:-:|:-:|
| 0.2527 | $-44.80\,\text{krad/s}$ | $-14.81\,\text{krad/s}$ | 0.33 |
| 0.4984 | $-111.09\,\text{krad/s}$ | $-16.74\,\text{krad/s}$ | 0.15 |

Level 0 **underestimates** $|\Omega_{tE}|$ because it omits the toroidal
rotation and neoclassical transport terms that contribute in the full force
balance. The sign is
correct (negative $\Omega_{tE}$, corresponding to inward-pointing $E_r$).

Test: [`test_diamagnetic_vs_neo2_sign_and_order_of_magnitude()`](../../test/test_compute_omte.py).

Using the NEO-2 pair `(Vphi, bcovar_tht)` gives a more complete Level 1
curve:

| $s_\text{tor}$ | NEO-2 $\Omega_{tE}$ | Level 1 $\Omega_{tE}$ | Ratio |
|:-:|:-:|:-:|:-:|
| 0.2527 | $-44.80\,\text{krad/s}$ | $-78.54\,\text{krad/s}$ | 1.75 |
| 0.4984 | $-111.09\,\text{krad/s}$ | $-82.83\,\text{krad/s}$ | 0.75 |

For this AUG reference, Level 1 roughly halves the mean absolute
error relative to Level 0, although it is still missing the poloidal-flow and
transport terms from the full NEO-2 solve.

Using the simple banana-regime Level 2 estimate with `K_i = -1.17`,
or equivalently the auto-selected low-collisionality branch of
[`compute_omte_neoclassical_poloidal_auto_k()`](compute_omte.py),
does not visibly move the curve for this reference:

| $s_\text{tor}$ | NEO-2 $\Omega_{tE}$ | Level 2 $\Omega_{tE}$ | Ratio |
|:-:|:-:|:-:|:-:|
| 0.2527 | $-44.80\,\text{krad/s}$ | $-78.54\,\text{krad/s}$ | 1.75 |
| 0.4984 | $-111.09\,\text{krad/s}$ | $-82.83\,\text{krad/s}$ | 0.75 |

That near-overlap is expected from the simple analytic estimate because
the `-v_\vartheta B_\varphi/c` contribution reduces to
`-K_i (dT_i/dr)/(Z_i e c)` and is very small here compared to the toroidal
rotation contribution.

Applying the strict reduced `isw_Vphi_loc=0` algebra gives a very different
curve:

| $s_\text{tor}$ | NEO-2 $\Omega_{tE}$ | Strict `Vphi` convention $\Omega_{tE}$ | Ratio |
|:-:|:-:|:-:|:-:|
| 0.2527 | $-44.80\,\text{krad/s}$ | $+17087.17\,\text{krad/s}$ | -381.47 |
| 0.4984 | $-111.09\,\text{krad/s}$ | $+11827.27\,\text{krad/s}$ | -106.46 |

This is not a better reduced model.  It is a useful diagnostic because it
shows that the bare `Vphi` term in the Fortran algebra is not sufficient by
itself: once the transport-coefficient pieces in eq. (10) are removed, the
remaining exact-convention term overshoots by two orders of magnitude and even
flips sign on this case.


## Algebraic consistency with NEO-2 Fortran

The Level 0 Python formula is algebraically identical to the NEO-2 Fortran
`compute_Er()` subroutine when the neoclassical transport coefficients and
toroidal rotation are set to zero
($D_{31} = D_{32} = D_{33} = 0$, $V_\varphi = 0$).

Starting from the Fortran (lines 2268--2272 of
[`ntv_mod.f90`](../../../NEO-2-QL/ntv_mod.f90)):

$$
\text{nom} = \frac{c\, T_i\, B_\vartheta}{Z_i e\, \sqrt{g} B^\varphi}
             \cdot \frac{1}{p_i}\frac{\mathrm{d}p_i}{\mathrm{d}r},
\qquad
\text{denom} = \frac{c\, B_\vartheta}{\sqrt{g} B^\varphi}
$$

$$
E_r = \frac{\text{nom}}{\text{denom}}
    = \frac{T_i}{Z_i e} \cdot \frac{1}{p_i}\frac{\mathrm{d}p_i}{\mathrm{d}r}
    = \frac{1}{Z_i e\, n_i}\frac{\mathrm{d}p_i}{\mathrm{d}r}
$$

which is exactly eq. (4).  The $B_\vartheta$ and $\sqrt{g} B^\varphi$ factors
cancel, so the Level 0 result depends only on thermodynamic profiles and
$\langle|\nabla s|\rangle$, $\iota$, $\sqrt{g} B^\varphi$.


## Internal code references

| File | Role |
|------|------|
| [`compute_omte.py`](compute_omte.py) | Python implementation of force balance models |
| [`neo2_output_omte.py`](neo2_output_omte.py) | Stored and reconstructed Om_tE from NEO-2 output |
| [`plot_omte_reference.py`](plot_omte_reference.py) | Reproducible AUG comparison plot helper |
| [`test_compute_omte.py`](../../test/test_compute_omte.py) | Unit + e2e tests |
| [`omte_reference_aug30835.npz`](../../test/data/omte_reference_aug30835.npz) | Reference fixture (AUG #30835) |
| [`generate_multispec_input.py`](generate_multispec_input.py) | Writes `Om_tE` profile to multispec HDF5 input |
| [`ntv_mod.f90`](../../../NEO-2-QL/ntv_mod.f90) | Fortran `compute_Er()` subroutine (lines 2038--2321) |
| [`neo2.f90`](../../../NEO-2-QL/neo2.f90) | Reads `isw_calc_Er`, `Om_tE_prof` from input |


## External references

[1] P. Helander and D. J. Sigmar,
*Collisional Transport in Magnetized Plasmas*,
Cambridge University Press (2002).
Chapter 12: Radial electric field and rotation.

[2] K. Ida and T. Fujita,
"Internal transport barrier in tokamak and helical plasmas,"
*Plasma Phys. Control. Fusion* **60**, 033001 (2018).
[doi:10.1088/1361-6587/aa9b03](https://doi.org/10.1088/1361-6587/aa9b03)

[3] Y. B. Kim, P. H. Diamond, and R. J. Groebner,
"Neoclassical poloidal and toroidal rotation in tokamaks,"
*Phys. Fluids B* **3**, 2050 (1991).
[doi:10.1063/1.859671](https://doi.org/10.1063/1.859671)

[4] E. Viezzer, T. Putterich, R. Dux, R. M. McDermott, and the ASDEX Upgrade Team,
"High-accuracy characterization of the edge radial electric field at ASDEX Upgrade,"
*Nucl. Fusion* **52**, 043011 (2012).
[doi:10.1088/0029-5515/52/4/043011](https://doi.org/10.1088/0029-5515/52/4/043011)

[5] W. Kernbichler, S. V. Kasilov, G. O. Leitold, V. V. Nemov, and K. Allmaier,
"Recent progress in NEO-2 -- A code for neoclassical transport computations based on field line tracing,"
*Plasma Fusion Res.* **3**, S1061 (2008).
[doi:10.1585/pfr.3.S1061](https://doi.org/10.1585/pfr.3.S1061)

[6] S. V. Kasilov, W. Kernbichler, A. F. Martitsch, H. Maassberg, and M. F. Heyn,
"Evaluation of the toroidal torque driven by external non-resonant non-axisymmetric
magnetic field perturbations in a tokamak,"
*Phys. Plasmas* **21**, 092506 (2014).
[doi:10.1063/1.4894479](https://doi.org/10.1063/1.4894479)
