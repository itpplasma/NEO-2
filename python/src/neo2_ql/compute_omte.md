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
$\text{G}\,\text{cm}^2$), and $E_r$ is the radial electric field with respect
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
[`compute_omte_diamagnetic()`](compute_omte.py) (lines 24--63).

**GitHub issue:**
[#72](https://github.com/itpplasma/NEO-2/issues/72).

**Inputs needed:** $n_i(s)$, $T_i(s)$, their derivatives $\mathrm{d}n_i/\mathrm{d}s$,
$\mathrm{d}T_i/\mathrm{d}s$, charge number $Z_i$, and the geometry quantities
$\iota$, $\sqrt{g} B^\varphi$, $\langle|\nabla s|\rangle$ from a Boozer-coordinate
equilibrium or a previous NEO-2 run.

### Level 1: Measured toroidal rotation (planned)

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

**GitHub issue:**
[#73](https://github.com/itpplasma/NEO-2/issues/73).

### Level 2: Neoclassical poloidal velocity estimate (planned)

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

**GitHub issue:**
[#74](https://github.com/itpplasma/NEO-2/issues/74).

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
\tag{9}
$$

where the neoclassical terms involve sums over species of
$D_{31}$, $D_{32}$, $D_{33}$ transport coefficients and the inductive electric
field.  This is activated with `isw_calc_Er=1` in the NEO-2 input.

The mode `isw_calc_Er=2` (added in the `om-te-profile-support` branch,
[#69](https://github.com/itpplasma/NEO-2/issues/69)) bypasses the internal
calculation and reads an externally provided $\Omega_{tE}(s)$ profile, which is
the intended consumer of the output from this module.


## Coordinate system and units

All computations use **Gaussian CGS** units, matching NEO-2 internals:

| Quantity | Symbol | Unit |
|----------|--------|------|
| Density | $n$ | $\text{cm}^{-3}$ |
| Temperature | $T$ | $\text{erg}$ |
| Elementary charge | $e$ | $4.803 \times 10^{-10}\,\text{statC}$ |
| Speed of light | $c$ | $2.998 \times 10^{10}\,\text{cm/s}$ |
| Magnetic field | $B$ | $\text{G}$ (Gauss) |
| Radial electric field | $E_r$ | $\text{statV/cm}$ |
| $\Omega_{tE}$ | | $\text{rad/s}$ |
| $\sqrt{g}\, B^\varphi$ | | $\text{G}\,\text{cm}^2$ |
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


## Validation against NEO-2

The Level 0 model is validated against full neoclassical NEO-2 output for
ASDEX Upgrade shot #30835 (2 flux surfaces with `isw_calc_Er=1`).

| $s_\text{tor}$ | NEO-2 $\Omega_{tE}$ | Level 0 $\Omega_{tE}$ | Ratio |
|:-:|:-:|:-:|:-:|
| 0.2527 | $-44.80\,\text{krad/s}$ | $-15.64\,\text{krad/s}$ | 0.35 |
| 0.4984 | $-111.09\,\text{krad/s}$ | $-16.74\,\text{krad/s}$ | 0.15 |

Level 0 captures 15--35% of the full neoclassical result.  The sign is
correct (negative $\Omega_{tE}$, corresponding to inward-pointing $E_r$).
The remainder is dominated by the toroidal rotation contribution (Level 1) and
neoclassical corrections (Levels 2--3).

The decreasing ratio at larger $s_\text{tor}$ is expected: the toroidal
rotation $V_\varphi$ enters the full NEO-2 formula (eq. 9) with a coefficient
$\iota B_\vartheta + B_\varphi$ that grows relative to the diamagnetic term
at larger minor radius.

Test: [`test_diamagnetic_vs_neo2_sign_and_order_of_magnitude()`](../../test/test_compute_omte.py).


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
| [`test_compute_omte.py`](../../test/test_compute_omte.py) | Unit + e2e tests |
| [`omte_reference_aug30835.npz`](../../test/data/omte_reference_aug30835.npz) | Reference fixture (AUG #30835) |
| [`generate_multispec_input.py`](generate_multispec_input.py) | Writes `Om_tE` profile to multispec HDF5 input |
| [`er_rotation_mod.f90`](../../../NEO-2-QL/er_rotation_mod.f90) | Fortran $\Omega_{tE} \leftrightarrow M_t/R$ conversion |
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
