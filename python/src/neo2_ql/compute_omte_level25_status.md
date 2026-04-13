# Level 2.5 Status, Literature Basis, and Sign Investigation

## Scope

This note documents what was implemented for the first Level 2.5 prototype,
what literature result is actually relevant, and why the current one-column
transport closure is not yet an acceptable physics model.

The work belongs to issue #75 and PR #76.


## What Was Implemented

The current branch adds a minimal transport-aware helper in
[neo2_output_omte.py](neo2_output_omte.py) that evaluates a reduced transport
ansatz on top of the exact NEO-2 `compute_Er()` algebra:

$$
D_{31}^{\mathrm{AX}} = \hat D_{31}\, D_{31,\mathrm{ref}},
$$

$$
D_{32}^{\mathrm{AX}} = \left(\frac{5}{2} - k_{\mathrm{cof}}\right) D_{31}^{\mathrm{AX}},
$$

$$
D_{33}^{\mathrm{AX}} = 0.
$$

The code path is:

- `compute_d31_reference_electron(...)`
- `decompose_neo2_er_k_cof_transport_model(...)`
- `compute_neo2_er_from_k_cof_transport_model(...)`
- `compute_neo2_omte_from_k_cof_transport_model(...)`

This prototype was also added to the AUG reference plotting helper in
[plot_omte_reference.py](plot_omte_reference.py), and focused tests were added in
[test_compute_omte.py](../../test/test_compute_omte.py).

Important: this prototype is a diagnostic scaffold, not a literature-derived
closure. It was useful for probing the exact NEO-2 algebra and for sweeping
$k_{\mathrm{cof}}$, but it should not be presented as a validated Level 2.5
model in its current form.


## Literature Result That Actually Matters

The key literature result is already closely aligned with the NEO-2 Fortran
implementation: Kasilov et al., Phys. Plasmas 21, 092506 (2014), section II,
their equations (5), (6), (10), and (11).

For tokamak neoclassical force balance and parallel flow they write

$$
V_{\parallel}
= -\frac{c B_\varphi}{g B^\vartheta B}
\left(\frac{1}{e_i n_i}\frac{d p_i}{d r} - E_r\right)
+ \frac{c k B B_\varphi}{e_i g B^\vartheta \langle B^2 \rangle}
\frac{dT_i}{dr},
$$

with the corresponding divergence-free poloidal and toroidal velocities

$$
V^\vartheta = \frac{c k B_\varphi}{e_i g \langle B^2 \rangle}\frac{dT_i}{dr},
$$

$$
V^\varphi = \frac{c}{g B^\vartheta}
\left(E_r - \frac{1}{e_i n_i}\frac{d p_i}{d r}\right) + q V^\vartheta.
$$

They also state that, for ion transport coefficients,

$$
D_{31} = \frac{c T_i B_\varphi}{e_i g B^\vartheta},
$$

and therefore

$$
k = \frac{5}{2} - \frac{D_{32}}{D_{31}}.
$$

This is the exact literature origin of the `k_cof = 2.5 - D32/D31` relation
used in the Fortran code and in the Python prototype.


## How This Maps to NEO-2

The relevant Fortran implementation is
[ntv_mod.f90](../../../NEO-2-QL/ntv_mod.f90). For the tested branch
`isw_Vphi_loc = 0`, NEO-2 computes

$$
E_r = \frac{N_E}{D_E},
$$

with

$$
D_E = c\,\frac{B^{\mathrm{cov}}_\vartheta}{\sqrt{g} B^\varphi}
+ \sum_j D_{31,ij}^{\mathrm{AX}} \frac{Z_j e}{T_j},
$$

and

$$
N_E = V^\varphi_i \left(\iota B^{\mathrm{cov}}_\vartheta + B^{\mathrm{cov}}_\varphi\right)
+ \frac{c T_i B^{\mathrm{cov}}_\vartheta}{Z_i e\, \sqrt{g} B^\varphi}
\frac{1}{p_i}\frac{dp_i}{dr}
+ \sum_j N_{31,ij} + \sum_j N_{32,ij} + \sum_j N_{33,ij},
$$

where

$$
N_{31,ij} = \langle |\nabla s| \rangle D_{31,ij}^{\mathrm{AX}}
\left(\frac{1}{n_j}\frac{dn_j}{ds} + \frac{1}{T_j}\frac{dT_j}{ds}\right),
$$

$$
N_{32,ij} = \langle |\nabla s| \rangle
\frac{1}{T_j}\frac{dT_j}{ds}
\left(D_{32,ij}^{\mathrm{AX}} - \frac{5}{2} D_{31,ij}^{\mathrm{AX}}\right).
$$

This is exactly what the Python replay function
`decompose_neo2_er_transport_terms(...)` mirrors.

The normalization for the dimensional transport coefficients is also explicit in
Fortran, but not fully uniform across all code paths. The HDF5 writer stores the
base reference scale

$$
D_{31,\mathrm{ref},00}
= \frac{v_{T0} \rho_0 B^{\mathrm{cov}}_{\varphi,\mathrm{hat}} B_{\mathrm{ref}}}
{2 \sqrt{g} B^\vartheta},
$$

and the species-pair re-normalization

$$
\frac{D_{31,\mathrm{ref},ab}}{D_{31,\mathrm{ref},00}}
= \frac{m_b}{m_0}
\frac{v_{Ta} v_{Tb}}{v_{T0}^2}
\frac{Z_0}{Z_b}.
$$

The dimensional coefficients are then reconstructed by multiplying the stored
normalized arrays with $D_{31,\mathrm{ref},00}$.

This is why reducing the model to a single scalar $\hat D_{31}$ is risky: the
actual NEO-2 transport closure has row and column species structure before the
sum over $j$ is taken.


## Geometry and Sign Caveat

Your suspicion about geometry is partly right, but it needs to be split into two
cases.

If Python is doing exact replay from the stored NEO-2 outputs, then the main
geometry factors are already supplied by the HDF5 or NPZ fixture:

$$
\iota, \quad \sqrt{g} B^\varphi, \quad B^{\mathrm{cov}}_\vartheta, \quad B^{\mathrm{cov}}_\varphi.
$$

In that situation, Python should not apply an extra ad hoc left-handed Boozer
sign correction on top of the stored values.

For the AUG fixtures inspected here, the signs are

$$
\sqrt{g} B^\varphi > 0, \qquad
\sqrt{g} B^\vartheta = \iota \sqrt{g} B^\varphi > 0, \qquad
B^{\mathrm{cov}}_\vartheta < 0, \qquad
B^{\mathrm{cov}}_\varphi < 0.
$$

With electron charge $Z_e = -1$, this gives a positive effective
$D_{31,\mathrm{ref}}$ in the present Python replay. That sign chain is internally
consistent.

If instead one tries to rebuild the geometry from simpler input quantities such
as $q$, poloidal flux, and a guessed Jacobian relation, then the left-handed AUG
branch becomes dangerous. In the Fortran code, NEO-2 explicitly does not use the
usual right-handed relation

$$
\sqrt{g} B^\vartheta \propto \iota\, \psi'_{\mathrm{Boozer}}
$$

for `lab_swi = 10`. Instead it reconstructs $\sqrt{g} B^\vartheta$ from the
direct metric and field components because the simple right-handed formula is not
safe in the left-handed AUG coordinate convention.

So:

$$
	ext{for exact replay: stored geometry should be trusted first,}
$$

but

$$
	ext{for reduced geometry reconstruction: left-handed sign errors are a real risk.}
$$


## D31 Reference Inconsistency

There is also a more concrete normalization issue in the Fortran code.

In one multispecies normalization path, `D31ref00` is formed as

$$
D_{31,\mathrm{ref},00}
= \frac{v_{T0} \rho_0 B^{\mathrm{cov}}_{\varphi,\mathrm{hat}} B_{\mathrm{ref}}^2}
{2 \sqrt{g} B^\vartheta},
$$

while the HDF5 writer later stores `D31ref0` using

$$
D_{31,\mathrm{ref},00}^{\mathrm{(written)}}
= \frac{v_{T0} \rho_0 B^{\mathrm{cov}}_{\varphi,\mathrm{hat}} B_{\mathrm{ref}}}
{2 \sqrt{g} B^\vartheta}.
$$

Numerically, in the committed axisymmetric fixture,

$$
\frac{D_{31}^{\mathrm{AX}}}{D_{31,\mathrm{AX}}/D_{31,\mathrm{ref}}}
= 1.09464 \times 10^{12},
$$

whereas the stored scalar dataset is

$$
D31ref0 = 5.68618 \times 10^7.
$$

Their ratio is exactly

$$
\frac{1.09464 \times 10^{12}}{5.68618 \times 10^7} = B_{\mathrm{ref}} \approx 1.92509 \times 10^4\,\mathrm{G}.
$$

So the mismatch is not random. It is one explicit factor of $B_{\mathrm{ref}}$.

This matters for amplitude normalization, but not for the AUG sign failure by
itself. Multiplying all transport terms by an extra positive factor of
$B_{\mathrm{ref}}$ changes the size of the response, but it does not explain why
the electron-column and ion-column reduced models have opposite signs.


## Full-Geometry $k$ Sweep Check

I also checked the current Level 2 model against the committed AUG fixture using
the actual NEO-2 geometry already stored in the file:

$$
\iota, \quad \sqrt{g} B^\varphi, \quad B^{\mathrm{cov}}_\vartheta, \quad B^{\mathrm{cov}}_\varphi,
\quad \langle |\nabla s| \rangle.
$$

For these data, the Boozer identity is satisfied numerically when interpreted
with respect to the effective radius used by NEO-2:

$$
\sqrt{g} B^\vartheta = \iota \sqrt{g} B^\varphi = \frac{\psi'_{\mathrm{tor}}(r_{\mathrm{eff}})}{q}
= \psi'_{\mathrm{pol}}(r_{\mathrm{eff}}).
$$

So for the stored geometry itself, the problem is not that this identity fails.

However, even if one sweeps the current Python Level 2 model over a very broad
range,

$$
k \in [-2.5,\,1.5],
$$

the resulting envelope is

$$
\Omega_{tE}^{\mathrm{Level\ 2}} \in
\left[-99.1,\,-66.2\right]\,\mathrm{krad/s}
$$

for the inner surface and

$$
\Omega_{tE}^{\mathrm{Level\ 2}} \in
\left[-100.5,\,-72.3\right]\,\mathrm{krad/s}
$$

for the outer surface, while NEO-2 gives

$$
\Omega_{tE}^{\mathrm{NEO2}} \approx
\left[-44.8,\,-111.1\right]\,\mathrm{krad/s}.
$$

Therefore the exact NEO-2 result is outside even this enlarged $k$ envelope on
both surfaces.

This is the strongest current argument for your proposed workflow:

$$
	ext{use the actual NEO-2 geometry first, then test the reduced model against it,}
$$

and only simplify further after that reduced model has been shown to bracket or
track the NEO-2 result in its expected regime.


## What Was Verified in Code

The focused Python test file passes with the current implementation.

The new tests verify:

1. The electron-side `D31` reference reconstructed in Python matches the stored
   raw-to-normalized ratio in the axisymmetric HDF5 fixture.
2. Setting $\hat D_{31} = 0$ collapses exactly to the no-transport replay of the
   reduced `compute_Er()` algebra.
3. Changing $k_{\mathrm{cof}}$ has no effect if $dT/ds = 0$, as required by the
   NEO-2 numerator structure.
4. On the single-surface axisymmetric fixture, the prototype moves the replayed
   $E_r$ in the right direction relative to the no-transport overshoot.

Those checks establish that the Python scaffold is algebraically consistent with
the tested Fortran branch. They do not validate the one-column model as a
physics closure for AUG.


## Sign Investigation on the AUG Reference Fixture

The AUG reference fixture is the real stress test. The NEO-2 target values are

$$
\Omega_{tE}^{\mathrm{NEO2}} \approx \left[-44.8,\,-111.1\right]\,\mathrm{krad/s}.
$$

Using the present one-column prototype with the default sweep values

$$
\hat D_{31} = -1, \qquad k_{\mathrm{cof}} = 0.565,
$$

gives very different results depending on which transport column is retained.

If the retained transport column is the electron column,

$$
\Omega_{tE}^{(j=e)} \approx \left[-310.8,\,-219.5\right]\,\mathrm{krad/s}.
$$

If the retained transport column is the measured-ion column,

$$
\Omega_{tE}^{(j=i)} \approx \left[+195.9,\,+168.7\right]\,\mathrm{krad/s}.
$$

So the sign problem is real, and it is not explained by a trivial statement like
`sqrtg has the wrong sign in Python`.

The stronger conclusion is this:

$$
\text{the sign and magnitude are highly sensitive to the retained transport column,}
$$

which means the one-column closure is too small to represent the species-summed
NEO-2 transport physics.


## Interpretation

There are two separate points.

First, the prototype does use a relation from the literature,

$$
k = \frac{5}{2} - \frac{D_{32}}{D_{31}},
$$

but that alone is not enough to define a valid Level 2.5 closure.

Second, the bad AUG sign is best interpreted as a model-structure failure, not
as proof that the exact replay code has a simple metric-sign mistake. The NEO-2
Fortran sums over transport columns, and those columns carry distinct species
normalizations before they enter the numerator and denominator of $E_r$.

Collapsing that structure to one coefficient pair,

$$
\left(D_{31}, D_{32}\right) \to \left(\hat D_{31}, k_{\mathrm{cof}}\right),
$$

throws away enough information that even the sign can change.


## Present Conclusion

The current one-column Level 2.5 prototype should be treated as a debugging and
exploration tool only.

What it is good for:

- checking that the Python replay matches the tested NEO-2 algebra,
- sweeping $k_{\mathrm{cof}}$ in a controlled way,
- exposing which transport terms and species columns matter for the sign.

What it is not good for:

- presenting a literature-grounded reduced model for AUG,
- claiming physical validation against NEO-2,
- claiming the correct sign is understood.

The next acceptable Level 2.5 step must be literature-first and must retain at
least enough transport structure to respect the species-column dependence that is
visible in the exact NEO-2 sums.


## References

1. S. V. Kasilov, W. Kernbichler, A. F. Martitsch, H. Maassberg, and M. F. Heyn,
   Evaluation of the toroidal torque driven by external non-resonant
   non-axisymmetric magnetic field perturbations in a tokamak,
   Phys. Plasmas 21, 092506 (2014).
   https://doi.org/10.1063/1.4894479
2. F. L. Hinton and R. D. Hazeltine, Theory of plasma transport in toroidal
   confinement systems, Rev. Mod. Phys. 48, 239 (1976).
3. P. Helander and D. J. Sigmar, Collisional Transport in Magnetized Plasmas,
   Cambridge University Press (2002).