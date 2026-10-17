# Direct magnetic work from the native kinetic response

`neo2_multispecies_out.h5` additionally reports `TphiEM_spec` and `TphiEM_tot`
for the Arnoldi Order 2 solver (`isw_ripple_solver=3`) in native lab 8/9 Boozer
coordinates with positive toroidal field. The datasets are present only when
`TphiEM_geometry_qualified=1`. This flag certifies the implemented coordinate
scope; it does not certify numerical convergence or agreement between models.
Existing `TphiNA_*`, particle fluxes, solves and collision operators retain
their definitions.

NEO-2 returns arrays ordered `(s,phi,theta)`. With `P=abs(Phi_edge)/(2*pi)`,
`C=G+iota*I`, its lab 8/9 geometry gives `sqrtg=-C*P/B^2`,
`isqrg=B^2/(C*P)=-1/sqrtg` and `hphi=B/C` before conversion to CGS.
Thus the supported positive `hphi`, positive `G` branch has a negative
returned Jacobian. The predicate requires these signs and positive `psi_pr_hat`,
positive `g,L_B` and nonzero `iota`. Other formats and reversed field orientations
remain outside the implemented scope. The line measure uses `abs(hphi)`;
neither it nor the physical work requires an additional Jacobian sign.

The new observable is the direct magnetic work

\[
T_{{\rm EM},a}=-\left\langle\int d^3v\,f_{11,a}\,
\partial_{\varphi_0}H_1\right\rangle,\qquad
H_1=T_a u_a^2(1+\lambda^2)b_n .
\]

Here `bnoverb0` is the complex amplitude of the relative field perturbation,
including its native field-line phase. Real physical fields are reconstructed
with `Re(amplitude*exp(i*n*phi0))`. Consequently the phase average is
`-n*Im(f11*conjg(H1))/2`. The remaining phase must not be applied twice.
This is the physical-velocity non-Maxwellian response of Kasilov et al.,
*Physics of Plasmas* **21**, 092506 (2014), equations (50), (59) and (68).
It is not an assertion that every full-domain native response is the same
canonical action-angle nonadiabatic variable in another code.

For native receiving species `a`, driving species `b` and source channels
`(1,2,3)=(A1,A3,A2)`, let `Y` denote the solved pitch-bin integral. The physical
distribution is `F_Ma` times its native energy/pitch reconstruction with weights

\[
\alpha_{1b}=\rho_b A_{1b}/g,\quad
\alpha_{2b}=B_{\rm ref} A_{3b},\quad
\alpha_{3b}=\rho_b A_{2b}/g,\quad
\rho_b=\frac{m_b v_{Tb}c}{Z_b e B_{\rm ref}},\quad g=\langle|\nabla s|\rangle.
\]

`A2=g*dT_ds/T`, `A1=g*dn_ds/n-1.5*A2-Z*e*Er/T`, and
`A3=Z*e*<Epar*B>/<B^2>/T`. The signed driving-species Larmor radius is required.
The forward-source primitive `|lambda|*(4/Bhat-eta)/(3*Bhat*hphi)` has derivative
`-(1+lambda^2)/(2*Bhat*|lambda|*hphi)`, fixing these weights independently of the
transport-output normalization. See also `ripple_solver_normalizations_and_output.tex`,
the multispecies DKE source and distribution reconstruction equations.

The native back-bra contains `i*n*scale*b_n*Bhat`, the pitch moment
`(1+lambda^2)/(Bhat^2*hphi)` and the energy moment
`integral du u^4 exp(-u^2) phi_m/sqrt(pi)`. Its contraction already includes
the real-field phase factor one half. With `L_B=y(6)=integral dl/Bhat` and
`scale=S/(iota*psi_pr_hat)`, where the native solver uses
`S=sign(psi_pr_hat*Bphi_hat*isqrg)=+1` in the admitted orientation, this gives

\[
T_{{\rm EM},a}=-\frac{2n_a T_a}{\mathrm{scale}\,L_B}
\sum_{b,k}\alpha_{kb}\,\mathrm{qflux}(1,k,a,b).
\]

`TphiEM_source_weight`, `TphiEM_qflux_row1`, `TphiEM_dl_over_b`,
`TphiEM_scale_kG`, and `TphiEM_epar_over_b2` retain these consumed quantities.
CGS output is `erg/cm^3`; multiply by 0.1 for `N m/m^3`.

No momentum-conservation assumption enters this work bra. A finite Krook
operator can exchange momentum with its relaxation reservoir. Therefore the
retained nonaxisymmetric radial moment is not automatically the complete
physical radial particle flux; omitted collisional momentum moments must be
retained when using a momentum balance to infer that flux.

The legacy radial transport conversion also contains `r_B=y(14)/y(13)` on
each radial leg. For radial drives only, `TphiNA=r_B^2*TphiEM` in the admitted
orientation. The parallel-drive contribution has no such factor. In Boozer
coordinates the converged ratio estimates `<B^3>/<B^2>/Bref=B00/Bref=1` for
the native `Bref=B00`. Finite field-line length and starting phase can bias
this double-integral estimator. `TphiEM_b00_estimate_ratio` reports it directly;
its deviation is a geometry-convergence diagnostic, not a collision correction.
