# Complete CGS covariant curl

The exercised `magfie_result=0` field interface uses coordinates
`(s,phi,theta)`. Its covariant unit field is `h=100*(C,G,I)/B`, where
`C=bcovar_s`, `G=curr_pol`, `I=curr_tor`, and B and the covariant field data
are in the existing SI convention before conversion. G and I are flux
functions. For the final signed CGS Jacobian J returned by the interface,
`curl(h)^i = epsilon^(ijk)*partial_j(h_k)/J` gives

```
curl1 = 100*(G*B_theta-I*B_phi)/(B**2*J)
curl2 = 100*(B*(C_theta-I_s)-C*B_theta+I*B_s)/(B**2*J)
curl3 = 100*(B*(G_s-C_phi)-G*B_s+C*B_phi)/(B**2*J)
```

The existing metric calculation supplies C_theta and C_phi. C_s is not needed.
For the CGS output, B_theta/B_phi and the metric angular partials come from
`eva2d_fd` applied to the SAME value coefficient arrays used for returned B,
gstb and gspb. The existing C chain rule composes these derivatives. Previously,
independent interpolants of analytic partials differed from derivatives of the
returned value polynomials. Primal values and radial partials are unchanged;
angular bder and curl can change. Consistency does not establish physical
accuracy or preserve historical scientific outputs at a fixed grid. Full
application golden qualification is a separate open gate.

The generated leaf uses the final signed J directly; the reduced formula's
extra orientation correction must not be applied twice. All other field units,
current conventions, and interpolants remain those of the existing provider.
Finite inputs, B>0 and J nonzero are required. See [the generator](../gen/README.md)
for the authoritative retained definition and optional regeneration recipe.

This statement covers the CGS interface only. The unused legacy
`magfie_result=1` output has a different coordinate order and explicitly zero
radial covariant component. It is unchanged and not qualified by these tests.
The radial derivative returned by `neo_magfie_b` is still unavailable; a curl
needs only the two angular derivatives of the radial covariant component.

The explicit provider remains a single retained global instance. Callers must
serialize its lifecycle and evaluations. `neo_field_load` accepts optional
`theta_points` and `phi_points` (each >=5) to select the existing cubic angular
interpolation grids before allocating caches. Omitting them preserves 1025/9.
Resolution is an interpolation choice, not an automatic accuracy guarantee.
Reloading a file/grid invalidates all old surface caches.

## Frozen independent regression

The complete-curl test uses seven radial surfaces s=.2,.3,...,.8, R0=5m,
a=.5m, B0=2T (reload B0=3T), r=a*sqrt(s), and the geometry

```
R=R0+r*cos(theta)+.18*r*cos(theta-phi)
Z=sigma_theta*(r*sin(theta)+.11*r*sin(theta-phi))
iota=sigma_iota*.2*(1+.1*s)
G=B0*R0*(1+.05*s)
I=B0*a**2*iota*s/R0
B=B0*sqrt(1+(r*iota/R0)**2)*
  (1-.12*sqrt(s)*cos(theta)+.06*sqrt(s)*sin(theta-phi))
```

The cylindrical angle is phi. The .bc format includes the required zero-padded
mode rectangle in addition to the three nonzero modes. Format9 uses
`m*theta+n*phi`, so the theta-phi harmonic is serialized with n=-1. The unequal geometric
harmonics make the radial covariant component and its partials material.
The fixture is a provider consistency oracle, not a force-balanced equilibrium.

At s=.5, the frozen (phi,theta) pairs are (.37,.79), (1.21,2.11), (4.33,5.27),
and their nearest angular grid nodes. Both chart orientations and both iota
signs are checked. The angular ladder is 65/129/257 in both directions.
Independent derivatives use the actual returned covariant field, fourth-order
centered differences with h=1e-3 and h/2, and an independent h=5e-4/h/2 pair.
Every component uses `1e-8*max(abs(expected),1e-8)` as its comparison bound.
The omitted radial-component terms must exceed 100 times that bound at a
frozen point. The existing axisymmetric 1e-8 oracle is unchanged.

One radial surface is cached at a time. At the largest grid, 26 double-precision
coefficient arrays of shape (4,4,257,257,1) require 219,811,072 bytes; temporary
angular arrays and the small source data remain additional. The test must stay
below 1GiB, with actual process high-water recorded externally during checks.
The point leaf has bounded scalar scratch and allocates nothing per call.

Qualification requires the final grid and FD-step consistency checks to pass.
If independent interpolated derivatives do not agree within that frozen gate,
report the unresolved interpolation/metric derivative bridge; do not change
the fixture, weaken the gate, or advertise arbitrary-spectrum accuracy.
No deterministic enclosure is supplied by these numerical consistency tests.

NEO-RT's declared zero-radial-field approximation and unavailable toroidal curl
component remain consumer choices. This raw-provider repair does not silently
change those choices or their capability flags.
