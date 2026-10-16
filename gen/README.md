# Covariant curl generator

The ordinary NEO-2 build compiles the committed generated Fortran module. It
needs no FortSym, CAS, generator executable, or Python interpreter.

The optional standalone generator retains the definition of the returned CGS
covariant unit field and its Levi-Civita curl in `(s,phi,theta)` order. It emits
the primal, point-jet JVP and transpose VJP together from that definition.
Its native symbolic engine must prove equality to the quotient-rule expansion
and the contracted JVP/VJP transpose identity;
an UNKNOWN result aborts generation. The domain is finite inputs, B>0 and a
nonzero final signed CGS Jacobian. This is a rational identity on that domain,
not a certified interpolation error bound.

Regenerate from the repository root:

```sh
cmake -S gen -B build-curl-gen -G Ninja
cmake --build build-curl-gen --target gen_covariant_curl -j2
build-curl-gen/gen_covariant_curl COMMON/generated/neo_covariant_curl_generated.f90
git diff -- COMMON/generated/neo_covariant_curl_generated.f90
```

FortSym is pinned to `7e8ef7abdd90c49bb5f75e2c119e7d181656bcc6`. For an
existing source checkout use CMake's `FETCHCONTENT_SOURCE_DIR_FORTSYM`; verify
its revision first. If overriding transitive source directories, they must
match the references in that FortSym revision (FortGen
`42a4c6cce9272a6bc92e0539b21f3ffc0c97a068`, FortPlot
`178e6fd2f0440d557faf01cc36cd6ec2a2ef5244`).

Generated files must be regenerated rather than edited by hand. The independent
point and end-to-end provider tests are registered with CTest as
`covariant_curl_leaf` and `field_provider_complete_curl`. Derivative oracles
are numerical checks; they do not replace the symbolic identity or constitute
a full equilibrium/field interpolation adjoint.
