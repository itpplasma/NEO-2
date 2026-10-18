# Small scalar-system capture

Set `NEO2_DUMP_LINEAR_SYSTEM` to an existing empty directory when running one
small Krook surface (at most 4096 unknowns per species). With the variable unset,
no capture is written. The writer refuses to overwrite files; use a fresh directory
for every case. This is a diagnostic of the assembled system, not a solver mode.

The Matrix Market files contain the native general complex matrix, the three
flux bras as columns, and each source RHS and resulting solution. Names identify
axisymmetric/nonaxisymmetric equations and responding species. The nonaxisymmetric
capture also records the geometric source before `add_f01_source` and the complete
source after it. Their difference includes both the collisional and pitch derivative
parts of the native perturbed-f01 source. It does not isolate either part separately.

Run the independent dense oracle on a cluster:

```
python TEST/dense_linear_system_oracle.py capture/na_species0 --output direct.json
```

The oracle solves the captured equation with NumPy's dense solver and compares
the native solution, equation backward error and the native flux contraction
`0.5*Re(conjg(bra)^T solution)`. It also solves the geometric and f01 source
pieces separately and reports their fluxes and linear closure. Fluxes are native
dimensionless qflux values before the caller's physical-unit reconstruction.
They are not already torques. No positivity assumption or physical inner product
is inferred from the coordinates of the native unknown.

A small residual with dense agreement excludes inaccurate solution of this captured
discrete equation as the cause of its flux. It does not validate the discretization,
source, bra, physical collision model, or torque normalization. A reproducible
difference in geometric versus f01 flux identifies where the next independent
physical balance must be tested, without claiming that either term is erroneous.

The focused writer/oracle regression can be run using the CMake-built
`test_linear_system_diagnostic` executable through
`NEO2_DIAGNOSTIC_TEST_EXE=/absolute/path/to/test_linear_system_diagnostic pytest
TEST/test_dense_linear_system_oracle.py`. It includes complex phases, a zero RHS,
and deliberate corruption of the native solution.
