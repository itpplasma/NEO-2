# Optional native work and orbit exports

The ArnoldiOrder2 solver can export the complex quantities consumed by its
non-axisymmetric work contraction. It is disabled when no request file exists.
No transport equation, normalization, resolution or convergence setting changes.

In a fresh run directory, create `NEO2_WORK.REQUEST` containing the selected
normalized flux surface and zero-based response species, for example:

```text
0.09 0
```

A negative species selects every response species on that surface. The surface
selection tolerance is `1e-10`. Axisymmetric-only runs do not open work exports.
Files use a species and serial suffix; existing outputs are never overwritten.
The exporter writes in the solver's process directory. Concurrent processes
must not export the same selection into a shared directory.

`NEO2_WORK_s<species>_j<serial>.txt` records the native vector length, source
channels, geometry, units, model switches, energy basis, pitch grid and Bloch
boundary factors. Its vector files contain:

- `RHS_k1_b*` and `RHS_k3_b*`: original particle/temperature-driving right sides,
  captured before the native iterator overwrites them;
- `FLUX_k1_b*`: the complex weighted particle-flux dual;
- `SOL_k1_b*` and `SOL_k3_b*`: full collision-coupled solutions before the native
  work contraction and deallocation.

For RHS/SOL the `b` suffix denotes the driving-species index; for FLUX it is
the response species, also selected by the prefix. `QFLUX` records the
corresponding native real contractions.
The basis stream records each four-by-four pitch reconstruction matrix in
Fortran column-major order and its actual stencil origin. Energy metadata
references the run's native `collop.h5` basis archive. Unused positive-lower and
negative-upper node delta slots are literal zero; active solver values are
preserved, including nonfinite values for diagnosis.

Binary streams carry magic/version, native byte order, floating-point widths
and record count. Complex vectors use real64/complex128; writes are chunked and
accept strided native rows. `ntv_work_reader.py` checks exact payload lengths,
byte order and layout; its iterators reject nonfinite active vector/basis data.
The text session ends with `END_NEO2_COMPLEX_WORK_V1` only after all expected
basis records and selected vectors have been written. Keep raw streams and the
native energy archive together.

An independent `NEO2_ORBIT.REQUEST` with the same two fields exports the live
magnetic/electric frequency coefficients and phase metadata before they are
freed. This export uses an exact species match; it does not accept the work
export's negative-species wildcard. It reports coefficients, not reconstructed
bounce frequencies or individual harmonic torques.

Registered tests `ntv_work_export_writer` and `ntv_work_export_reader` exercise
disabled operation, chunk boundaries, actual strided rows, preserved right
sides, independent closed-form complex work, pitch layout, inactive NaN slots,
and rejection of truncated/nonfinite data. The reader uses Python's standard
library only. The exports do not define a per-harmonic torque decomposition;
that requires separate physical orbit, dual, phase and convergence checks.
