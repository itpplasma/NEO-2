# Owned CI controller handoff — 21 October 2026

The source candidate is [PR 236](https://github.com/itpplasma/NEO-2/pull/236),
6a1a266e550c2b4f9ca3fd6912f7294d653971f4, branch
`ci/owned-capacity-20261021`. This documentation branch does not change that
candidate or its running qualification. The controller owns promotion.

All PR events, including labelled PRs, skip required jobs before allocating
runners. Merge groups and direct default-branch nonqueue pushes validate on
owned capacity; queue-bot/tag publication remains enabled. Required contexts
remain `repo-hygiene / repo-hygiene`, `unit-tests`, `run-golden-record`.
The unconditional hygiene caller pins both shared fields to
104634c01aa669d542e8a26a34a50b6471fc19ce.

Native golden preparation uses the trusted owned `itp-dispatch` runner. The
Slurm native tier reserves eight CPUs, 48 GiB and one hour. Existing native
comparators, numerical tolerances, stable-before-base order, performance and
PAR coverage are preserved. The golden-only toolchain restores the original
initial compiler flags; unit flags are separate. Two MPI ranks remain within
the eight reserved CPUs. Private golden inputs stay in temporary runtime
bundles outside Git and public artifacts; retain only allowed receipts and
scalar summaries. Never publish the private bundle.

At handoff:

- Native GNU 12 unit job 3321711 passed all seven CTest tests (55 s, 422 MiB).
  Its unit/scientific inputs are unchanged in the candidate; do not rerun it
  without a changed input or concrete failure.
- MPI topology probe 3321759 passed both native launchers.
- Exact-candidate golden job 3322336 is running. Stable/base Lorentz and QL
  comparisons and QL performance passed; PAR remains in progress. This is not
  complete golden qualification.
- Previous job 3321775 remains running. Its QL comparisons failed legacy
  `CMake_Flags` metadata, while current/reference numerical comparisons and
  Lorentz/performance passed. Preserve its terminal outcome; do not replace it
  with a success claim. Earlier prebuild/Git-object/MPI startup failures were
  fixed in the harness without weakening native checks.
- Candidate PR workflows 39643381882, 39643382556 and 39643381713 skipped.
  No authoritative merge-group run has occurred yet.

The existing bounded local receipt watcher is collecting the running Slurm
job without promotion. The controller must inspect the job ID and revision
in each receipt: concurrent historical watchers can overwrite the convenience
latest-state file. Per-job terminal artifacts remain distinct.

Resume priorities: collect 3322336; require every applicable native golden
row and scheduler exit to pass; verify the exact candidate and body; then use
the normal squash merge queue under the controller's conditional approval.
Do not bypass the native queue or launch replacement runs from this handoff. Source export
[PR 231](https://github.com/itpplasma/NEO-2/pull/231) has a separate owner and
may queue only after this CI migration lands. Its reviewed head is
277feaf5d77946c17a9c3b115a34b7e1313f211d, ready with no autoqueue. Preserve healthy native jobs.

Verification already completed: actionlint with owned-label configuration,
Bash syntax, warning-level ShellCheck, configured hooks, self-contained bundle
Git connectivity/LFS checks, and synthetic aggregation of all eight possible
comparison commands. Scientific source and thresholds are unchanged.
