# Contributing to NEO-2

## Tests and coverage

Build and run the test suite from the repository root:

    make test

To generate a code coverage report:

    make coverage

This rebuilds with coverage instrumentation, runs all tests and prints a
summary. Coverage data are written to the `build` directory. For a detailed
report:

    lcov --summary build/coverage_filtered.info
    cd build && genhtml coverage_filtered.info --output-directory coverage_html
    # then open build/coverage_html/index.html in a browser

## Merging

- Pull requests, including every layer of a stacked PR, run the fast CI tier.
  Label a PR `full-ci` to run the full tier on it before queueing.
- Review a stack layer by layer. When it is ready, collapse it into one PR into
  `main` (retarget the top layer to `main`); the lower PRs then show as merged.
- Click "Merge when ready" (or `gh pr merge --merge --auto`). The PR enters the
  merge queue, the full CI tier runs on the exact merge result, and on success
  it lands on `main` as a merge commit. `main` itself is not re-tested on push.
- A branch need not be up to date with `main`. Merge `main` into a PR only
  when it conflicts.

## Pull requests that change results on purpose

The merge queue compares the PR build with the merge-queue base (QL and PAR
tolerance 1e-15), with the latest release (QL 1e-10, Lorentz 1e-14) and with
the stored records in the plasma/data repository (`TESTS/NEO-2/golden_record`).
How a PR that changes results is handled depends on how far they move:

- Bit-identical results (refactoring, removed code): nothing to do.
- Changes within the tolerances of the latest release and of the stored
  records, such as a different factorization ordering that moves results by
  rounding: label the PR `golden-change`, put the measured differences in the
  PR body and add a `CHANGELOG` entry. The base and stable comparisons are
  then report-only, and the stored records of plasma/data `main` stay strict
  and unchanged.
- Larger, intended changes: open a plasma/data branch with the same name as
  the PR's head branch, regenerate its records on the cluster from the PR head
  and record the NEO-2 commit and the cluster job in the data commit message.
  Label the PR `golden-change`, put the measured differences against the old
  records in the PR body, merge the data merge request, then queue the PR, and
  cut a release soon after, so that the stable comparison is green again for
  other PRs. In a stack, each layer that changes results brings its own data
  branch and records.
- Every `golden-change` queue run measures run time, peak memory per MPI rank
  and the relative change of every output against the base; once the PR
  lands, this record joins the history site on the `gh-pages` branch.

Without the label every comparison stays strict, and the stored records come
from the data branch named like the PR's head branch if it exists, otherwise
from data `main`.

## Golden-record history

One command builds a NEO-2 revision with a libneo revision (by default its
pin, or before the pin the libneo `main` of the commit date), runs the
golden-record cases and writes one record:

    python -m neo2_util.golden_replay run --neo2 REV --work WORK \
        --data DATA/TESTS/NEO-2/golden_record --out OUT \
        [--paired-build MAIN_TREE --concurrent]

`OUT/record.json` holds the run times, peak memory per MPI rank and, after
`golden_replay compare OUT OTHER --label previous|main`, the relative changes
of every compared output. It contains no absolute outputs and may be
published. `OUT/<case>/summary.h5` holds the outputs and stays with the
golden-record data. `python -m neo2_util.golden_changes events` lists the
record refreshes and tolerance changes, and `golden_history render` writes
the site.
