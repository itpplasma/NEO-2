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
A PR that changes results on purpose fails the first two by design. For such a
PR:

- Open a plasma/data branch with the same name as the PR's head branch, and
  regenerate its records on the cluster from the PR head. Record the NEO-2
  commit and the cluster job in the data commit message.
- Put the measured differences against the old records in the PR body.
- Label the PR `golden-change`. The queue then runs the base and stable
  comparisons in report-only mode: their differences appear in the job
  summary, marked as report-only, but do not fail the check. The comparison
  with the stored records stays strict and uses the data branch; the run fails
  if that branch is missing.
- Merge the data merge request, then queue the PR.
- Cut a release soon after, so that the stable comparison is green again for
  other PRs.
- In a stack, each layer that changes results brings its own data branch and
  records.

Without the label every comparison stays strict, and the stored records come
from the data branch named like the PR's head branch if it exists, otherwise
from data `main`.
