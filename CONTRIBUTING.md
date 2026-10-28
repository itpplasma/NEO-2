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
