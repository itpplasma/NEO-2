#!/bin/bash
# Bundle this revision with its comparison revisions and run the golden-record
# comparisons on sCluster: ci/golden_record.sh BASE_REF
# BASE_REF is the revision for exact and performance comparisons; the stable
# reference is the latest release. Receipts land in $SCI_OUT.
set -euo pipefail
base_ref=${1:?usage: ci/golden_record.sh BASE_REF}
base=$(git rev-parse --verify --end-of-options "$base_ref^{commit}" ||
  git rev-parse --verify --end-of-options "refs/remotes/origin/$base_ref^{commit}")
stable=$(curl -fsS ${GH_TOKEN:+-H "Authorization: Bearer $GH_TOKEN"} \
  https://api.github.com/repos/itpplasma/NEO-2/releases/latest | jq -er .tag_name)
test -n "$stable" && test "$stable" != null
git rev-parse --verify "$stable^{commit}"
export SCI_BASE_SHA=$base SCI_STABLE_REF=$stable SCI_TIER=golden-record
export SCI_OUT=${SCI_OUT:-${RUNNER_TEMP:-/tmp}/golden-receipts}
bundle=${RUNNER_TEMP:-/tmp}/neo2-golden.tar.zst
trap 'rm -f "$bundle"' EXIT
TMPDIR=${RUNNER_TEMP:-/tmp} ci/bundle.sh "$bundle"
sci-dispatch golden "$bundle"
