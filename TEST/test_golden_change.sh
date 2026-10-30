#!/bin/bash
# Dry run of the golden-change path with stubbed comparisons and API.
set -euo pipefail
command -v jq > /dev/null || exit 77
ci=$(cd "$(dirname "${BASH_SOURCE[0]}")/../ci" && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
fail() { echo "FAIL: $*" >&2; exit 1; }

# The API stub serves a PR on branch fix/x, labelled when LABEL is set.
mkdir "$work/bin"
cat > "$work/bin/curl" <<'EOF'
#!/bin/bash
[[ ${!#} == */repos/itpplasma/NEO-2/pulls/255 ]] || exit 22
extra=
[ -z "${LABEL:-}" ] || extra=",{\"name\":\"$LABEL\"}"
echo "{\"head\":{\"ref\":\"fix/x\"},\"labels\":[{\"name\":\"bug\"}$extra]}"
EOF
chmod +x "$work/bin/curl"
resolve() {
    : > "$work/env"
    env PATH="$work/bin:$PATH" GITHUB_ENV="$work/env" GH_TOKEN=t \
        GITHUB_REPOSITORY=itpplasma/NEO-2 GITHUB_HEAD_REF= "$@" "$ci/resolve-pr.sh"
    tr '\n' ' ' < "$work/env"
}
queue=gh-readonly-queue/main/pr-255-707a2698ca85bf14b6c5fdf2d27a39f35019401a
got=$(resolve GITHUB_EVENT_NAME=merge_group GITHUB_REF_NAME=$queue LABEL=golden-change)
[ "$got" = 'SCI_DATA_BRANCH=fix/x SCI_GOLDEN_CHANGE=1 ' ] || fail "labelled queue: $got"
got=$(resolve GITHUB_EVENT_NAME=merge_group GITHUB_REF_NAME=$queue)
[ "$got" = 'SCI_DATA_BRANCH=fix/x SCI_GOLDEN_CHANGE=0 ' ] || fail "unlabelled queue: $got"
got=$(resolve GITHUB_EVENT_NAME=push GITHUB_REF_NAME=main)
[ "$got" = 'SCI_DATA_BRANCH=main SCI_GOLDEN_CHANGE=0 ' ] || fail "push: $got"

# The case Makefiles mirror the golden-record targets; the base comparison
# differs, and the legacy comparison follows LEGACY_RC.
test_dir=$work/golden
mkdir -p "$test_dir/lorentz"
printf '.PHONY: lorentz\nlorentz:\n\tmake -C lorentz\n' > "$test_dir/Makefile"
cat > "$test_dir/lorentz/Makefile" <<'EOF'
all: comparison legacy_comparison
comparison:
	@echo python compare_h5.py reference/summary.h5 current/summary.h5 $(TOL)
	@echo 'Difference in D11: 0.0016'
	@exit 1
legacy_comparison:
	@exit $(LEGACY_RC)
EOF
run() {
    local src=$work/src status
    rm -rf "$src"
    mkdir -p "$src/ci-out"
    export LEGACY_RC=$2 SCI_GOLDEN_CHANGE=$1
    . "$ci/golden-checks.sh"
    compare base-lorentz lorentz TOL=1e-14 > /dev/null 2>&1
    echo "$status"
}
[ "$(run 0 0)" = 1 ] || fail 'unlabelled base difference passed'
[ "$(run 1 0)" = 0 ] || fail 'labelled base difference failed'
grep -qx 'Difference in D11: 0.0016' "$work/src/ci-out/golden-change.txt" ||
    fail 'difference missing from report'
grep -qx $'base-lorentz\t2\treport-only' "$work/src/ci-out/checks.tsv" ||
    fail 'report-only row missing'
[ "$(run 1 2)" = 1 ] || fail 'labelled legacy difference passed'
echo 'golden-change dry run passed'
