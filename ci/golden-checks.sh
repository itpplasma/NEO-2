# Golden check helpers, sourced by ci/golden.sh with test_dir and src set.
# shellcheck shell=bash
status=0
check() {
    local name=$1
    shift
    local rc=0
    (cd "$test_dir" && "$@") || rc=$?
    printf '%s\t%s\n' "$name" "$rc" | tee -a "$src/ci-out/checks.tsv"
    [ "$rc" -eq 0 ] || status=1
}
# A PR labelled golden-change alters results on purpose: its differences to
# the stable and base builds are reported, not enforced.
report() {
    local name=$1
    shift
    local rc=0
    (cd "$test_dir" && "$@") 2>&1 | tee "$src/ci-out/$name.log" || rc=$?
    printf '%s\t%s\treport-only\n' "$name" "$rc" | tee -a "$src/ci-out/checks.tsv"
    {
        printf '%s (exit %s)\n' "$name" "$rc"
        grep -E '^(python compare_h5|Difference in |Key .* only found|Warning: number of group keys|ERROR: hdf5)' \
            "$src/ci-out/$name.log" || true
    } >> "$src/ci-out/golden-change.txt"
    rm "$src/ci-out/$name.log"
}
# compare NAME CASE VAR=VALUE...: the stored legacy records stay strict.
compare() {
    local name=$1 case=$2
    shift 2
    if [ "${SCI_GOLDEN_CHANGE:-0}" = 1 ]; then
        report "$name" env "$@" make -C "$case" -k comparison
        check "$name-legacy" env "$@" make -C "$case" legacy_comparison
    else
        check "$name" env "$@" make -k "$case"
    fi
}
