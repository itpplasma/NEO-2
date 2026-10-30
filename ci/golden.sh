#!/bin/bash
# Native golden authority, run on the bounded Slurm reference CPU tier.
set -euo pipefail
src=$(pwd)
root=${SCI_ROOT:?Slurm runtime root missing}
export MPI_HOME=/usr BLA_VENDOR=OpenBLAS
data=$root/golden-data
test_dir=$data/TESTS/NEO-2/golden_record
python3 -m venv "$root/venv"
. "$root/venv/bin/activate"
if [ "${SCI_CPUS:-0}" -lt 8 ]; then
    echo 'Golden comparisons require eight allocated CPUs' >&2
    exit 2
fi
# The profile reserves eight CPUs as one Slurm task. Existing cases launch
# two MPI ranks with up to four OpenMP threads, confined by the Slurm cpuset.
printf 'localhost slots=2 max_slots=2\n' > "$root/golden-mpi-hosts"
export OMPI_MCA_ras='^slurm'
export OMPI_MCA_orte_default_hostfile=$root/golden-mpi-hosts
export OMPI_MCA_hwloc_base_binding_policy=none
printf 'ranks=2\nallocated_cpus=%s\n' "$SCI_CPUS" > ci-out/mpi-topology.txt
grep Cpus_allowed_list /proc/self/status >> ci-out/mpi-topology.txt
python -m pip install --upgrade pip
libneo_ref=$(python -c 'import json; print(json.load(open("cmake/dependency_versions.json"))["libneo"])')
test -n "$libneo_ref" && test "$libneo_ref" != null
python -m pip install "git+https://github.com/itpplasma/libneo.git@$libneo_ref"
python -m pip install -e python
. "$data/activate.sh"
(cd "$data" && sha256sum -c inputs.sha256 > /dev/null)
cp "$data/revision.txt" ci-out/data-revision.txt
cp "$data/inputs.sha256" ci-out/inputs.sha256
printf 'source=%s\nbase=%s\nstable=%s\n' "$(git rev-parse HEAD)" \
    "${SCI_BASE_SHA:?base missing}" "${SCI_STABLE_REF:?stable missing}" > ci-out/revisions.txt

make
export CURRENT_NEO2_PAR=$src/build/NEO-2-PAR/neo_2_par.x
export CURRENT_NEO2_QL=$src/build/NEO-2-QL/neo_2_ql.x
for reference in stable base; do
    git clone --no-hardlinks --dissociate "$src" "$root/$reference"
done
git -C "$root/stable" checkout --detach "$SCI_STABLE_REF"
make -C "$root/stable" LIBNEO_REF="$SCI_STABLE_REF"
status=0
check() {
    local name=$1
    shift
    local rc=0
    (cd "$test_dir" && "$@") || rc=$?
    printf '%s\t%s\n' "$name" "$rc" | tee -a "$src/ci-out/checks.tsv"
    [ "$rc" -eq 0 ] || status=1
}
export REFERENCE_NEO2_PAR=$root/stable/build/NEO-2-PAR/neo_2_par.x
export REFERENCE_NEO2_QL=$root/stable/build/NEO-2-QL/neo_2_ql.x
check stable-lorentz env TOL=1e-14 LEGACY_TOL=1e-10 MAGNETICS_TOL=1e-12 make -k lorentz
check stable-ql env TOL=1e-10 LEGACY_TOL=1e-8 MAGNETICS_TOL=1e-12 make -k ql

# Keep the native order: a freshly built base reference invalidates the
# timestamp-based case outputs from the stable-reference comparisons.
git -C "$root/base" checkout --detach "$SCI_BASE_SHA"
make -C "$root/base"
export REFERENCE_NEO2_PAR=$root/base/build/NEO-2-PAR/neo_2_par.x
export REFERENCE_NEO2_QL=$root/base/build/NEO-2-QL/neo_2_ql.x
check base-lorentz env TOL=1e-14 LEGACY_TOL=1e-10 MAGNETICS_TOL=1e-12 make -k lorentz
check base-ql env TOL=1e-15 LEGACY_TOL=1e-8 MAGNETICS_TOL=1e-12 make -k ql
check performance-ql env TOL=1e-15 make performance_ql
check par env TOL=1e-15 LEGACY_TOL=1e-8 make par
if [ -f "$test_dir/performance_ql/performance_report.txt" ]; then
    cp "$test_dir/performance_ql/performance_report.txt" ci-out/
fi

old_ref=$(git -C "$root/base/build/_deps/libneo-src" rev-parse HEAD)
new_ref=$(git -C build/_deps/libneo-src rev-parse HEAD)
printf 'libneo_base=%s\nlibneo_current=%s\n' "$old_ref" "$new_ref" >> ci-out/revisions.txt
if [ "$old_ref" != "$new_ref" ]; then
    mkdir "$root/previous-libneo"
    git archive HEAD | tar -x -C "$root/previous-libneo"
    make -C "$root/previous-libneo" LIBNEO_REF="$old_ref"
    export REFERENCE_NEO2_PAR=$root/previous-libneo/build/NEO-2-PAR/neo_2_par.x
    export REFERENCE_NEO2_QL=$root/previous-libneo/build/NEO-2-QL/neo_2_ql.x
    check libneo-lorentz env TOL=1e-14 MAGNETICS_TOL=1e-12 make -k lorentz
    check libneo-ql env TOL=1e-15 make -k ql
fi
# ci-out contains receipts and scalar performance summaries, no private raw data.
exit "$status"
