#!/usr/bin/env bash
set -euo pipefail

# Regenerates the NEO-2-QL golden-record (ql) output used by
# python/test/test_force_balance.py. The force-balance replay additionally
# needs dn_spec_ov_ds, dT_spec_ov_ds, Vphi, species_tag_Vphi and isw_Vphi_loc
# in neo2_multispecies_out.h5, which write_multispec_output_a on main does not
# write yet (issue #75). The committed fixture was produced by commit e827e1c
# on branch omte-force-balance-level0 (PR #76), whose NEO-2-QL/ntv_mod.f90
# adds these h5_add calls, from the golden-record ql input deck below.
# load_neo2_force_balance_inputs raises KeyError if any dataset is missing.

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
fixture_path="$repo_root/python/test/data/neo2_ql_axisymmetric_multispecies_out.h5"
golden_dir="${NEO2_GOLDEN_QL_DIR:-/home/ert/data/TESTS/NEO-2/golden_record/ql}"
neo2_exe="${NEO2_QL_EXE:-$repo_root/build/NEO-2-QL/neo_2_ql.x}"
mpi_ranks="${NEO2_MPI_RANKS:-2}"
omp_threads="${OMP_NUM_THREADS:-2}"

if [[ ! -f "$golden_dir/reference/neo.in" ]]; then
    echo "Missing golden-record input deck in $golden_dir" >&2
    exit 1
fi

if [[ ! -x "$neo2_exe" ]]; then
    echo "Missing NEO-2-QL executable at $neo2_exe" >&2
    exit 1
fi

tmpdir=$(mktemp -d /tmp/neo2_ql_fixture.XXXXXX)
trap 'rm -rf "$tmpdir"' EXIT

cp "$golden_dir/reference/neo.in" "$tmpdir/neo.in"
cp "$golden_dir/reference/neo2.in" "$tmpdir/neo2.in"
cp "$golden_dir/test_axi.bc" "$tmpdir/test_axi.bc"
cp "$golden_dir/test_pert.bc" "$tmpdir/test_pert.bc"
ln -sf "$neo2_exe" "$tmpdir/neo_2.x"

(
    cd "$tmpdir"
    OMP_NUM_THREADS="$omp_threads" mpiexec \
        -mca orte_tmpdir_base /tmp \
        -x OMP_NUM_THREADS \
        -np "$mpi_ranks" \
        ./neo_2.x > /tmp/neo2_ql_axisymmetric_fixture.log 2>&1
)

# Refuse to replace the fixture with output the replay cannot use.
PYTHONPATH="$repo_root/python/src" python3 -c \
    "import sys; from neo2_ql.force_balance import load_neo2_force_balance_inputs as f; f(sys.argv[1])" \
    "$tmpdir/neo2_multispecies_out.h5"
cp "$tmpdir/neo2_multispecies_out.h5" "$fixture_path"
echo "Wrote $fixture_path"
echo "Run log: /tmp/neo2_ql_axisymmetric_fixture.log"
