"""ci/golden_measure.py: wrapper timing and peak memory of golden-case runs."""
import json
import os
import subprocess
import sys

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "..", "ci", "golden_measure.py")


def test_wrapper_records_case_time_memory_and_exit(tmp_path):
    solver = tmp_path / "fake_solver.py"
    solver.write_text("import sys\nbuf = bytearray(64 * 2**20)\nsys.exit(int(sys.argv[1]))\n")
    binary = tmp_path / "neo_2_ql.x"
    binary.write_text(f'#!/bin/sh\nexec "{sys.executable}" "{solver}" "$@"\n')
    binary.chmod(0o755)
    log = tmp_path / "runs" / "current.ranks"
    log.parent.mkdir()
    wrapper = subprocess.run([sys.executable, SCRIPT, "wrap", str(log), str(binary)],
                             check=True, capture_output=True, text=True).stdout.strip()
    run_dir = tmp_path / "ql" / "current"
    run_dir.mkdir(parents=True)
    assert subprocess.run([wrapper, "0"], cwd=run_dir).returncode == 0
    assert subprocess.run([wrapper, "3"], cwd=run_dir).returncode == 3
    entries = [json.loads(line) for line in log.read_text().splitlines()]
    assert [e["case"] for e in entries] == ["ql", "ql"]
    assert [e["exit"] for e in entries] == [0, 3]
    assert all(e["peak_rss_mib"] >= 64 for e in entries)
    table = subprocess.run([sys.executable, SCRIPT, "summary", str(log)], check=True,
                           capture_output=True, text=True).stdout.splitlines()
    assert table[0].split("\t") == ["case", "build", "runs", "median_wall_s", "peak_rss_mib"]
    case, build, runs, _, peak = table[1].split("\t")
    assert (case, build, runs) == ("ql", "current", "2") and float(peak) >= 64
