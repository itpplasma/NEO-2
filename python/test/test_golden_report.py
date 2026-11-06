"""ci/golden_report.py: relative differences of a golden case for the CI artifact."""
import os
import subprocess
import sys

import h5py
import numpy as np

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "..", "ci", "golden_report.py")


def write(path, values):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with h5py.File(path, "w") as f:
        for name, value in values.items():
            f[name] = value


def test_relative_differences_and_plot(tmp_path):
    case = tmp_path / "ql"
    write(str(case / "current" / "summary.h5"),
          {"out/D11": np.array([1.0, 2.0]), "out/only_current": np.array([3.0])})
    write(str(case / "reference" / "summary.h5"),
          {"out/D11": np.array([1.0, 2.002])})
    write(str(case / "legacy" / "summary.h5"), {"out/D11": np.array([1.0, 2.0])})
    env = dict(os.environ, TOL="1e-10", LEGACY_TOL="1e-8", MPLBACKEND="Agg")
    subprocess.run([sys.executable, SCRIPT, str(case), "base-ql", str(tmp_path / "out")],
                   check=True, env=env)
    rows = [line.split("\t") for line in
            (tmp_path / "out" / "base-ql.tsv").read_text().splitlines()[1:]]
    values = {(r[0], r[1]): float(r[2]) for r in rows}
    assert abs(values[("reference", "out/D11")] / (0.002 / 2.002) - 1) < 1e-3
    assert values[("legacy", "out/D11")] == 0.0
    assert ("reference", "out/only_current") not in values
    assert (tmp_path / "out" / "base-ql.png").stat().st_size > 0


def test_missing_current_writes_nothing(tmp_path):
    subprocess.run([sys.executable, SCRIPT, str(tmp_path / "none"), "x",
                    str(tmp_path / "out")], check=True)
    assert not (tmp_path / "out").exists()
