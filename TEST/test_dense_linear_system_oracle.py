"""Behavioral checks for native capture and the independent dense oracle."""
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from dense_linear_system_oracle import analyze


def test_native_diagnostic_fixture(tmp_path, monkeypatch):
    """A real Fortran writer, complex phases, and a zero RHS survive export."""
    executable = __import__("os").environ.get("NEO2_DIAGNOSTIC_TEST_EXE")
    if not executable:
        __import__("pytest").skip("Set NEO2_DIAGNOSTIC_TEST_EXE to native CTest fixture")
    monkeypatch.setenv("NEO2_DUMP_LINEAR_SYSTEM", str(tmp_path))
    subprocess.run([executable], check=True)
    prefix = tmp_path / "axi_species0"
    result = analyze(prefix)
    assert max(c["equation_backward_error"] for c in result["channels"]) < 1e-15
    assert max(c["relative_solution_error"] for c in result["channels"]) < 1e-14
    assert np.array_equal(result["channels"][1]["native_flux"], [0.0, 0.0, 0.0])
    decomposition = analyze(tmp_path / "na_species0")["source_decomposition"]
    assert np.max(decomposition["closure_absolute_error"]) < 1e-15
    assert np.max(abs(np.asarray(decomposition["f01_flux"]))) > 0.1
    # The oracle must detect an inaccurate native solve, including its flux.
    Path(f"{prefix}_source1_solution.mtx").write_text(
        "%%MatrixMarket matrix array complex general\n2 1\n0 0\n0 0\n")
    corrupt = analyze(prefix)["channels"][0]
    assert corrupt["equation_backward_error"] == 1.0
    assert corrupt["relative_solution_error"] == 1.0
    assert max(corrupt["absolute_flux_error"]) > 0.1
