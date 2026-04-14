import json
import os
from pathlib import Path

import numpy as np

from neo2_ql import compute_transport_k

FIXTURE_REL = "TESTS/NEO-2/golden_record/transport_k/aug_30835_sauter_reference.json"


def _fixture_path() -> Path | None:
    data_dir = os.environ.get("NEO2_DATA_DIR")
    if data_dir is None:
        return None
    p = Path(data_dir) / FIXTURE_REL
    return p if p.is_file() else None


def test_compute_transport_k_formula():
    d31 = np.array([2.0, 4.0])
    d32 = np.array([3.0, 6.0])
    assert np.allclose(compute_transport_k(d31, d32), np.array([1.0, 1.0]))


def test_transport_k_vs_sauter():
    path = _fixture_path()
    if path is None:
        print("SKIP: NEO2_DATA_DIR not set or fixture not found")
        return
    fixture = json.loads(path.read_text(encoding="utf-8"))
    tol = fixture["tolerance"]
    for surface in fixture["surfaces"]:
        k = float(compute_transport_k(surface["d31_ii"], surface["d32_ii"]))
        k_sauter = surface["k_sauter"]
        diff = abs(k - k_sauter)
        assert diff < tol, (
            f"rho_pol={surface['rho_pol']:.4f}: "
            f"|k_neo2 - k_sauter| = {diff:.4e} >= {tol}"
        )


if __name__ == "__main__":
    test_compute_transport_k_formula()
    test_transport_k_vs_sauter()
    print("All transport-k tests passed.")
