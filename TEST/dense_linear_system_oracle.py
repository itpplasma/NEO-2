"""Independent dense solve of an opt-in small native scalar-system capture.

This checks the native solve and flux contraction, not continuum correctness,
positivity in an invented metric, or convergence of the discretization.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def read_array(path: Path) -> np.ndarray:
    with path.open() as stream:
        header = next(stream).strip().split()
        if len(header) != 5 or header[:2] != ["%%MatrixMarket", "matrix"]:
            raise ValueError(f"Invalid Matrix Market header: {path}")
        if header[3:] != ["complex", "general"]:
            raise ValueError(f"Expected a general complex capture: {path}")
        for line in stream:
            if not line.startswith("%"):
                dimensions = list(map(int, line.split()))
                break
        else:
            raise ValueError(f"Missing dimensions: {path}")
        nrow, ncol = dimensions[:2]
        if not 1 <= nrow <= 4096 or not 1 <= ncol <= 4096:
            raise ValueError(f"Capture dimension outside diagnostic bound: {path}")
        records = np.loadtxt(stream, ndmin=2)
        if header[2] == "coordinate":
            if len(dimensions) != 3 or records.shape != (dimensions[2], 4):
                raise ValueError(f"Truncated coordinate capture: {path}")
            indices = records[:, :2].astype(np.int64)
            if not np.array_equal(indices, records[:, :2]):
                raise ValueError(f"Noninteger coordinate index: {path}")
            if np.any(indices < 1) or np.any(indices[:,0] > nrow) or np.any(indices[:,1] > ncol):
                raise ValueError(f"Coordinate index outside matrix: {path}")
            value = np.zeros((nrow, ncol), dtype=np.complex128)
            np.add.at(value, (indices[:,0]-1, indices[:,1]-1), records[:,2]+1j*records[:,3])
        elif header[2] == "array":
            if len(dimensions) != 2 or records.shape != (nrow*ncol, 2):
                raise ValueError(f"Truncated array capture: {path}")
            value = (records[:,0]+1j*records[:,1]).reshape(ncol,nrow).T
        else:
            raise ValueError(f"Unsupported Matrix Market format: {path}")
    if not np.isfinite(value).all():
        raise ValueError(f"Nonfinite capture: {path}")
    return value


def analyze(prefix: Path) -> dict:
    matrix = read_array(Path(f"{prefix}_matrix.mtx"))
    nrow, ncol = matrix.shape
    if nrow != ncol:
        raise ValueError("Oracle requires a square matrix")
    bra = read_array(Path(f"{prefix}_flux_bra.mtx"))
    rhs = np.column_stack([
        read_array(Path(f"{prefix}_source{k}_rhs.mtx")).reshape(nrow)
        for k in range(1, 4)
    ])
    native = np.column_stack([
        read_array(Path(f"{prefix}_source{k}_solution.mtx")).reshape(nrow)
        for k in range(1, 4)
    ])
    if bra.shape != (nrow, 3):
        raise ValueError("Expected three native flux bras")
    direct = np.linalg.solve(matrix, rhs)
    if not np.isfinite(direct).all():
        raise ValueError("Dense oracle returned nonfinite values")
    matrix_norm = np.linalg.norm(matrix, ord=np.inf)
    tiny = np.finfo(float).tiny
    native_flux = 0.5 * (bra.conj().T @ native).real
    direct_flux = 0.5 * (bra.conj().T @ direct).real
    channels = []
    for channel in range(3):
        x, b = native[:, channel], rhs[:, channel]
        residual = matrix @ x - b
        scale = matrix_norm * np.linalg.norm(x, ord=np.inf)
        scale += np.linalg.norm(b, ord=np.inf)
        relative_solution_error = np.linalg.norm(x - direct[:, channel], ord=np.inf)
        relative_solution_error /= max(np.linalg.norm(direct[:, channel], ord=np.inf), tiny)
        channels.append({
            "source": channel + 1,
            "equation_backward_error": float(np.linalg.norm(residual, ord=np.inf) / max(scale, tiny)),
            "relative_solution_error": float(relative_solution_error),
            "native_flux": native_flux[:, channel].tolist(),
            "direct_flux": direct_flux[:, channel].tolist(),
            "absolute_flux_error": abs(native_flux[:, channel] - direct_flux[:, channel]).tolist(),
        })
    result = {"schema": 1, "prefix": str(prefix), "unknowns": nrow, "channels": channels}
    # Native source decomposition is captured independently of the solver.
    if prefix.name.startswith("na_species"):
        species = int(prefix.name.removeprefix("na_species"))
        geometric = read_array(Path(f"{prefix}_drive{species}_geometric_source.mtx"))
        total = read_array(Path(f"{prefix}_drive{species}_total_source.mtx"))
        if not np.allclose(total, rhs, rtol=0.0, atol=0.0):
            raise ValueError("Captured total source differs from the solved RHS")
        pieces = np.concatenate([geometric, total - geometric], axis=1)
        solutions = np.linalg.solve(matrix, pieces)
        fluxes = 0.5 * (bra.conj().T @ solutions).real
        result["source_decomposition"] = {
            "geometric_flux": fluxes[:, :3].tolist(),
            "f01_flux": fluxes[:, 3:].tolist(),
            "closure_absolute_error": abs(fluxes[:, :3] + fluxes[:, 3:] - direct_flux).tolist(),
        }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prefix", type=Path, help="e.g. capture/na_species0")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.prefix)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
