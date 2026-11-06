#!/usr/bin/env python3
"""Relative differences of a golden case, for the CI artifact.

Usage: golden_report.py CASE_DIR NAME OUT_DIR

Compares CASE_DIR/current/summary.h5 with CASE_DIR/reference/summary.h5 and
CASE_DIR/legacy/summary.h5 where present. For every dataset present in both
files it writes the maximum relative difference max|a - b| / max|b| to
OUT_DIR/NAME.tsv and plots the largest ones against the tolerances TOL and
LEGACY_TOL from the environment to OUT_DIR/NAME.png. Only relative
differences leave the run: the golden inputs are not public.
"""
import os
import sys

import h5py
import numpy as np


def datasets(path):
    found = {}

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset) and obj.dtype.kind in "fciu":
            found[name] = np.asarray(obj[()], dtype=complex)

    with h5py.File(path, "r") as f:
        f.visititems(visit)
    return found


def relative_differences(current, other):
    rows = []
    for name, a in sorted(current.items()):
        b = other.get(name)
        if b is None or a.shape != b.shape:
            continue
        scale = np.max(np.abs(b)) if b.size else 0.0
        diff = np.max(np.abs(a - b)) if a.size else 0.0
        rel = diff / scale if scale > 0 else diff
        rows.append((name, float(rel)))
    return rows


def main():
    case_dir, name, out_dir = sys.argv[1:4]
    current_file = os.path.join(case_dir, "current", "summary.h5")
    if not os.path.isfile(current_file):
        return
    os.makedirs(out_dir, exist_ok=True)
    current = datasets(current_file)
    results = {}
    for label, tol_name in (("reference", "TOL"), ("legacy", "LEGACY_TOL")):
        path = os.path.join(case_dir, label, "summary.h5")
        if os.path.isfile(path):
            results[label] = (relative_differences(current, datasets(path)),
                              os.environ.get(tol_name))
    with open(os.path.join(out_dir, name + ".tsv"), "w") as out:
        out.write("against\tdataset\tmax_relative_difference\n")
        for label, (rows, _) in results.items():
            for dataset, rel in rows:
                out.write(f"{label}\t{dataset}\t{rel:.3e}\n")
    plot(results, name, os.path.join(out_dir, name + ".png"))


def plot(results, name, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panels = [(label, rows, tol) for label, (rows, tol) in results.items() if rows]
    if not panels:
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(7 * len(panels), 8),
                             squeeze=False)
    for ax, (label, rows, tol) in zip(axes[0], panels):
        rows = sorted(rows, key=lambda r: r[1], reverse=True)[:40]
        values = [max(r[1], 1e-18) for r in rows]
        ax.barh([r[0] for r in rows][::-1], values[::-1], color="tab:blue")
        ax.set_xscale("log")
        ax.set_xlabel("max |current - %s| / max |%s|" % (label, label))
        if tol:
            ax.axvline(float(tol), color="tab:red", label=f"tolerance {tol}")
            ax.legend(loc="lower right")
        ax.set_title(f"{name}: current vs {label}")
        ax.tick_params(axis="y", labelsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
