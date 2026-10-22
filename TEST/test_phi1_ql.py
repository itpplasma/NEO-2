"""Bounded synthetic two-species Phi_1 closure and superposition test."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import h5py
import numpy as np

executable, fixtures, scratch, mpiexec = sys.argv[1:]
scratch = Path(scratch)
scratch.mkdir(parents=True, exist_ok=True)
env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MPLBACKEND="Agg")
env.pop("DISPLAY", None)
results = {}
for name in ("driven", "zero", "unequal"):
    run = scratch / name
    run.mkdir(exist_ok=True)
    for filename in ("neo.in", "neo2.in", "background.bc", "pert.bc"):
        shutil.copyfile(Path(fixtures) / filename, run / filename)
    if name == "unequal":
        p = run / "neo2.in"
        p.write_text(p.read_text().replace("t_vec=1.6d-9,1.6d-9", "t_vec=1.6d-9,8d-10"))
    if name == "zero":
        # Zero every magnetic perturbation coefficient, leaving the operator,
        # profiles and the independently assembled Phi responses unchanged.
        p = run / "pert.bc"
        lines = p.read_text().splitlines()
        for i, line in enumerate(lines):
            fields = line.split()
            if len(fields) == 10:
                try:
                    int(fields[0]); int(fields[1])
                except ValueError:
                    continue
                fields[8:10] = ["0", "0"]
                lines[i] = " ".join(fields)
        p.write_text("\n".join(lines) + "\n")
    with (run / "stdout.log").open("w") as log:
        subprocess.run([mpiexec, "-np", "2", executable], cwd=run,
                       env=env, stdout=log, stderr=subprocess.STDOUT,
                       check=True, timeout=25)
    with h5py.File(run / "neo2_multispecies_out.h5") as f:
        g = f["phi1"]
        for key in g:
            assert np.isfinite(g[key][()]).all(), (name, key)
        potential = g["Phi_re"][()] + 1j*g["Phi_im"][()]
        residual = float(g["residual"][()])
        assert residual < 1e-10, residual
        collision = np.max(g["collision_source_norm"][()])
        if name == "unequal":
            assert collision > 1e-12, "missing unequal-T collision source"
        else:
            assert collision == 0, "common-T S_coll"
        terms = 4.8032e-10*np.asarray(g["Gamma_ExB"][()])
        remainder = float(g["chi_Phi1_Ur"][()])
        check = float(g["qn_flux_check"][()])
        scale = abs(terms[0]-terms[1]) + abs(remainder)
        relative = abs(check)/max(scale, 1e-100)
        assert relative < 1e-2, ("D9", relative, check, scale)
        if name == "zero":
            assert np.max(np.abs(potential)) < 1e-20, potential
            assert np.max(np.abs(g["Gamma_total"][()])) < 1e-20
        else:
            assert np.max(np.abs(potential)) > 1e-12, "missing magnetic response"
        results[name] = {"Phi_re": potential.real.tolist(), "Phi_im": potential.imag.tolist(),
                         "residual": residual, "D9_relative": relative,
                         "collision_source_norm": float(collision),
                         "Gamma_total": g["Gamma_total"][()].tolist()}
print(json.dumps(results, indent=2))
print("All tests passed!")
