"""Synthetic response superposition, D10 channel and default-off isolation."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import h5py
import numpy as np
from test_phi_ma_response import identical_legacy

exe, fixture, profiles, work, mpi = sys.argv[1:]
root = Path(tempfile.mkdtemp(prefix="phi1-response-", dir=work))
print(f"Run evidence: {root}", flush=True)
env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MPLBACKEND="Agg")
env.pop("DISPLAY", None)

def run(name, modes=None, unequal=False, executable=exe):
    target = root / name
    target.mkdir()
    for source in Path(fixture).glob("*"):
        if source.suffix in (".in", ".bc"):
            shutil.copyfile(source, target/source.name)
    # Reflect L1's synthetic embedding to right-handed (s,theta,phi).
    source = target/"background.bc"
    lines = source.read_text().splitlines()
    for i,line in enumerate(lines):
        fields = line.split()
        if len(fields) == 10 and fields[0] == "1":
            fields[5] = str(-float(fields[5]))
            lines[i] = " ".join(fields)
    source.write_text("\n".join(lines)+"\n")
    config = Path(profiles).read_text()
    if modes is None:
        config = config.replace(" isw_phi1=1, phi1_num_m=3, phi1_m=-1,0,1,", "")
    else:
        config = config.replace("phi1_num_m=3, phi1_m=-1,0,1",
                                f"phi1_num_m={len(modes)}, phi1_m=" + ",".join(map(str,modes)))
    if unequal:
        config = config.replace("t_vec=1.6d-9,1.6d-9", "t_vec=1.6d-9,8d-10")
    (target/"neo2.in").write_text(config)
    with (target/"run.log").open("w") as log:
        subprocess.run([mpi, "-np", "2", str(executable)], cwd=target, env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True, timeout=60)
    return target

def data(path):
    with h5py.File(path/"neo2_multispecies_out.h5") as f:
        g = f["phi1"]
        return {k:g[k][()] for k in g}

def z(d, key):
    return d[key+"_re"] + 1j*d[key+"_im"]

off = run("off")
baseline = os.environ.get("NEO2_PHI1_BASELINE")
if not baseline:
    raise RuntimeError("NEO2_PHI1_BASELINE must name the L1 baseline executable")
identical_legacy(run("baseline", executable=baseline), off,
                 unordered_text=("phi_placement_problem.dat",))
a = data(run("basis", [-1,0,1]))
b = data(run("single", [-1]))
np.testing.assert_allclose(z(a,"density_response")[:,0,0],
                           z(b,"density_response")[:,0,0], rtol=1e-10, atol=1e-9)
assert np.max(np.abs(z(a,"density_response"))) > 0
assert np.all(z(a,"density_response")[:,1,:] == 0), "aligned m=0 ExB null"
assert np.all(z(a,"density_dB") == 0), "zero magnetic drive null"
assert np.all(a["collision_source_norm"] == 0), "common-T D10 null"
c = data(run("unequal", [-1], unequal=True))
assert np.max(c["collision_source_norm"]) > 0, "missing D10 source"
for d in (a,b,c):
    assert all(np.isfinite(v).all() for v in d.values())
print("PASS: L1/off bitwise, per-basis superposition, aligned/zero-drive null, D10")
