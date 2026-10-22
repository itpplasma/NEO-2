"""Coupled synthetic QN, D9 pairing, null and native-reader round trip."""
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import h5py
import numpy as np

exe, reader, fixture, profiles, work, mpi = sys.argv[1:]
root = Path(tempfile.mkdtemp(prefix="phi1-closure-", dir=work))
print(f"Run evidence: {root}", flush=True)
env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MPLBACKEND="Agg")
env.pop("DISPLAY", None)
summary = {}
for name in ("driven", "zero", "unequal"):
    run = root/name
    run.mkdir()
    for source in Path(fixture).glob("*"):
        if source.suffix in (".in", ".bc"):
            shutil.copyfile(source, run/source.name)
    # Same right-handed reflection of L1's synthetic embedding as L2.
    source = run/"background.bc"
    lines = source.read_text().splitlines()
    for i,line in enumerate(lines):
        fields = line.split()
        if len(fields) == 10 and fields[0] == "1":
            fields[5] = str(-float(fields[5]))
            lines[i] = " ".join(fields)
    source.write_text("\n".join(lines)+"\n")
    config = Path(profiles).read_text().replace("isw_phi1=1", "isw_phi1=2")
    if name == "unequal":
        config = config.replace("t_vec=1.6d-9,1.6d-9", "t_vec=1.6d-9,8d-10")
    (run/"neo2.in").write_text(config)
    if name != "zero":
        source = run/"pert.bc"
        lines = source.read_text().splitlines()
        for i,line in enumerate(lines):
            fields = line.split()
            if len(fields) != 10 or not fields[0].lstrip("-").isdigit():
                continue
            fields[-2:] = ["0.0002", "0"] if fields[0] == "0" else ["0.0001", "0.00005"]
            lines[i] = " ".join(fields)
        source.write_text("\n".join(lines)+"\n")
    with (run/"run.log").open("w") as log:
        subprocess.run([mpi,"-np","2",exe],cwd=run,env=env,stdout=log,
                       stderr=subprocess.STDOUT,check=True,timeout=90)
    with h5py.File(run/"neo2_multispecies_out.h5") as f:
        g = {k:v[()] for k,v in f["phi1"].items()}
    assert all(np.isfinite(v).all() for v in g.values())
    phi = g["Phi_re"]+1j*g["Phi_im"]
    matrix = g["M_re"]+1j*g["M_im"]
    rhs = g["rhs_re"]+1j*g["rhs_im"]
    np.testing.assert_allclose(matrix.T@phi,rhs,rtol=1e-11,atol=1e-12)
    assert g["residual"] < 1e-12
    # D9 compares the signed raw ftilde ExB current against its Boltzmann term.
    terms = 4.8032e-10*g["Gamma_ExB"]
    remainder = g["chi_Phi1_Ur"]
    defect = terms[0]-terms[1]-remainder
    scale = abs(terms[0]-terms[1])+abs(remainder)
    assert abs(defect) < 2e-9*max(scale,1e-15), (defect,scale)
    if name == "zero":
        assert np.all(phi == 0) and np.all(g["Gamma_total"] == 0)
    else:
        assert np.max(abs(phi)) > 1e-12
    lines = (run/"phi1.bc").read_text().splitlines()
    assert lines[0] == "CC quantity: Phi_1 [statvolt] in bmnc/bmns columns"
    assert len(lines[5].split()) == 7 and int(lines[5].split()[2]) == 1
    table = np.loadtxt(io.StringIO("\n".join(lines[10:])), ndmin=2)
    assert np.all(table[:,2:8] == 0)
    np.testing.assert_array_equal(table[:,0],g["m"])
    np.testing.assert_array_equal(table[:,1]*int(lines[5].split()[3]),
                                  np.full(len(phi),g["n"]))
    np.testing.assert_allclose(table[:,8]-1j*table[:,9],phi,rtol=2e-16,atol=0)
    # The native perturbation spline requires multiple surfaces. Repeat the
    # emitted spectrum and sample exactly at its original s=.5 knot.
    header = lines[5].split()
    header[2] = "4"
    replicated = lines[:5]+[" ".join(header)]
    for radius in (.25,.5,.75,1.0):
        surface = lines[8].split()
        surface[0] = str(radius)
        replicated += lines[6:8]+[" ".join(surface)]+lines[9:]
    native = run/"phi1-reader.bc"
    native.write_text("\n".join(replicated)+"\n")
    text = subprocess.check_output([reader,str(native)],env=env,text=True)
    values = np.loadtxt(io.StringIO(text),ndmin=2)
    for j, row in enumerate(values,1):
        expected = np.sum(phi*np.exp(1j*(g["m"]*.31*j+g["n"]*.79*j)))
        np.testing.assert_allclose(row[0]+1j*row[1],expected,rtol=2e-13,atol=1e-20)
    assert np.allclose([float(v) for v in lines[8].split()], [.5,.7,3e7,0,0,1])
    theta, zeta = .31, .79
    field = np.sum(table[:,8]*np.cos(table[:,0]*theta+g["n"]*zeta)
                   + table[:,9]*np.sin(table[:,0]*theta+g["n"]*zeta))
    expected = np.real(np.sum(phi*np.exp(1j*(g["m"]*theta+g["n"]*zeta))))
    np.testing.assert_allclose(field,expected,rtol=1e-14,atol=1e-20)
    summary[name] = dict(residual=float(g["residual"]),rcond=float(g["rcond"]),
                         D9_defect=float(defect),Phi_re=phi.real.tolist(),Phi_im=phi.imag.tolist())
(root/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2))
print("PASS: QN residual, D9, zero-Phi null and native reader round trip")
