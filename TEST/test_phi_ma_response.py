"""Exercise the production local O2 source, solve and current quadrature."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import h5py
import numpy as np


def run(executable, fixture, root, name, potential=None, delta=None, conjugate=False):
    directory = root / name
    directory.mkdir()
    for source in fixture.iterdir():
        if source.suffix in (".in", ".bc"):
            shutil.copy2(source, directory / source.name)
    additions = ""
    for key, value in (("phi_m", potential), ("delta_m", delta)):
        if value is not None:
            additions += f" {key}=({value.real:.17e},{value.imag:.17e}),\n"
    if potential is not None or delta is not None:
        additions += f" phi_ma_m={2 if conjugate else -2},\n"
    config = directory / "neo2.in"
    config.write_text(config.read_text().replace("&ntv_input\n", "&ntv_input\n" + additions))
    if conjugate:
        source = directory / "pert.bc"
        source.write_text(source.read_text().replace("0 1 ", "0 -1 ").replace("1 1 ", "1 -1 "))
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MPLBACKEND="Agg")
    env.pop("DISPLAY", None)
    with (directory / "run.log").open("w") as output:
        result = subprocess.run([str(executable)], cwd=directory, env=env,
                                stdout=output, stderr=subprocess.STDOUT, timeout=45)
    assert result.returncode == 0, directory / "run.log"
    assert "solving completed" in (directory / "run.log").read_text(), directory
    return directory


def read(directory):
    with h5py.File(directory / "phi_ma_current.h5") as source:
        data = {key: source[key][()] for key in source}
    for value in data.values():
        assert np.all(np.isfinite(value))
    return data


def complex_data(data, key):
    return data[key + "_real"] + 1j * data[key + "_imag"]


def identical_legacy(left, right, unordered_text=()):
    count = 0
    for name in ("collop.h5", "fulltransp.h5", "neo2_multispecies_out.h5", "neo2_config.h5"):
        with h5py.File(left / name) as a, h5py.File(right / name) as b:
            assert set(a.keys()) == set(b.keys()), name
            def compare(key, item):
                nonlocal count
                if isinstance(item, h5py.Dataset):
                    # Revision identities and run metadata differ between builds.
                    if name == "neo2_config.h5" and key in (
                            "metadata/Run_0", "metadata/Timestamp_start_0",
                            "metadata/NEO-2 Version", "metadata/NEO-2 Version Additional",
                            "metadata/MPILib Version"):
                        return
                    other = b[key]
                    assert item.shape == other.shape and item.dtype == other.dtype, key
                    assert item[()].tobytes() == other[()].tobytes(), (name, key)
                    count += 1
            a.visititems(compare)
    for source in left.glob("*.dat"):
        other = right / source.name
        if source.name in unordered_text:
            # MPI ranks append this diagnostic without a shared ordering.
            assert sorted(source.read_bytes().splitlines()) == sorted(other.read_bytes().splitlines()), source.name
        else:
            assert source.read_bytes() == other.read_bytes(), source.name
    print(f"PASS: {count} deterministic legacy datasets/text bitwise identical; MPI diagnostics checked as line multisets")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("fixture", type=Path)
    parser.add_argument("work", type=Path)
    parser.add_argument("--switch-off", action="store_true")
    args = parser.parse_args()
    baseline = os.environ.get("NEO2_PHI_MA_BASELINE")
    if args.switch_off and not baseline:
        print("NEO2_PHI_MA_BASELINE must name an origin/main binary")
        raise SystemExit(77)
    root = Path(tempfile.mkdtemp(prefix="phi-ma-", dir=args.work))
    print(f"Run evidence: {root}", flush=True)
    off = run(args.executable, args.fixture, root, "off")
    assert not (off / "phi_ma_current.h5").exists()
    if args.switch_off:
        reference = run(Path(baseline), args.fixture, root, "main")
        identical_legacy(reference, off)
        return
    phi = 1.0e-6 + 4.0e-7j
    positive_path = run(args.executable, args.fixture, root, "positive", phi)
    positive = read(positive_path)
    identical_legacy(off, positive_path)
    j = complex_data(positive, "j_par")
    assert np.max(np.abs(j)) > 0
    negative = read(run(args.executable, args.fixture, root, "negative", -phi))
    np.testing.assert_allclose(complex_data(negative, "j_par"), -j, rtol=2e-8, atol=1e-9)
    print("PASS: drive reversal")
    conjugate = read(run(args.executable, args.fixture, root, "conjugate",
                         phi.conjugate(), conjugate=True))
    np.testing.assert_allclose(complex_data(conjugate, "j_par"), j.conjugate(),
                               rtol=2e-8, atol=1e-9)
    print("PASS: conjugate harmonic")
    aligned = read(run(args.executable, args.fixture, root, "aligned",
                       complex(-positive["Er"], 0), 1.0 + 0j))
    assert complex_data(aligned, "phi_ma") == 0
    assert np.all(complex_data(aligned, "j_par") == 0)
    print("PASS: aligned-potential null")
    weights = positive["surface_weight"]
    baverage = np.sum(weights * positive["bhat"] * j) / np.sum(weights)
    np.testing.assert_allclose(baverage, complex_data(positive, "current_row"), rtol=1e-12)
    average = np.sum(weights * j) / np.sum(weights)
    np.testing.assert_allclose(average, complex_data(positive, "j_par_line_average"), rtol=1e-12)
    assert complex_data(positive, "j_par_fsa") == 0
    print("PASS: local-to-current-row closure and line/full-surface averages")


if __name__ == "__main__":
    main()
