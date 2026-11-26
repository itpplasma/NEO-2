"""Local parallel flow: its flux-surface average must equal D3j and ParFlow."""
import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from test_phi_ma_response import identical_legacy  # noqa: E402

# Synthetic n=1 |B| perturbation (T) and a temperature gradient, so that
# both problems and all three drives carry a parallel flow.
B_PERTURBATION = "1e-3"


def run(executable, fixture, root, name, switch):
    directory = root / name
    directory.mkdir()
    for source in fixture.iterdir():
        if source.suffix in (".in", ".bc"):
            shutil.copy2(source, directory / source.name)
    pert = directory / "pert.bc"
    pert.write_text(re.sub(r"^1 1 0 0 0 0 0 0 0 0$", f"1 1 0 0 0 0 0 0 {B_PERTURBATION} 0",
                           pert.read_text(), flags=re.M))
    config = directory / "neo2.in"
    text = config.read_text().replace("dt_vec_ov_ds=0", "dt_vec_ov_ds=-1d-9")
    if switch:
        text = text.replace("&ntv_input\n", "&ntv_input\n lsw_local_parflow=.true.,\n")
    config.write_text(text)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MPLBACKEND="Agg")
    env.pop("DISPLAY", None)
    with (directory / "run.log").open("w") as output:
        result = subprocess.run([str(executable)], cwd=directory, env=env,
                                stdout=output, stderr=subprocess.STDOUT, timeout=60)
    assert result.returncode == 0, directory / "run.log"
    return directory


def check_closure(directory):
    with h5py.File(directory / "local_parallel_flow.h5") as source:
        local = {key: source[key][()] for key in source}
    with h5py.File(directory / "neo2_multispecies_out.h5") as source:
        legacy = {key: source[key][()] for key in source}
    weight = local["surface_weight"]
    volume = local["surface_volume"]
    print(f"sum(surface_weight)/surface_volume - 1 = {weight.sum() / volume - 1:.3e}")
    # Existing transport rows: AX keeps Re, NA the Fourier half product.
    for problem, factor in (("AX", 1.0), ("NA", 0.5)):
        d3 = local[f"D3j_{problem}_real"] + 1j * local[f"D3j_{problem}_imag"]
        # h5py reverses the Fortran order: (flux species, drive species, j, step).
        nspec = d3.shape[0]
        average = factor * np.real(np.einsum("abjs,s->abj", d3, weight)) / volume
        for j in range(3):
            expected = legacy[f"D3{j + 1}_{problem}"].reshape(nspec, nspec)
            np.testing.assert_allclose(average[:, :, j], expected, rtol=1e-12,
                                       atol=1e-14 * np.max(np.abs(expected)))
        assert np.max(np.abs(average)) > 0, problem
        vparb = local[f"VparB_{problem}_real"] + 1j * local[f"VparB_{problem}_imag"]
        flow = factor * np.real(vparb @ weight) / volume
        # ParFlow = -sum_bj D3j*Aj cancels strongly for ambipolar Er, so the
        # roundoff bound scales with the summands, not with the sum.
        scale = np.abs(average * local["A_forces"][np.newaxis]).sum(axis=(1, 2))
        np.testing.assert_array_less(np.abs(flow - legacy[f"ParFlow_{problem}_spec"]),
                                     1e-13 * scale + 1e-300)
        print(f"PASS: {problem} local D3j and VparB average to D3j_{problem} "
              f"and ParFlow_{problem}_spec = {legacy[f'ParFlow_{problem}_spec']}")
    # The real axisymmetric solve leaves only complex-arithmetic roundoff.
    for key in ("D3j_AX", "VparB_AX"):
        assert np.max(np.abs(local[key + "_imag"])) < 1e-12 * np.max(np.abs(local[key + "_real"]))
    jpar = local["j_par_NA_real"] + 1j * local["j_par_NA_imag"]
    print(f"NA helical j_par: max |j| = {np.max(np.abs(jpar)):.6e} statA/cm^2, "
          f"line average = {np.sum(weight * jpar) / weight.sum():.6e}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("fixture", type=Path)
    parser.add_argument("work", type=Path)
    args = parser.parse_args()
    root = Path(tempfile.mkdtemp(prefix="local-parflow-", dir=args.work))
    print(f"Run evidence: {root}", flush=True)
    off = run(args.executable, args.fixture, root, "off", switch=False)
    assert not (off / "local_parallel_flow.h5").exists()
    on = run(args.executable, args.fixture, root, "on", switch=True)
    identical_legacy(off, on)
    check_closure(on)
    baseline = os.environ.get("NEO2_LOCAL_PARFLOW_BASELINE")
    if baseline:
        identical_legacy(run(Path(baseline), args.fixture, root, "base", switch=False), off)
        print("PASS: switch off is bitwise identical to the baseline binary")
    else:
        print("SKIP: set NEO2_LOCAL_PARFLOW_BASELINE for the baseline comparison")


if __name__ == "__main__":
    main()
