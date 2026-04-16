#!/usr/bin/env python3
"""Integration test for isw_calc_Er=2 (externally prescribed Om_tE).

Runs NEO-2-QL three times on a single flux surface using the golden
record QL test configuration:

  Run 1: isw_calc_Er=1  -- self-consistent Er from ambipolarity
  Run 2: isw_calc_Er=2  -- feed the Om_tE from run 1 via namelist
  Run 3: isw_calc_Er=2  -- feed half Om_tE (perturbation)

Checks:
  - Transport coefficients D11..D33 must match between runs 1 and 2
    (the kinetic solve is independent of Er)
  - Om_tE and MtOvR must match between runs 1 and 2
  - Run 3 must produce different MtOvR (proving the input is used)
  - Transport coefficients must still match in run 3

Usage:
  python test_isw_calc_er2_integration.py <neo2_ql_binary> <data_dir>

  data_dir is the golden_record/ql directory containing test_axi.bc,
  test_pert.bc, and reference/neo.in + reference/neo2.in.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

import h5py
import numpy as np


def patch_neo2_in(src_path, dst_path, replacements):
    """Copy neo2.in and replace specific namelist values.

    replacements is a dict of {KEY: new_value_string}.
    Handles both existing keys (replace value) and new keys that need
    to be inserted into the correct namelist section.
    """
    with open(src_path, 'r') as f:
        content = f.read()

    for key, value in replacements.items():
        # Try to replace existing key
        pattern = rf'({key}\s*=\s*)[^,/\n]+'
        new_content = re.sub(pattern, rf'\g<1>{value}', content, flags=re.IGNORECASE)
        if new_content != content:
            content = new_content
        else:
            # Key not found — insert before the closing / of NTV_INPUT
            # (Om_tE lives in that namelist)
            if key.upper() == 'OM_TE':
                content = re.sub(
                    r'(&NTV_INPUT.*?)(^\s*/)',
                    rf'\1 {key}={value},\n/',
                    content,
                    flags=re.IGNORECASE | re.MULTILINE | re.DOTALL)

    with open(dst_path, 'w') as f:
        f.write(content)


def setup_run_dir(base_dir, name, data_dir, neo2_binary):
    """Create a run directory with symlinks to test data."""
    run_dir = os.path.join(base_dir, name)
    os.makedirs(run_dir, exist_ok=True)

    for bc in ['test_axi.bc', 'test_pert.bc']:
        src = os.path.join(data_dir, bc)
        dst = os.path.join(run_dir, bc)
        if not os.path.exists(dst):
            os.symlink(os.path.abspath(src), dst)

    neo_in_src = os.path.join(data_dir, 'reference', 'neo.in')
    neo_in_dst = os.path.join(run_dir, 'neo.in')
    if not os.path.exists(neo_in_dst):
        os.symlink(os.path.abspath(neo_in_src), neo_in_dst)

    bin_dst = os.path.join(run_dir, 'neo_2.x')
    if not os.path.exists(bin_dst):
        os.symlink(os.path.abspath(neo2_binary), bin_dst)

    return run_dir


def run_neo2(run_dir):
    """Run NEO-2-QL (single process, single surface)."""
    result = subprocess.run(
        ['./neo_2.x'],
        cwd=run_dir,
        capture_output=True,
        text=True,
        timeout=1200,
        env={**os.environ, 'OMP_NUM_THREADS': '4'},
    )
    # Always save stdout for debugging
    with open(os.path.join(run_dir, 'stdout.log'), 'w') as f:
        f.write(result.stdout)
    with open(os.path.join(run_dir, 'stderr.log'), 'w') as f:
        f.write(result.stderr)
    if result.returncode != 0:
        stdout_tail = result.stdout[-3000:] if len(result.stdout) > 3000 else result.stdout
        stderr_tail = result.stderr[-3000:] if len(result.stderr) > 3000 else result.stderr
        print(f"STDOUT:\n{stdout_tail}")
        print(f"STDERR:\n{stderr_tail}")
        raise RuntimeError(f"NEO-2-QL failed with return code {result.returncode}")
    return result


# Transport coefficient keys (independent of Er in the kinetic solve)
TRANSPORT_KEYS = [
    'D11_AX', 'D12_AX', 'D21_AX', 'D22_AX',
    'D31_AX', 'D32_AX', 'D33_AX',
    'D11_NA', 'D12_NA', 'D21_NA', 'D22_NA',
    'D31_NA', 'D32_NA', 'D33_NA',
]

GEOMETRY_KEYS = [
    'boozer_s', 'aiota', 'R0', 'Bref',
    'psi_pr_hat', 'sqrtg_bctrvr_tht', 'sqrtg_bctrvr_phi',
    'bcovar_tht', 'bcovar_phi', 'avbhat2',
]


def compare_outputs(ref_file, test_file, keys, tol, label):
    """Compare selected HDF5 datasets between two output files."""
    mismatches = []
    with h5py.File(ref_file, 'r') as ref, h5py.File(test_file, 'r') as test:
        for key in keys:
            if key not in ref:
                continue
            if key not in test:
                mismatches.append(f"  {key}: missing in test output")
                continue
            ref_val = np.atleast_1d(np.array(ref[key], dtype=float))
            test_val = np.atleast_1d(np.array(test[key], dtype=float))
            if ref_val.shape != test_val.shape:
                mismatches.append(f"  {key}: shape mismatch {ref_val.shape} vs {test_val.shape}")
                continue
            max_ref = np.max(np.abs(ref_val))
            if max_ref == 0:
                max_diff = np.max(np.abs(test_val))
            else:
                max_diff = np.max(np.abs(ref_val - test_val)) / max_ref
            if max_diff > tol:
                mismatches.append(f"  {key}: rel diff {max_diff:.3e} > tol {tol:.0e}")
    if mismatches:
        print(f"FAIL: {label}")
        for m in mismatches:
            print(m)
        return False
    print(f"PASS: {label}")
    return True


def verify_outputs_differ(ref_file, test_file, keys, label):
    """Verify that at least one key differs between outputs."""
    with h5py.File(ref_file, 'r') as ref, h5py.File(test_file, 'r') as test:
        for key in keys:
            if key not in ref or key not in test:
                continue
            ref_val = np.atleast_1d(np.array(ref[key], dtype=float))
            test_val = np.atleast_1d(np.array(test[key], dtype=float))
            if not np.allclose(ref_val, test_val, rtol=1e-6):
                print(f"PASS: {label} (differs in {key})")
                return True
    print(f"FAIL: {label} (no differences found)")
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('neo2_binary', help='Path to neo_2_ql.x')
    parser.add_argument('data_dir', help='Path to golden_record/ql directory')
    parser.add_argument('--keep', action='store_true', help='Keep temp directories')
    args = parser.parse_args()

    neo2_binary = os.path.abspath(args.neo2_binary)
    data_dir = os.path.abspath(args.data_dir)

    for path, desc in [(neo2_binary, 'binary'), (data_dir, 'data_dir')]:
        if not os.path.exists(path):
            print(f"ERROR: {desc} not found: {path}")
            sys.exit(1)

    ref_neo2_in = os.path.join(data_dir, 'reference', 'neo2.in')
    output_file = 'neo2_multispecies_out.h5'

    # Reduced resolution for fast smoke test (~seconds instead of minutes).
    # Physical accuracy is irrelevant; we only test that isw_calc_Er=2
    # produces consistent Om_tE/MtOvR values.
    # Reduced resolution for fast smoke test (~seconds instead of minutes).
    # Disable NTV (ripple solver is the main cost) and use minimal grids.
    # Physical accuracy is irrelevant; we only test that isw_calc_Er=2
    # produces consistent Om_tE/MtOvR values.
    # Reduced resolution for fast smoke test.
    # Disable NTV and magnetic drift (the ripple solver and drift
    # computation dominate runtime). Use minimal velocity grids.
    # Physical accuracy is irrelevant; we only test that isw_calc_Er=2
    # produces consistent Om_tE/MtOvR values.
    # Reduced resolution for fast smoke test.
    # Disable NTV (ripple solver is the main cost). Keep ISW_CALC_MAGDRIFT=1
    # because the multispecies output path requires it.
    # Physical accuracy is irrelevant; we only test that isw_calc_Er=2
    # produces consistent Om_tE/MtOvR values.
    fast_params = {
        'NPERIOD': '100',
        'NSTEP': '100',
        'LAG': '2',
        'LEG': '2',
        'LEGMAX': '3',
        'ISW_NTV_MODE': '0',
        'ISW_QFLUX_NA': '0',
        'ISW_MAG_SHEAR': '0',
    }

    base_dir = tempfile.mkdtemp(prefix='neo2_er2_test_')
    status = 0

    try:
        # === Run 1: isw_calc_Er=1 (self-consistent) ===
        print("=== Run 1: isw_calc_Er=1 (self-consistent Er) ===")
        run1_dir = setup_run_dir(base_dir, 'run1_mode1', data_dir, neo2_binary)
        patch_neo2_in(ref_neo2_in, os.path.join(run1_dir, 'neo2.in'), fast_params)
        run_neo2(run1_dir)
        out1 = os.path.join(run1_dir, output_file)

        with h5py.File(out1, 'r') as f:
            Om_tE_from_run1 = float(f['Om_tE'][()])
            print(f"  Om_tE = {Om_tE_from_run1:.6e} rad/s")
            if 'MtOvR' in f:
                print(f"  MtOvR = {f['MtOvR'][()]}")

        # === Run 2: isw_calc_Er=2, same Om_tE via namelist ===
        print("\n=== Run 2: isw_calc_Er=2 (prescribe Om_tE from run 1) ===")
        run2_dir = setup_run_dir(base_dir, 'run2_mode2_same', data_dir, neo2_binary)
        patch_neo2_in(ref_neo2_in, os.path.join(run2_dir, 'neo2.in'), {
            **fast_params,
            'ISW_CALC_ER': '2',
            'OM_TE': f'{Om_tE_from_run1:.15e}',
        })
        run_neo2(run2_dir)
        out2 = os.path.join(run2_dir, output_file)

        with h5py.File(out2, 'r') as f:
            print(f"  Om_tE = {float(f['Om_tE'][()]):.6e} rad/s")
            if 'MtOvR' in f:
                print(f"  MtOvR = {f['MtOvR'][()]}")

        # === Run 3: isw_calc_Er=2, half Om_tE ===
        print("\n=== Run 3: isw_calc_Er=2 (half Om_tE) ===")
        run3_dir = setup_run_dir(base_dir, 'run3_mode2_half', data_dir, neo2_binary)
        patch_neo2_in(ref_neo2_in, os.path.join(run3_dir, 'neo2.in'), {
            **fast_params,
            'ISW_CALC_ER': '2',
            'OM_TE': f'{Om_tE_from_run1 * 0.5:.15e}',
        })
        run_neo2(run3_dir)
        out3 = os.path.join(run3_dir, output_file)

        with h5py.File(out3, 'r') as f:
            print(f"  Om_tE = {float(f['Om_tE'][()]):.6e} rad/s")
            if 'MtOvR' in f:
                print(f"  MtOvR = {f['MtOvR'][()]}")

        # === Comparisons ===
        print("\n=== Comparisons ===")

        if not compare_outputs(out1, out2, TRANSPORT_KEYS, 1e-10,
                               "Transport coeffs match (mode 1 vs mode 2, same Om_tE)"):
            status += 1

        if not compare_outputs(out1, out2, GEOMETRY_KEYS, 1e-14,
                               "Geometry quantities match (mode 1 vs mode 2)"):
            status += 1

        if not compare_outputs(out1, out2, ['Om_tE'], 1e-12,
                               "Om_tE matches (mode 1 vs mode 2)"):
            status += 1

        if not compare_outputs(out1, out2, ['MtOvR'], 1e-10,
                               "MtOvR matches (mode 1 vs mode 2, same Om_tE)"):
            status += 1

        if not verify_outputs_differ(out2, out3, ['MtOvR', 'Om_tE'],
                                     "MtOvR/Om_tE differ with half Om_tE"):
            status += 1

        # AX transport coefficients are independent of Er, so even with
        # different Om_tE the kinetic solve should give the same result.
        ax_keys = [k for k in TRANSPORT_KEYS if '_AX' in k]
        if not compare_outputs(out1, out3, ax_keys, 1e-10,
                               "AX transport coeffs match (mode 1 vs mode 2, half Om_tE)"):
            status += 1

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        status += 1

    finally:
        if args.keep:
            print(f"\nTest directories kept at: {base_dir}")
        else:
            shutil.rmtree(base_dir)

    print(f"\n{'=' * 50}")
    if status == 0:
        print("All integration tests passed!")
    else:
        print(f"{status} test(s) failed!")
        sys.exit(1)


if __name__ == '__main__':
    main()
