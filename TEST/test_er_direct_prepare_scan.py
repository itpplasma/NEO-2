#!/usr/bin/env python3
"""Prepare-scan test: isw_calc_Er=2 must not require or read V_phi.

Runs the NEO-2-QL multi-species scan preparation (isw_multispecies_init=1)
on a synthetic multi-species HDF5 input and checks the per-surface neo2.in
files it writes. The oracle is the synthetic input itself: every surface
directory must carry the Om_tE, T, n and boozer_s values of its own radial
point.

Cases:
  1. isw_calc_Er=2, input without any V_phi dataset: preparation succeeds and
     each surface gets the prescribed Om_tE.
  2. isw_calc_Er=2, input with a V_phi profile and invalid V_phi settings in
     the top-level neo2.in: identical per-surface Om_tE and species data as
     case 1, and neutral V_phi settings in every surface (nothing of the
     V_phi input is propagated).
  3. isw_calc_Er=1, input with V_phi: V_phi of each surface is propagated
     (the V_phi-driven mode keeps its behaviour).
  4. isw_calc_Er=1, input without V_phi: preparation must fail, since that
     mode needs V_phi.

Usage: test_er_direct_prepare_scan.py <neo_2_ql.x> <neo.in fixture>
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

import h5py
import numpy as np

NRAD = 3
BOOZER_S = np.array([0.25, 0.5, 0.75])
OM_TE = np.array([-1.5e4, 2.0e3, 3.25e4])
VPHI = np.array([7.0e3, -9.0e3, 1.1e4])
Z = np.array([-1.0, 1.0])
M = np.array([9.1093836e-28, 3.3436e-24])
T = np.array([[2.0e-9, 1.5e-9, 1.0e-9], [2.5e-9, 1.75e-9, 0.75e-9]])
N = np.array([[3.0e13, 2.0e13, 1.0e13], [2.9e13, 1.8e13, 0.7e13]])

NEO2_IN = """&multi_spec
 lsw_multispecies = .true.
 isw_multispecies_init = 1
 fname_multispec_in = 'multispec.h5'
 isw_calc_er = {isw_calc_er}
{extra}/
&settings
/
&collision
/
&binsplit
/
&propagator
/
&plotting
/
&ntv_input
 in_file_pert = 'test_pert.bc'
/
"""


def write_multispec(path, with_vphi):
    with h5py.File(path, 'w') as f:
        f['num_radial_pts'] = np.array([NRAD], dtype=np.int32)
        f['num_species'] = np.array([2], dtype=np.int32)
        f['species_tag'] = np.array([1, 2], dtype=np.int32)
        species_def = np.zeros((2, 2, NRAD))
        species_def[0] = Z[:, None]
        species_def[1] = M[:, None]
        f['species_def'] = species_def
        f['boozer_s'] = BOOZER_S
        f['rho_pol'] = np.sqrt(BOOZER_S)
        f['rel_stages'] = np.full(NRAD, 2, dtype=np.int32)
        f['T_prof'] = T
        f['dT_ov_ds_prof'] = -T
        f['n_prof'] = N
        f['dn_ov_ds_prof'] = -N
        f['kappa_prof'] = np.full((2, NRAD), 1.0e-5)
        f['Om_tE'] = OM_TE
        if with_vphi:
            f['Vphi'] = VPHI
            f['species_tag_Vphi'] = np.array([2], dtype=np.int32)
            f['isw_Vphi_loc'] = np.array([0], dtype=np.int32)


def run_prepare(binary, neo_in, workdir, isw_calc_er, with_vphi, extra=''):
    os.makedirs(workdir)
    shutil.copy(neo_in, os.path.join(workdir, 'neo.in'))
    with open(os.path.join(workdir, 'neo2.in'), 'w') as f:
        f.write(NEO2_IN.format(isw_calc_er=isw_calc_er, extra=extra))
    write_multispec(os.path.join(workdir, 'multispec.h5'), with_vphi)
    return subprocess.run([binary], cwd=workdir, capture_output=True,
                          text=True, timeout=120)


def surface_dir(workdir, s):
    whole, frac = f'{s:10.5f}'.strip().split('.')
    return os.path.join(workdir, f'es_{whole}p{frac}')


def namelist_values(text, key):
    """Values of KEY in gfortran namelist output, with r*x repeats expanded."""
    match = re.search(rf'^\s*{key}\s*=(.*?)(?=^\s*[A-Z_0-9]+\s*=|^\s*/)',
                      text, flags=re.IGNORECASE | re.MULTILINE | re.DOTALL)
    if match is None:
        raise AssertionError(f'{key} not found in generated namelist')
    values = []
    for item in match.group(1).replace('\n', ' ').split(','):
        item = item.strip()
        if not item:
            continue
        count, value = item.split('*') if '*' in item else ('1', item)
        if value.upper() in ('T', '.TRUE.', 'F', '.FALSE.'):
            value = '1' if value.upper().strip('.') == 'T' else '0'
        values += [float(value.upper().replace('D', 'E'))] * int(count)
    return np.array(values)


def read_surfaces(workdir):
    out = []
    for s in BOOZER_S:
        with open(os.path.join(surface_dir(workdir, s), 'neo2.in')) as f:
            out.append(f.read())
    return out


def same(actual, expected):
    expected = np.atleast_1d(expected)
    return actual.shape == expected.shape and np.allclose(
        actual, expected, rtol=1e-14, atol=0.0)


def check_species_and_omte(texts, label):
    for k, text in enumerate(texts):
        assert same(namelist_values(text, 'BOOZER_S'), BOOZER_S[k]), \
            f'{label}: boozer_s at surface {k}'
        assert same(namelist_values(text, 'OM_TE'), OM_TE[k]), \
            f'{label}: Om_tE at surface {k}'
        assert same(namelist_values(text, 'T_VEC'), T[:, k]), \
            f'{label}: t_vec at surface {k}'
        assert same(namelist_values(text, 'N_VEC'), N[:, k]), \
            f'{label}: n_vec at surface {k}'
        assert np.all(namelist_values(text, 'ISW_CALC_ER') == 2), label


def main():
    binary = os.path.abspath(sys.argv[1])
    neo_in = os.path.abspath(sys.argv[2])
    base = tempfile.mkdtemp(prefix='er_direct_')
    failures = 0

    def case(name, fn):
        nonlocal failures
        try:
            fn()
            print(f'PASS {name}')
        except AssertionError as err:
            failures += 1
            print(f'FAIL {name}: {err}')

    def mode2(with_vphi, extra=''):
        workdir = os.path.join(base, f'mode2_vphi{int(with_vphi)}')
        res = run_prepare(binary, neo_in, workdir, 2, with_vphi, extra)
        assert res.returncode == 0, res.stdout[-2000:] + res.stderr[-2000:]
        return read_surfaces(workdir)

    def case1():
        check_species_and_omte(mode2(False), 'no Vphi')

    def case2():
        poison = (' isw_vphi_loc = 3\n species_tag_vphi = 1\n'
                  ' vphi = -1.0e12\n boozer_theta_vphi = 2.0\n')
        with_vphi = mode2(True, poison)
        check_species_and_omte(with_vphi, 'with Vphi')
        for k, text in enumerate(with_vphi):
            for key in ('VPHI', 'ISW_VPHI_LOC', 'SPECIES_TAG_VPHI',
                        'BOOZER_THETA_VPHI', 'R_VPHI', 'Z_VPHI'):
                assert same(namelist_values(text, key), 0.0), \
                    f'{key} of surface {k} not neutral in isw_calc_Er=2'

    def case3():
        workdir = os.path.join(base, 'mode1_vphi1')
        res = run_prepare(binary, neo_in, workdir, 1, True)
        assert res.returncode == 0, res.stdout[-2000:] + res.stderr[-2000:]
        for k, text in enumerate(read_surfaces(workdir)):
            assert same(namelist_values(text, 'VPHI'), VPHI[k]), \
                f'Vphi at surface {k}'
            assert np.all(namelist_values(text, 'SPECIES_TAG_VPHI') == 2)
            assert np.all(namelist_values(text, 'ISW_CALC_ER') == 1)

    def case4():
        workdir = os.path.join(base, 'mode1_vphi0')
        res = run_prepare(binary, neo_in, workdir, 1, False)
        created = os.path.isdir(surface_dir(workdir, BOOZER_S[0]))
        assert res.returncode != 0 or not created, \
            'isw_calc_Er=1 prepared a scan without V_phi'

    case('mode 2 without V_phi', case1)
    case('mode 2 ignores V_phi', case2)
    case('mode 1 propagates V_phi', case3)
    case('mode 1 requires V_phi', case4)

    if failures == 0:
        shutil.rmtree(base)
        print('All tests passed!')
    else:
        print(f'{failures} case(s) FAILED; work dirs kept in {base}')
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
