#!/usr/bin/env python3
"""Per-surface solve test: isw_calc_Er=2 runs without V_phi and ignores it.

Uses the golden-record QL case (private data, gitlab.tugraz.at/plasma/data,
TESTS/NEO-2/golden_record/ql) at reduced resolution. Runs neo_2_ql.x four
times on one flux surface with two MPI ranks:

  m1      isw_calc_Er=1, Er from force balance with the measured V_phi
  m2      isw_calc_Er=2, Om_tE of m1, every V_phi entry removed from neo2.in
  m2_bad  as m2, but with invalid/poisoned V_phi settings
  m2_neg  as m2 with -Om_tE

Checks (oracles independent of the code path under test):
  - m2 and m2_bad complete and all output datasets are bitwise identical,
    so no V_phi setting is read in the solve.
  - Er = Om_tE*aiota*sqrtg_bctrvr_phi/c and MtOvR_a = Om_tE/sqrt(2 T_a/m_a),
    evaluated here from the output geometry and the input T_a, m_a.
  - m2 reproduces m1's Er, MtOvR and transport coefficients.
  - The V_phi of the measured species in m2 plus m1's Ware-pinch part equals
    the V_phi that m1 was given: a prescribed Er consistent with V_phi
    reproduces V_phi.
  - Flipping the sign of Om_tE flips Er and MtOvR and leaves the AX
    coefficients unchanged.

Usage: test_er_direct_solve.py <neo_2_ql.x> <golden_record/ql directory>
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

import h5py
import numpy as np

C_CGS = 2.9979e10  # speed of light as in NEO-2 ntv_mod (c = 2.9979e10)
OUT = 'neo2_multispecies_out.h5'
FAST = {'NPERIOD': '100', 'NSTEP': '100', 'LAG': '2', 'LEG': '2',
        'LEGMAX': '3'}
VPHI_KEYS = ('VPHI', 'SPECIES_TAG_VPHI', 'ISW_VPHI_LOC', 'R_VPHI', 'Z_VPHI',
             'BOOZER_THETA_VPHI')
POISON = {'VPHI': '-1.0e12', 'SPECIES_TAG_VPHI': '1', 'ISW_VPHI_LOC': '3',
          'BOOZER_THETA_VPHI': '2.0'}
D_KEYS = [f'D{i}{j}_{k}' for k in ('AX', 'NA') for i in (1, 2, 3)
          for j in (1, 2, 3)]


def namelist_value(text, key):
    m = re.search(rf'^\s*{key}\s*=\s*([^,/!\n]+)', text, re.I | re.M)
    return m.group(1).strip()


def namelist_array(text, key):
    m = re.search(rf'^\s*{key}\s*=(.*?)(?=^\s*[A-Z_0-9]+\s*=|^\s*/)', text,
                  re.I | re.M | re.S)
    values = []
    for item in m.group(1).replace('\n', ' ').split(','):
        item = item.strip()
        if item:
            count, value = item.split('*') if '*' in item else ('1', item)
            values += [float(value.upper().replace('D', 'E'))] * int(count)
    return np.array(values)


def patch(text, replace, drop=()):
    for key in drop:
        text = re.sub(rf'^\s*{key}\s*=.*\n', '', text, flags=re.I | re.M)
    for key, value in replace.items():
        new = re.sub(rf'^(\s*{key}\s*=\s*)[^,/!\n]+', rf'\g<1>{value} ',
                     text, flags=re.I | re.M)
        if new == text:  # key absent: add to the namelist that owns it
            group = 'NTV_INPUT' if key == 'OM_TE' else 'MULTI_SPEC'
            new = re.sub(rf'(&{group}.*?)(^\s*/)', rf'\1 {key}={value},\n/',
                         text, count=1, flags=re.I | re.M | re.S)
        text = new
    return text


def run(base, name, data_dir, binary, text):
    d = os.path.join(base, name)
    os.makedirs(d)
    for f in ('test_axi.bc', 'test_pert.bc'):
        os.symlink(os.path.join(data_dir, f), os.path.join(d, f))
    os.symlink(os.path.join(data_dir, 'reference', 'neo.in'),
               os.path.join(d, 'neo.in'))
    with open(os.path.join(d, 'neo2.in'), 'w') as f:
        f.write(text)
    res = subprocess.run(['mpiexec', '-np', '2', binary], cwd=d,
                         capture_output=True, text=True, timeout=1800,
                         env={**os.environ, 'OMP_NUM_THREADS': '1'})
    if res.returncode != 0 or not os.path.exists(os.path.join(d, OUT)):
        raise AssertionError(f'{name} failed:\n{res.stdout[-1500:]}'
                             f'\n{res.stderr[-1500:]}')
    with h5py.File(os.path.join(d, OUT), 'r') as f:
        return {k: np.atleast_1d(f[k][()]) for k in f
                if isinstance(f[k], h5py.Dataset)}


def rel(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300)


def main():
    binary = os.path.abspath(sys.argv[1])
    data_dir = os.path.abspath(sys.argv[2])
    ref = open(os.path.join(data_dir, 'reference', 'neo2.in')).read()
    ref = patch(ref, {**FAST, 'ISW_CALC_ER': '1'})
    vphi_in = float(namelist_value(ref, 'VPHI'))
    tag_vphi = int(namelist_value(ref, 'SPECIES_TAG_VPHI'))
    t_in = namelist_array(ref, 'T_VEC')
    m_in = namelist_array(ref, 'M_VEC')

    base = tempfile.mkdtemp(prefix='er_direct_solve_')
    checks = []

    def check(name, ok):
        checks.append(ok)
        print(('PASS ' if ok else 'FAIL ') + name, flush=True)

    m1 = run(base, 'm1', data_dir, binary, ref)
    om_te = float(m1['Om_tE'][0])
    mode2 = {'ISW_CALC_ER': '2', 'OM_TE': f'{om_te:.17e}'}
    m2 = run(base, 'm2', data_dir, binary, patch(ref, mode2, VPHI_KEYS))
    m2_bad = run(base, 'm2_bad', data_dir, binary,
                 patch(ref, {**mode2, **POISON}))
    m2_neg = run(base, 'm2_neg', data_dir, binary,
                 patch(ref, {'ISW_CALC_ER': '2', 'OM_TE': f'{-om_te:.17e}'},
                       VPHI_KEYS))

    check('mode 2 output independent of V_phi settings',
          m2.keys() == m2_bad.keys()
          and all(np.array_equal(m2[k], m2_bad[k]) for k in m2))
    er_expected = om_te * m2['aiota'][0] * m2['sqrtg_bctrvr_phi'][0] / C_CGS
    check('Er = Om_tE*iota*sqrt(g)B^phi/c', rel(m2['Er'], er_expected) < 1e-12)
    check('MtOvR_a = Om_tE/v_Ta',
          rel(m2['MtOvR'], om_te / np.sqrt(2.0 * t_in / m_in)) < 1e-12)
    check('mode 2 reproduces mode-1 Er and MtOvR',
          rel(m2['Er'], m1['Er']) < 1e-12 and rel(m2['MtOvR'], m1['MtOvR'])
          < 1e-12)
    missing = [k for k in D_KEYS for o in (m1, m2, m2_neg) if k not in o]
    check(f'transport coefficients present {missing}', not missing)
    check('transport coefficients unchanged',
          not missing and all(rel(m2[k], m1[k]) < 1e-12 for k in D_KEYS))
    ispec = list(m1['species_tag']).index(tag_vphi)
    vphi_rebuilt = m2['VphiB_spec'][ispec] + m1['VphiB_Ware_spec'][ispec]
    check(f'V_phi rebuilt from prescribed Er ({vphi_rebuilt:.6e} vs '
          f'{vphi_in:.6e})', rel(vphi_rebuilt, vphi_in) < 1e-9)
    check('sign of Om_tE flips Er and MtOvR',
          rel(m2_neg['Er'], -m2['Er']) < 1e-14
          and rel(m2_neg['MtOvR'], -m2['MtOvR']) < 1e-14
          and not missing
          and all(rel(m2_neg[k], m2[k]) < 1e-12 for k in D_KEYS
                  if k.endswith('_AX')))

    if all(checks):
        shutil.rmtree(base)
        print('All tests passed!')
        return 0
    print(f'FAIL: {checks.count(False)} check(s); run dirs kept in {base}')
    return 1


if __name__ == '__main__':
    try:
        sys.exit(main())
    except AssertionError as err:
        print(f'FAIL: {err}')
        sys.exit(1)
