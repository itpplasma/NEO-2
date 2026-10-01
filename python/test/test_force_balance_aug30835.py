"""Force-balance ramp on two AUG 30835 surfaces from an independent NEO-2 run.

The reference data ``data/omte_reference_aug30835.npz`` was extracted from two
NEO-2-QL runs (e + D, axisymmetric AUG 30835 geometry at s = 0.253 and 0.498)
made with an earlier NEO-2 revision than the one used for
``neo2_ql_axisymmetric_multispecies_out.h5``; see
``data/regenerate_omte_reference_aug30835_fixture.py``. Both runs use the same
local plasma parameters (n, T, their gradients and Vphi); only the geometry
differs. Unlike the axisymmetric fixture used by ``test_force_balance.py``,
the ion row has sizeable cross-species coefficients and <E_par B> is nonzero.

Oracles are quantities written by the Fortran code: the stored ``Er`` and the
toroidal Mach number ``MtOvR``. The physics checks restate the AUG 30835
numbers quoted in ``DOC/ExtraDocuments/omte-force-balance-levels.tex``; run
this file with ``--table`` to print them.

Runs under pytest or directly as a script (used by ctest).
"""
import sys
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

from neo2_ql.force_balance import (
    POLOIDAL_ROTATION_K_LIMITS, er_level0_diamagnetic,
    er_level1_toroidal_rotation, er_level2_poloidal_rotation,
    er_level3_neo2_multispecies, omte_from_er,
    poloidal_rotation_coefficient_from_neo2, rigid_rotation_defect)

REFERENCE = (Path(__file__).resolve().parent / 'data'
             / 'omte_reference_aug30835.npz')


def _surfaces():
    """Return one keyword dict for ``er_level3_neo2_multispecies`` per surface.

    The npz stores per-surface arrays along axis 0 and species tags in
    ``row_ind_spec``/``col_ind_spec``; tags are mapped to zero-based indices.
    """
    d = np.load(REFERENCE)
    tags = [int(t) for t in d['species_tag']]
    index = {t: i for i, t in enumerate(tags)}
    spec_i = index[int(d['species_tag_Vphi'])]
    surfaces = []
    for s in range(d['boozer_s'].size):
        surfaces.append(dict(
            spec_i=spec_i,
            n=d['n_prof'][s], T=d['T_prof'][s],
            dn_ds=d['dn_ov_ds_prof'][s], dT_ds=d['dT_ov_ds_prof'][s],
            z=d['z_spec'],
            row=np.array([index[int(t)] for t in d['row_ind_spec'][s]]),
            col=np.array([index[int(t)] for t in d['col_ind_spec'][s]]),
            D31=d['D31_AX'][s], D32=d['D32_AX'][s], D33=d['D33_AX'][s],
            avEparB_ov_avb2=float(d['avEparB_ov_avb2'][s]),
            vphi=float(d['Vphi'][s]), aiota=float(d['aiota'][s]),
            av_nabla_stor=float(d['av_nabla_stor'][s]),
            sqrtg_bctrvr_tht=float(d['sqrtg_bctrvr_tht'][s]),
            bcovar_tht=float(d['bcovar_tht'][s]),
            bcovar_phi=float(d['bcovar_phi'][s]),
        ))
    extra = dict(boozer_s=d['boozer_s'], Er=d['Er'], MtOvR=d['MtOvR'],
                 m_spec=d['m_spec'])
    return surfaces, extra


def _ion(d):
    i = d['spec_i']
    return dict(n=d['n'][i], T=d['T'][i], dn_ds=d['dn_ds'][i],
                dT_ds=d['dT_ds'][i], z=d['z'][i],
                av_nabla_stor=d['av_nabla_stor'])


def _levels(d):
    """E_r [statV/cm] of every level on one surface, keyed by name."""
    ion = _ion(d)
    geom = dict(aiota=d['aiota'], bcovar_tht=d['bcovar_tht'],
                bcovar_phi=d['bcovar_phi'])
    k_neo2 = poloidal_rotation_coefficient_from_neo2(
        d['spec_i'], d['row'], d['col'], d['D31'], d['D32'])
    level1 = dict(ion, vphi=d['vphi'], sqrtg_bctrvr_tht=d['sqrtg_bctrvr_tht'])
    er3, _ = er_level3_neo2_multispecies(**d)
    return {
        'k_neo2': k_neo2,
        'level0': er_level0_diamagnetic(**ion),
        'level1': er_level1_toroidal_rotation(**level1),
        'level2_banana': er_level2_poloidal_rotation(
            **level1, **geom, k=POLOIDAL_ROTATION_K_LIMITS['banana']),
        'level2_neo2': er_level2_poloidal_rotation(**level1, **geom,
                                                   k=k_neo2),
        'level3': er3,
    }


def test_level3_replays_fortran_er_on_aug_surfaces():
    # Independent NEO-2 revision, two species with cross-species D31/D32 and
    # nonzero <E_par B>: the replay must still give the Fortran Er.
    surfaces, extra = _surfaces()
    for d, er_fortran in zip(surfaces, extra['Er']):
        assert_allclose(_levels(d)['level3'], er_fortran, rtol=1e-8)


def test_omte_matches_fortran_mach_number():
    # NEO-2 writes MtOvR = Om_tE / v_th,a with v_th,a = sqrt(2 T_a / m_a) for
    # every species a; Om_tE = c E_r / psi_pr must reproduce it.
    surfaces, extra = _surfaces()
    for s, d in enumerate(surfaces):
        omte = omte_from_er(extra['Er'][s], d['sqrtg_bctrvr_tht'])
        vth = np.sqrt(2.0 * d['T'] / extra['m_spec'])
        assert_allclose(omte / vth, extra['MtOvR'][s], rtol=1e-4)


def test_level2_with_neo2_k_is_close_and_banana_limit_fails():
    # With k = 5/2 - D32_ii/D31_ii from NEO-2, Level 2 overestimates |E_r| of
    # the full closure by 9-12 % (the rest is cross-species coupling and the
    # inductive drive). The large-aspect-ratio banana value k = 1.17 is
    # about twice NEO-2's k here (finite trapped fraction) and gives E_r of
    # the wrong sign; Levels 0 and 1 are off by factors 2-5.
    surfaces, _ = _surfaces()
    for d in surfaces:
        e = _levels(d)
        assert 0.5 < e['k_neo2'] < 0.65
        assert 0.0 < e['level2_neo2'] / e['level3'] - 1.0 < 0.15
        assert np.sign(e['level2_banana']) != np.sign(e['level3'])
        assert e['level1'] / e['level3'] > 2.0
        assert e['level0'] / e['level3'] > e['level1'] / e['level3']


def test_level2_is_level3_restricted_to_ion_ion_entry():
    # Without cross-species entries and D33, Level 3 differs from Level 2
    # with NEO-2's k only through the momentum defect of the D31_ii entry
    # (1.3 % and 0.4 % here).
    surfaces, _ = _surfaces()
    for d in surfaces:
        ii = (d['row'] == d['spec_i']) & (d['col'] == d['spec_i'])
        er3_ii, _ = er_level3_neo2_multispecies(**dict(
            d, row=d['row'][ii], col=d['col'][ii], D31=d['D31'][ii],
            D32=d['D32'][ii], D33=None))
        assert_allclose(_levels(d)['level2_neo2'], er3_ii, rtol=0.02)


def test_ion_row_nearly_conserves_momentum():
    # The momentum identity sum_b D31_ib Z_b e psi_pr/(c T_b) = B_phi holds
    # to better than 1 % on these runs, so the D31 term in the denominator
    # is the geometric factor discussed in the derivation.
    surfaces, _ = _surfaces()
    for d in surfaces:
        defect = rigid_rotation_defect(
            d['spec_i'], d['T'], d['z'], d['row'], d['col'], d['D31'],
            d['sqrtg_bctrvr_tht'], d['bcovar_phi'])
        assert abs(defect) < 1e-2


def print_table():
    surfaces, extra = _surfaces()
    names = ['level0', 'level1', 'level2_banana', 'level2_neo2', 'level3']
    rows = [_levels(d) for d in surfaces]
    print('boozer_s   ' + '  '.join(f'{s:10.4f}' for s in extra['boozer_s']))
    print('Er NEO-2   ' + '  '.join(f'{er:10.4f}' for er in extra['Er']))
    print('k NEO-2    ' + '  '.join(f"{r['k_neo2']:10.4f}" for r in rows))
    for name in names:
        print(f'{name:<11}' + '  '.join(
            f'{r[name]:10.4f} ({r[name] / er - 1.0:+7.1%})'
            for r, er in zip(rows, extra['Er'])))


if __name__ == '__main__':
    if '--table' in sys.argv:
        print_table()
        sys.exit(0)
    tests = [obj for name, obj in sorted(globals().items())
             if name.startswith('test_') and callable(obj)]
    for test in tests:
        test()
    print(f'All tests passed! ({len(tests)} AUG 30835 force-balance checks)')
