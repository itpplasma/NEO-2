"""Behavioral tests for the E_r / Om_tE force-balance ramp (issue #75).

Oracles are independent of the implementation: SI hand calculations,
laboratory-frame textbook formulas, asymptotic limits of published fits,
analytic limits in which two levels must coincide, and quantities written by
the Fortran code into a committed NEO-2-QL output file.

Runs under pytest or directly as a script (used by ctest).
"""
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

from neo2_ql.force_balance import (
    POLOIDAL_ROTATION_K_LIMITS, STATV_PER_CM_TO_V_PER_M,
    er_level0_diamagnetic, er_level1_toroidal_rotation,
    er_level2_poloidal_rotation, er_level3_neo2_multispecies,
    load_neo2_force_balance_inputs, omte_from_er,
    poloidal_rotation_coefficient_from_neo2,
    poloidal_rotation_coefficient_sauter, rigid_rotation_defect)

FIXTURE = (Path(__file__).resolve().parent / 'data'
           / 'neo2_ql_axisymmetric_multispecies_out.h5')
EV_TO_ERG = 1.602176634e-12
AV_NABLA_STOR = 0.02  # 1/cm; d/ds = d/dr / AV_NABLA_STOR

# Deuterium case in SI: n = 5e19 m^-3, T = 2 keV, dn/dr = -1e20 m^-4,
# dT/dr = -10 keV/m. By hand: E_r = (T dn/dr + n dT/dr) / (Z n) in V when T is
# in eV, i.e. (2e3 * -1e20 + 5e19 * -1e4) / 5e19 = -14000 V/m.
ION_SI = dict(
    n=5e13, T=2e3 * EV_TO_ERG, z=1.0, av_nabla_stor=AV_NABLA_STOR,
    dn_ds=-1e12 / AV_NABLA_STOR, dT_ds=-100.0 * EV_TO_ERG / AV_NABLA_STOR)
ER_DIA_SI = -14000.0
# The CGS charge 4.8032e-10 statC used by NEO-2 differs from the exact value
# by 1e-6 relative, c by 9e-6.
RTOL_CONST = 2e-5


def to_si(er_cgs):
    return er_cgs * STATV_PER_CM_TO_V_PER_M


def test_level0_matches_si_hand_calculation():
    assert_allclose(to_si(er_level0_diamagnetic(**ION_SI)), ER_DIA_SI,
                    rtol=RTOL_CONST)


def test_level0_sign_follows_charge():
    # Peaked ion pressure gives an inward (negative) E_r; for electrons with
    # the same profiles the diamagnetic E_r reverses sign.
    electron = dict(ION_SI, z=-1.0)
    assert er_level0_diamagnetic(**ION_SI) < 0.0
    assert_allclose(er_level0_diamagnetic(**electron),
                    -er_level0_diamagnetic(**ION_SI), rtol=1e-14)


def test_level1_matches_lab_frame_vphi_btheta():
    # Lab frame: E_r - p'/(Zen) = v_phi B_theta = Omega R B_theta.
    # Omega = 2e4 rad/s, R = 1.65 m, B_theta = 0.3 T -> 9900 V/m.
    omega, R_cm, btheta_gauss = 2e4, 165.0, 3000.0
    er = er_level1_toroidal_rotation(**ION_SI, vphi=omega,
                                     sqrtg_bctrvr_tht=R_cm * btheta_gauss)
    assert_allclose(to_si(er), ER_DIA_SI + 9900.0, rtol=RTOL_CONST)


def test_level2_large_aspect_ratio_limit():
    # For B_phi >> iota B_tht the poloidal term is -k dT/dr / (Z e):
    # -1.17 * (-1e4 V/m) = +11700 V/m in the banana regime.
    k = POLOIDAL_ROTATION_K_LIMITS['banana']
    er = er_level2_poloidal_rotation(**ION_SI, vphi=0.0,
                                     sqrtg_bctrvr_tht=4.95e5, aiota=0.5,
                                     bcovar_tht=1.0, bcovar_phi=3.3e7, k=k)
    assert_allclose(to_si(er), ER_DIA_SI + 11700.0, rtol=RTOL_CONST)


def test_levels_nest():
    geo = dict(sqrtg_bctrvr_tht=4.95e5, aiota=0.4, bcovar_tht=-1.2e5,
               bcovar_phi=-2.9e6)
    l0 = er_level0_diamagnetic(**ION_SI)
    l1 = er_level1_toroidal_rotation(**ION_SI, vphi=0.0,
                                     sqrtg_bctrvr_tht=geo['sqrtg_bctrvr_tht'])
    l2 = er_level2_poloidal_rotation(**ION_SI, vphi=0.0, k=0.0, **geo)
    assert_allclose([l1, l2], [l0, l0], rtol=1e-14)


def test_sauter_coefficient_reaches_regime_limits():
    # Banana limit (f_t -> 0, nu* -> 0) and Pfirsch-Schlueter limit
    # (nu* -> inf) of the Sauter fit must reproduce the asymptotic values.
    assert_allclose(poloidal_rotation_coefficient_sauter(1e-8, 0.0),
                    POLOIDAL_ROTATION_K_LIMITS['banana'], rtol=1e-6)
    assert_allclose(poloidal_rotation_coefficient_sauter(0.5, 1e8),
                    POLOIDAL_ROTATION_K_LIMITS['pfirsch-schlueter'],
                    rtol=1e-3)


def test_sauter_coefficient_at_finite_collisionality():
    # Hand evaluation of Sauter et al. (1999, erratum 2002), f_t = 0.5,
    # nu*_i = 1: alpha_0 = -0.585 / 0.8425 = -0.694362;
    # (alpha_0 + 0.1875) / 1.5 = -0.337908; + 0.315/64 = -0.332986;
    # / (1 + 0.15/64) = -0.332208, so k = 0.332208.
    assert_allclose(poloidal_rotation_coefficient_sauter(0.5, 1.0), 0.332208,
                    rtol=2e-6)


def _single_ion_coefficients(k, T, z, psi_pr, bcovar_phi):
    # Momentum-conserving single ion species: a rigid rotation carries
    # <V_par B> = omega B_phi, i.e. D31 Z e psi_pr / (c T) = B_phi.
    from neo2_ql.force_balance import C_CGS, E_CGS
    d31 = bcovar_phi * C_CGS * T / (z * E_CGS * psi_pr)
    return np.array([d31]), np.array([(2.5 - k) * d31])


def test_level3_single_ion_reduces_to_level2():
    geo = dict(aiota=0.46, bcovar_tht=-1.23e5, bcovar_phi=-2.93e6,
               sqrtg_bctrvr_tht=4.17e5)
    k, vphi = 0.57, 2.24e4
    d31, d32 = _single_ion_coefficients(k, ION_SI['T'], 1.0,
                                        geo['sqrtg_bctrvr_tht'],
                                        geo['bcovar_phi'])
    er3, terms = er_level3_neo2_multispecies(
        0, [ION_SI['n']], [ION_SI['T']], [ION_SI['dn_ds']],
        [ION_SI['dT_ds']], [1.0], [0], [0], d31, d32, vphi,
        av_nabla_stor=AV_NABLA_STOR, **geo)
    er2 = er_level2_poloidal_rotation(**ION_SI, vphi=vphi, k=k, **geo)
    assert_allclose(er3, er2, rtol=1e-12)
    assert_allclose(poloidal_rotation_coefficient_from_neo2(
        0, [0], [0], d31, d32), k, rtol=1e-14)
    # Denominator amplification is then purely geometric.
    assert_allclose((terms['den_base'] + terms['den_d31']) / terms['den_base'],
                    1.0 + geo['bcovar_phi']
                    / (geo['aiota'] * geo['bcovar_tht']), rtol=1e-12)


# --- Benchmarks against a committed NEO-2-QL run ----------------------------
# Fixture: two species (e, D), one surface, isw_calc_Er = 1, isw_Vphi_loc = 0,
# written by a NEO-2-QL build that also stores dn/dT_spec_ov_ds and Vphi.


def _fixture():
    data = load_neo2_force_balance_inputs(FIXTURE)
    return data, data.pop('Er_stored')


def _ion_kwargs(d):
    i = d['spec_i']
    return dict(n=d['n'][i], T=d['T'][i], dn_ds=d['dn_ds'][i],
                dT_ds=d['dT_ds'][i], z=d['z'][i],
                av_nabla_stor=d['av_nabla_stor'])


def test_level3_replays_fortran_er():
    # Fortran stores Er from the last fixed-point iteration of
    # (Er, <E_par B>/<B^2>) and the converged <E_par B>/<B^2>; observed
    # replay error is 2e-11.
    d, er_stored = _fixture()
    er, _ = er_level3_neo2_multispecies(**d)
    assert_allclose(er, er_stored, rtol=1e-9)


def test_loader_maps_species_tags_not_positions():
    # Relabel the fixture species with non-contiguous tags (1 -> 7, 2 -> 3);
    # the physics, and hence E_r, must not change.
    import shutil
    import tempfile
    import h5py
    d, er_stored = _fixture()
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / 'relabelled.h5'
        shutil.copy(FIXTURE, path)
        relabel = {1: 7, 2: 3}
        with h5py.File(path, 'r+') as f:
            for key in ('species_tag', 'row_ind_spec', 'col_ind_spec'):
                f[key][...] = [relabel[int(t)] for t in f[key][()]]
            f['species_tag_Vphi'][...] = relabel[int(f['species_tag_Vphi'][()])]
        relabelled = load_neo2_force_balance_inputs(path)
    relabelled.pop('Er_stored')
    assert relabelled['spec_i'] == d['spec_i']
    er, _ = er_level3_neo2_multispecies(**relabelled)
    assert_allclose(er, er_stored, rtol=1e-9)


def test_omte_matches_fortran_mach_number():
    # NEO-2 writes MtOvR_spec = Om_tE / sqrt(2 T / m) (er_rotation_mod).
    import h5py
    d, er_stored = _fixture()
    with h5py.File(FIXTURE, 'r') as f:
        mtovr, m = f['MtOvR'][()], f['m_spec'][()]
    omte = omte_from_er(er_stored, d['sqrtg_bctrvr_tht'])
    assert_allclose(omte, mtovr * np.sqrt(2.0 * d['T'] / m), rtol=1e-10)


def test_fixture_coefficients_conserve_momentum():
    # Physics check on NEO-2's D31 row of the rotation species, and with it
    # the size of the "D31 denominator correction" (52x on this surface).
    d, _ = _fixture()
    defect = rigid_rotation_defect(d['spec_i'], d['T'], d['z'], d['row'],
                                   d['col'], d['D31'], d['sqrtg_bctrvr_tht'],
                                   d['bcovar_phi'])
    assert abs(defect) < 1e-3
    _, terms = er_level3_neo2_multispecies(**d)
    ratio = (terms['den_base'] + terms['den_d31']) / terms['den_base']
    assert_allclose(ratio, 1.0 + d['bcovar_phi']
                    / (d['aiota'] * d['bcovar_tht']), rtol=1e-3)


def test_level2_with_neo2_k_tracks_level3_without_inductive_drive():
    # Level 2 with k = 5/2 - D32_ii/D31_ii neglects the electron cross
    # coefficients D31_ie, D32_ie and the momentum defect; on this surface
    # they change E_r by 3.3 %. The inductive (D33) term is excluded from the comparison
    # because Level 2 has no parallel electric field.
    d, _ = _fixture()
    k = poloidal_rotation_coefficient_from_neo2(d['spec_i'], d['row'],
                                                d['col'], d['D31'], d['D32'])
    er2 = er_level2_poloidal_rotation(
        **_ion_kwargs(d), vphi=d['vphi'], k=k,
        sqrtg_bctrvr_tht=d['sqrtg_bctrvr_tht'], aiota=d['aiota'],
        bcovar_tht=d['bcovar_tht'], bcovar_phi=d['bcovar_phi'])
    er3, _ = er_level3_neo2_multispecies(**dict(d, D33=None))
    assert_allclose(er2, er3, rtol=0.05)
    # Without poloidal rotation the reduced levels are far off: Level 0 has
    # the wrong sign, Level 1 misses 40 %.
    er0 = er_level0_diamagnetic(**_ion_kwargs(d))
    er1 = er_level1_toroidal_rotation(**_ion_kwargs(d), vphi=d['vphi'],
                                      sqrtg_bctrvr_tht=d['sqrtg_bctrvr_tht'])
    assert np.sign(er0) != np.sign(er3)
    assert abs(er1 / er3 - 1.0) > 0.3


if __name__ == '__main__':
    tests = [obj for name, obj in sorted(globals().items())
             if name.startswith('test_') and callable(obj)]
    for test in tests:
        test()
    print(f'All tests passed! ({len(tests)} force-balance checks)')
