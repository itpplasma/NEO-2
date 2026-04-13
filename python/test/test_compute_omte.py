"""Tests for Om_tE force balance computation."""

import os
import numpy as np
from neo2_ql.compute_omte import (
    compute_omte_diamagnetic,
    compute_omte_force_balance,
    compute_omte_toroidal_rotation,
    compute_poloidal_rotation_neoclassical,
    compute_omte_neoclassical_poloidal,
    C_CGS,
    E_CGS,
)


# --- Unit tests with analytic profiles ---

def test_uniform_profiles_give_zero():
    """Uniform density and temperature means zero pressure gradient -> Om_tE=0."""
    n = 1e13
    T = 1e-9
    dn_ds = 0.0
    dT_ds = 0.0
    z = 1.0
    aiota = 0.5
    sqrtg_phi = 1e6
    av_nabla = 0.01

    Om_tE, Er = compute_omte_diamagnetic(
        n, T, dn_ds, dT_ds, z, aiota, sqrtg_phi, av_nabla
    )
    assert Om_tE == 0.0
    assert Er == 0.0


def test_negative_density_gradient_gives_negative_omte():
    """Negative dn/ds with positive geometry gives negative Om_tE (inward E_r)."""
    n = 1e13
    T = 1e-9
    dn_ds = -1e13   # falling density
    dT_ds = 0.0
    z = 1.0
    aiota = 0.5
    sqrtg_phi = 1e6  # positive
    av_nabla = 0.01

    Om_tE, Er = compute_omte_diamagnetic(
        n, T, dn_ds, dT_ds, z, aiota, sqrtg_phi, av_nabla
    )
    assert Om_tE < 0.0
    assert Er < 0.0


def test_analytic_linear_profiles():
    """Analytic test: linear n(s), constant T -> known dp/ds = T * dn/ds."""
    n = 2e13         # density at evaluation point
    T = 5e-9         # constant temperature
    dn_ds = -3e13    # linear density gradient
    dT_ds = 0.0
    z = 1.0
    aiota = 0.4
    sqrtg_phi = 5e5
    av_nabla = 0.02

    Om_tE, Er = compute_omte_diamagnetic(
        n, T, dn_ds, dT_ds, z, aiota, sqrtg_phi, av_nabla
    )

    # Manual calculation
    dp_ds = T * dn_ds  # = 5e-9 * (-3e13) = -1.5e5
    dp_dr = dp_ds * av_nabla  # = -1.5e5 * 0.02 = -3000
    Er_expected = dp_dr / (n * z * E_CGS)  # = -3000 / (2e13 * 4.803e-10)
    Om_tE_expected = C_CGS * Er_expected / (aiota * sqrtg_phi)

    assert np.isclose(Er, Er_expected, rtol=1e-12)
    assert np.isclose(Om_tE, Om_tE_expected, rtol=1e-12)


def test_temperature_gradient_contribution():
    """Both dn/ds and dT/ds contribute to dp/ds = T*dn/ds + n*dT/ds."""
    n = 1e13
    T = 2e-9
    dn_ds = -1e13
    dT_ds = -1e-9
    z = 1.0
    aiota = 0.5
    sqrtg_phi = 1e6
    av_nabla = 0.01

    Om_tE, Er = compute_omte_diamagnetic(
        n, T, dn_ds, dT_ds, z, aiota, sqrtg_phi, av_nabla
    )

    dp_ds = T * dn_ds + n * dT_ds  # = 2e-9*(-1e13) + 1e13*(-1e-9) = -3e4
    dp_dr = dp_ds * av_nabla
    Er_expected = dp_dr / (n * z * E_CGS)
    Om_tE_expected = C_CGS * Er_expected / (aiota * sqrtg_phi)

    assert np.isclose(Om_tE, Om_tE_expected, rtol=1e-12)


def test_force_balance_without_vphi_matches_diamagnetic():
    """Unified force-balance API must preserve Level 0 behavior."""
    inputs = {
        'n': 1e13,
        'T': 2e-9,
        'dn_ds': -1e13,
        'dT_ds': -1e-9,
        'z': 1.0,
        'aiota': 0.5,
        'sqrtg_bctrvr_phi': 1e6,
        'av_nabla_stor': 0.01,
    }

    om_dia, er_dia = compute_omte_diamagnetic(**inputs)
    om_full, er_full = compute_omte_force_balance(**inputs)

    assert np.isclose(om_full, om_dia, rtol=1e-12)
    assert np.isclose(er_full, er_dia, rtol=1e-12)


def test_charge_number_scaling():
    """Om_tE scales as 1/Z for the same pressure gradient."""
    n = 1e13
    T = 2e-9
    dn_ds = -1e13
    dT_ds = 0.0
    aiota = 0.5
    sqrtg_phi = 1e6
    av_nabla = 0.01

    Om_tE_z1, _ = compute_omte_diamagnetic(
        n, T, dn_ds, dT_ds, 1.0, aiota, sqrtg_phi, av_nabla
    )
    Om_tE_z2, _ = compute_omte_diamagnetic(
        n, T, dn_ds, dT_ds, 2.0, aiota, sqrtg_phi, av_nabla
    )

    assert np.isclose(Om_tE_z1 / Om_tE_z2, 2.0, rtol=1e-12)


def test_array_input():
    """Function works element-wise on arrays (multiple surfaces)."""
    n = np.array([1e13, 2e13, 3e13])
    T = np.array([3e-9, 2e-9, 1e-9])
    dn_ds = np.array([-1e13, -2e13, -3e13])
    dT_ds = np.array([-1e-9, -1e-9, -1e-9])
    z = 1.0
    aiota = np.array([0.4, 0.5, 0.6])
    sqrtg_phi = np.array([5e5, 7e5, 9e5])
    av_nabla = np.array([0.01, 0.015, 0.02])

    Om_tE, Er = compute_omte_diamagnetic(
        n, T, dn_ds, dT_ds, z, aiota, sqrtg_phi, av_nabla
    )

    assert Om_tE.shape == (3,)
    assert Er.shape == (3,)

    # Check each element individually
    for i in range(3):
        Om_i, Er_i = compute_omte_diamagnetic(
            n[i], T[i], dn_ds[i], dT_ds[i], z,
            aiota[i], sqrtg_phi[i], av_nabla[i]
        )
        assert np.isclose(Om_tE[i], Om_i, rtol=1e-12)
        assert np.isclose(Er[i], Er_i, rtol=1e-12)


def test_toroidal_rotation_contribution():
    """Level 1 adds the toroidal-rotation term v_phi * B_theta / c."""
    om_tE, er = compute_omte_toroidal_rotation(
        n=1e13,
        T=1e-9,
        dn_ds=0.0,
        dT_ds=0.0,
        z=1.0,
        aiota=0.5,
        sqrtg_bctrvr_phi=1e6,
        av_nabla_stor=0.01,
        v_phi=2e7,
        b_theta=-1e3,
    )

    er_expected = 2e7 * (-1e3) / C_CGS
    om_expected = C_CGS * er_expected / (0.5 * 1e6)

    assert np.isclose(er, er_expected, rtol=1e-12)
    assert np.isclose(om_tE, om_expected, rtol=1e-12)


def test_neoclassical_poloidal_rotation_formula():
    """The Level 2 poloidal estimate should match its analytic definition."""
    v_theta = compute_poloidal_rotation_neoclassical(
        dT_ds=-2e-9,
        z=1.0,
        b_phi=-2e4,
        av_nabla_stor=0.02,
        k_i=-1.17,
    )
    expected = -1.17 * (-2e-9 * 0.02) / (E_CGS * -2e4)
    assert np.isclose(v_theta, expected, rtol=1e-12)


def test_invalid_zero_density_raises():
    """Public API should fail explicitly on zero-density input."""
    try:
        compute_omte_diamagnetic(
            n=0.0,
            T=1e-9,
            dn_ds=-1e13,
            dT_ds=0.0,
            z=1.0,
            aiota=0.5,
            sqrtg_bctrvr_phi=1e6,
            av_nabla_stor=0.01,
        )
    except ValueError as exc:
        assert str(exc) == 'n must be nonzero to compute E_r and Om_tE'
    else:
        raise AssertionError('Expected ValueError for zero density')


def test_missing_toroidal_rotation_pair_raises():
    """v_phi and b_theta must be supplied together."""
    try:
        compute_omte_force_balance(
            n=1e13,
            T=1e-9,
            dn_ds=0.0,
            dT_ds=0.0,
            z=1.0,
            aiota=0.5,
            sqrtg_bctrvr_phi=1e6,
            av_nabla_stor=0.01,
            v_phi=1e7,
        )
    except ValueError as exc:
        assert str(exc) == 'v_phi and b_theta must be provided together'
    else:
        raise AssertionError('Expected ValueError for incomplete Level 1 inputs')


def test_missing_poloidal_rotation_pair_raises():
    """v_theta and b_phi must be supplied together."""
    try:
        compute_omte_force_balance(
            n=1e13,
            T=1e-9,
            dn_ds=0.0,
            dT_ds=0.0,
            z=1.0,
            aiota=0.5,
            sqrtg_bctrvr_phi=1e6,
            av_nabla_stor=0.01,
            v_theta=1e5,
        )
    except ValueError as exc:
        assert str(exc) == 'v_theta and b_phi must be provided together'
    else:
        raise AssertionError('Expected ValueError for incomplete Level 2 inputs')


# --- E2E test against NEO-2 reference ---

FIXTURE = os.path.join(os.path.dirname(__file__), 'data', 'omte_reference_aug30835.npz')


def test_diamagnetic_vs_neo2_sign_and_order_of_magnitude():
    """Level 0 must match NEO-2 in sign and stay within one order of magnitude."""
    if not os.path.isfile(FIXTURE):
        raise FileNotFoundError(
            f'Reference fixture not found: {FIXTURE}\n'
            'Generate it from a NEO-2 run with isw_calc_Er=1.'
        )

    ref = np.load(FIXTURE)

    ion_idx = np.where(ref['species_tag'] == ref['species_tag_Vphi'])[0][0]
    z_i = ref['species_def'][0, ion_idx, 0]

    # The multispecies output stores scalar species state in n_spec/T_spec.
    # The force-balance estimate must use the actual per-surface profile input.
    assert not np.allclose(
        ref['n_prof'][:, ion_idx], ref['n_spec'][:, ion_idx], rtol=0.0, atol=0.0
    )
    assert not np.allclose(
        ref['T_prof'][:, ion_idx], ref['T_spec'][:, ion_idx], rtol=0.0, atol=0.0
    )

    Om_tE_dia, Er_dia = compute_omte_diamagnetic(
        n=ref['n_prof'][:, ion_idx],
        T=ref['T_prof'][:, ion_idx],
        dn_ds=ref['dn_ov_ds_prof'][:, ion_idx],
        dT_ds=ref['dT_ov_ds_prof'][:, ion_idx],
        z=z_i,
        aiota=ref['aiota'],
        sqrtg_bctrvr_phi=ref['sqrtg_bctrvr_phi'],
        av_nabla_stor=ref['av_nabla_stor'],
    )

    # Regression values from verified computation (AUG #30835)
    # using the actual radial profile input consumed by the run
    Er_expected = np.array([-0.20159416, -0.23357900])
    Om_tE_expected = np.array([-14809.7343, -16736.3417])
    assert np.allclose(Er_dia, Er_expected, rtol=1e-5), (
        f'Er regression: got {Er_dia}, expected {Er_expected}'
    )
    assert np.allclose(Om_tE_dia, Om_tE_expected, rtol=1e-5), (
        f'Om_tE regression: got {Om_tE_dia}, expected {Om_tE_expected}'
    )

    # NEO-2 reference Om_tE (computed from stored Er)
    Om_tE_neo2 = C_CGS * ref['Er_neo2'] / (ref['aiota'] * ref['sqrtg_bctrvr_phi'])

    # Sign must match at every surface
    for i in range(len(Om_tE_neo2)):
        assert np.sign(Om_tE_dia[i]) == np.sign(Om_tE_neo2[i]), (
            f'Sign mismatch at s={ref["boozer_s"][i]:.4f}: '
            f'dia={Om_tE_dia[i]:.1f}, neo2={Om_tE_neo2[i]:.1f}'
        )

    # Diamagnetic alone underestimates because it omits toroidal rotation
    # and the neoclassical transport terms that enter the full force balance.
    ratio = np.abs(Om_tE_dia / Om_tE_neo2)
    for i in range(len(ratio)):
        assert 0.1 < ratio[i] < 1.0, (
            f'Unexpected ratio at s={ref["boozer_s"][i]:.4f}: '
            f'ratio={ratio[i]:.4f}'
        )


def test_toroidal_rotation_proxy_reduces_aug_reference_error():
    """Level 1 proxy should improve the AUG reference fit over Level 0."""
    ref = np.load(FIXTURE)

    ion_idx = np.where(ref['species_tag'] == ref['species_tag_Vphi'])[0][0]
    z_i = ref['species_def'][0, ion_idx, 0]

    om_dia, _ = compute_omte_diamagnetic(
        n=ref['n_prof'][:, ion_idx],
        T=ref['T_prof'][:, ion_idx],
        dn_ds=ref['dn_ov_ds_prof'][:, ion_idx],
        dT_ds=ref['dT_ov_ds_prof'][:, ion_idx],
        z=z_i,
        aiota=ref['aiota'],
        sqrtg_bctrvr_phi=ref['sqrtg_bctrvr_phi'],
        av_nabla_stor=ref['av_nabla_stor'],
    )

    # The HDF5 input stores a toroidal rotation frequency Vphi [rad/s].
    # For the current reference plot we reconstruct a surface-averaged
    # velocity/field pair via v_phi = R0 * Vphi and B_theta ≈ bcovar_tht / R0.
    om_lvl1, er_lvl1 = compute_omte_toroidal_rotation(
        n=ref['n_prof'][:, ion_idx],
        T=ref['T_prof'][:, ion_idx],
        dn_ds=ref['dn_ov_ds_prof'][:, ion_idx],
        dT_ds=ref['dT_ov_ds_prof'][:, ion_idx],
        z=z_i,
        aiota=ref['aiota'],
        sqrtg_bctrvr_phi=ref['sqrtg_bctrvr_phi'],
        av_nabla_stor=ref['av_nabla_stor'],
        v_phi=ref['R0'] * ref['Vphi'],
        b_theta=ref['bcovar_tht'] / ref['R0'],
    )

    er_expected = np.array([-1.06905889, -1.15607731])
    om_expected = np.array([-78536.3915, -82834.9494])
    assert np.allclose(er_lvl1, er_expected, rtol=1e-5), (
        f'Level 1 Er regression: got {er_lvl1}, expected {er_expected}'
    )
    assert np.allclose(om_lvl1, om_expected, rtol=1e-5), (
        f'Level 1 Om_tE regression: got {om_lvl1}, expected {om_expected}'
    )

    om_neo2 = C_CGS * ref['Er_neo2'] / (ref['aiota'] * ref['sqrtg_bctrvr_phi'])
    dia_mae = np.mean(np.abs(om_dia - om_neo2))
    lvl1_mae = np.mean(np.abs(om_lvl1 - om_neo2))
    assert lvl1_mae < dia_mae, (
        f'Expected Level 1 to improve MAE: level0={dia_mae}, level1={lvl1_mae}'
    )


def test_neoclassical_poloidal_proxy_regression():
    """Level 2 banana-regime proxy should be stable on the AUG reference."""
    ref = np.load(FIXTURE)

    ion_idx = np.where(ref['species_tag'] == ref['species_tag_Vphi'])[0][0]
    z_i = ref['species_def'][0, ion_idx, 0]

    om_lvl2, er_lvl2 = compute_omte_neoclassical_poloidal(
        n=ref['n_prof'][:, ion_idx],
        T=ref['T_prof'][:, ion_idx],
        dn_ds=ref['dn_ov_ds_prof'][:, ion_idx],
        dT_ds=ref['dT_ov_ds_prof'][:, ion_idx],
        z=z_i,
        aiota=ref['aiota'],
        sqrtg_bctrvr_phi=ref['sqrtg_bctrvr_phi'],
        av_nabla_stor=ref['av_nabla_stor'],
        v_phi=ref['R0'] * ref['Vphi'],
        b_theta=ref['bcovar_tht'] / ref['R0'],
        b_phi=ref['bcovar_phi'] / ref['R0'],
        k_i=-1.17,
    )

    er_expected = np.array([-1.06905889, -1.15607731])
    om_expected = np.array([-78536.3915, -82834.9494])
    assert np.allclose(er_lvl2, er_expected, rtol=1e-5)
    assert np.allclose(om_lvl2, om_expected, rtol=1e-5)


if __name__ == '__main__':
    test_uniform_profiles_give_zero()
    test_negative_density_gradient_gives_negative_omte()
    test_analytic_linear_profiles()
    test_temperature_gradient_contribution()
    test_force_balance_without_vphi_matches_diamagnetic()
    test_charge_number_scaling()
    test_array_input()
    test_toroidal_rotation_contribution()
    test_neoclassical_poloidal_rotation_formula()
    test_invalid_zero_density_raises()
    test_missing_toroidal_rotation_pair_raises()
    test_missing_poloidal_rotation_pair_raises()
    test_diamagnetic_vs_neo2_sign_and_order_of_magnitude()
    test_toroidal_rotation_proxy_reduces_aug_reference_error()
    test_neoclassical_poloidal_proxy_regression()
    print('\nAll tests passed.')
