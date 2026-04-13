"""Tests for Om_tE force balance computation."""

import os
import numpy as np
from neo2_ql.compute_omte import compute_omte_diamagnetic, C_CGS, E_CGS


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


# --- E2E test against NEO-2 reference ---

FIXTURE = os.path.join(os.path.dirname(__file__), 'data', 'omte_reference_aug30835.npz')


def test_diamagnetic_vs_neo2_sign_and_order_of_magnitude():
    """Level 0 must match NEO-2 in sign and be within one order of magnitude."""
    if not os.path.isfile(FIXTURE):
        raise FileNotFoundError(
            f'Reference fixture not found: {FIXTURE}\n'
            'Generate it from a NEO-2 run with isw_calc_Er=1.'
        )

    ref = np.load(FIXTURE)

    # Ion species index (species_tag_Vphi points to the ion species)
    isp = int(ref['species_tag_Vphi']) - 1
    z_i = ref['z_spec'][0, abs(int(ref['z_spec'][0].argmax()))]  # positive charge species
    # Find the ion species (z > 0)
    z_spec = ref['z_spec'][0]  # same for all surfaces
    ion_idx = np.where(z_spec > 0)[0][0]
    z_i = z_spec[ion_idx]

    Om_tE_dia, Er_dia = compute_omte_diamagnetic(
        n=ref['n_prof'][ion_idx],
        T=ref['T_prof'][ion_idx],
        dn_ds=ref['dn_ov_ds_prof'][ion_idx],
        dT_ds=ref['dT_ov_ds_prof'][ion_idx],
        z=z_i,
        aiota=ref['aiota'],
        sqrtg_bctrvr_phi=ref['sqrtg_bctrvr_phi'],
        av_nabla_stor=ref['av_nabla_stor'],
    )

    # NEO-2 reference Om_tE (computed from stored Er)
    Om_tE_neo2 = C_CGS * ref['Er_neo2'] / (ref['aiota'] * ref['sqrtg_bctrvr_phi'])

    # Sign must match at every surface
    for i in range(len(Om_tE_neo2)):
        assert np.sign(Om_tE_dia[i]) == np.sign(Om_tE_neo2[i]), (
            f'Sign mismatch at s={ref["boozer_s"][i]:.4f}: '
            f'dia={Om_tE_dia[i]:.1f}, neo2={Om_tE_neo2[i]:.1f}'
        )

    # Order of magnitude: ratio must be between 0.01 and 100
    ratio = np.abs(Om_tE_dia / Om_tE_neo2)
    for i in range(len(ratio)):
        assert 0.01 < ratio[i] < 100, (
            f'Order of magnitude mismatch at s={ref["boozer_s"][i]:.4f}: '
            f'ratio={ratio[i]:.4f}'
        )

    # Diamagnetic alone should underestimate (ratio < 1) because it
    # misses rotation and neoclassical terms.
    for i in range(len(ratio)):
        assert ratio[i] < 1.0, (
            f'Diamagnetic exceeds full neoclassical at s={ref["boozer_s"][i]:.4f}: '
            f'ratio={ratio[i]:.4f} (expected < 1)'
        )

    # Print comparison for manual inspection
    print('\nOm_tE comparison (Level 0 diamagnetic vs NEO-2):')
    print(f'{"s_tor":>8s}  {"NEO-2 [krad/s]":>16s}  {"Level 0 [krad/s]":>16s}  {"ratio":>8s}')
    for i in range(len(Om_tE_neo2)):
        print(f'{ref["boozer_s"][i]:8.4f}  {Om_tE_neo2[i]/1e3:16.2f}  '
              f'{Om_tE_dia[i]/1e3:16.2f}  {ratio[i]:8.4f}')


if __name__ == '__main__':
    test_uniform_profiles_give_zero()
    test_negative_density_gradient_gives_negative_omte()
    test_analytic_linear_profiles()
    test_temperature_gradient_contribution()
    test_charge_number_scaling()
    test_array_input()
    test_diamagnetic_vs_neo2_sign_and_order_of_magnitude()
    print('\nAll tests passed.')
