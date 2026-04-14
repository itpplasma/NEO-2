"""Tests for Om_tE force balance computation."""

import os
import tempfile
from pathlib import Path

import h5py
import numpy as np
from neo2_ql.compute_omte import (
    compute_omte_diamagnetic,
    compute_omte_force_balance,
    compute_omte_neo2_single_ion_limit,
    compute_omte_neoclassical_poloidal,
    compute_omte_neoclassical_poloidal_auto_k,
    compute_omte_toroidal_rotation,
    compute_omte_toroidal_rotation_neo2_convention,
    compute_poloidal_rotation_neoclassical,
    select_poloidal_rotation_coefficient,
    C_CGS,
    E_CGS,
)
from neo2_ql.neo2_output_omte import (
    compute_d31_reference_electron,
    decompose_neo2_er_transport_terms,
    compute_neo2_er_from_transport_coefficients,
    compute_neo2_er_from_k_cof_transport_model,
    compute_neo2_omte_from_transport_coefficients,
    compute_omte_from_neo2_output,
)
from neo2_ql.plot_omte_reference import (
    get_level2_k_scan_reference,
    get_omte_reference_models,
    get_transport_reference_decomposition,
)

FIXTURE_DIR = Path(__file__).with_name('data')
AXISYMMETRIC_OUTPUT_FIXTURE = FIXTURE_DIR / 'neo2_ql_axisymmetric_multispecies_out.h5'


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


def test_poloidal_rotation_pair_product_identity():
    """Physical and Boozer b_phi must give the same v_theta * b_phi product.

    The pair-product identity guarantees that the force balance contribution
    -v_theta * B_phi / c is coordinate-independent: the metric factors in
    v_theta and B_phi cancel exactly.
    """
    dT_ds = -1.0e-9
    z = 1.0
    av_nabla_stor = 0.02
    k_i = -1.17
    r0 = 170.0
    b_phi_phys = -3.0e4
    bcovar_phi = r0 * b_phi_phys

    v_theta_phys = compute_poloidal_rotation_neoclassical(
        dT_ds=dT_ds, z=z, b_phi=b_phi_phys, av_nabla_stor=av_nabla_stor, k_i=k_i,
    )
    v_theta_boozer = compute_poloidal_rotation_neoclassical(
        dT_ds=dT_ds, z=z, b_phi=bcovar_phi, av_nabla_stor=av_nabla_stor, k_i=k_i,
    )
    er_pol_phys = -v_theta_phys * b_phi_phys / C_CGS
    er_pol_boozer = -v_theta_boozer * bcovar_phi / C_CGS
    assert np.isclose(er_pol_phys, er_pol_boozer, rtol=1e-12)
    expected_er_pol = -k_i * dT_ds * av_nabla_stor / (z * E_CGS)
    assert np.isclose(er_pol_phys, expected_er_pol, rtol=1e-12)


def test_toroidal_rotation_neo2_convention_regression():
    """Exact NEO-2 Vphi convention should match the reduced Fortran formula."""
    om_tE, er = compute_omte_toroidal_rotation_neo2_convention(
        n=1e13,
        T=2e-9,
        dn_ds=-1e13,
        dT_ds=0.0,
        z=1.0,
        aiota=0.5,
        sqrtg_bctrvr_phi=1e6,
        av_nabla_stor=0.01,
        vphi=2e5,
        bcovar_tht=-2e3,
        bcovar_phi=-4e4,
    )
    dp_dr = 2e-9 * (-1e13) * 0.01
    pressure = 1e13 * 2e-9
    denom = C_CGS * (-2e3) / 1e6
    nom = 2e5 * (0.5 * (-2e3) - 4e4) + (
        C_CGS * 2e-9 * (-2e3) / (E_CGS * 1e6)
    ) * (dp_dr / pressure)
    er_expected = nom / denom
    om_expected = C_CGS * er_expected / (0.5 * 1e6)
    assert np.isclose(er, er_expected, rtol=1e-12)
    assert np.isclose(om_tE, om_expected, rtol=1e-12)


def test_toroidal_rotation_neo2_convention_without_vphi_matches_diamagnetic():
    """The reduced NEO-2 convention must collapse to Level 0 when Vphi=0."""
    inputs = {
        'n': np.array([2.1e13, 3.2e13]),
        'T': np.array([2.0e-9, 3.0e-9]),
        'dn_ds': np.array([-1.0e13, -2.0e13]),
        'dT_ds': np.array([-1.0e-9, -2.0e-9]),
        'z': 1.0,
        'aiota': np.array([0.4, 0.5]),
        'sqrtg_bctrvr_phi': np.array([5.0e5, 7.0e5]),
        'av_nabla_stor': np.array([0.01, 0.02]),
    }
    om_dia, er_dia = compute_omte_diamagnetic(**inputs)
    om_exact, er_exact = compute_omte_toroidal_rotation_neo2_convention(
        **inputs,
        vphi=np.zeros(2),
        bcovar_tht=np.array([-1.0e3, -2.0e3]),
        bcovar_phi=np.array([-4.0e4, -5.0e4]),
    )
    assert np.allclose(er_exact, er_dia, rtol=1e-12)
    assert np.allclose(om_exact, om_dia, rtol=1e-12)


def test_neoclassical_poloidal_rotation_formula():
    """The Level 2 poloidal estimate should match its analytic definition."""
    v_theta = compute_poloidal_rotation_neoclassical(
        dT_ds=-2e-9,
        z=1.0,
        b_phi=-2e4,
        av_nabla_stor=0.02,
        k_i=-1.17,
    )
    expected = -1.17 * C_CGS * (-2e-9 * 0.02) / (E_CGS * -2e4)
    assert np.isclose(v_theta, expected, rtol=1e-12)


def test_neoclassical_poloidal_rotation_preserves_force_balance_product():
    """Physical and Boozer inputs must give the same v_theta * B_phi product."""
    dT_ds = -2e-9
    z = 1.0
    av_nabla_stor = 0.02
    k_i = -1.17
    r0 = 170.0
    b_phi_phys = -2e4
    bcovar_phi = r0 * b_phi_phys

    v_theta_phys = compute_poloidal_rotation_neoclassical(
        dT_ds=dT_ds,
        z=z,
        b_phi=b_phi_phys,
        av_nabla_stor=av_nabla_stor,
        k_i=k_i,
    )
    v_theta_boozer = compute_poloidal_rotation_neoclassical(
        dT_ds=dT_ds,
        z=z,
        b_phi=bcovar_phi,
        av_nabla_stor=av_nabla_stor,
        k_i=k_i,
    )

    assert np.isclose(v_theta_phys * b_phi_phys, v_theta_boozer * bcovar_phi, rtol=1e-12)


def test_select_poloidal_rotation_coefficient():
    """Auto-K selection should follow the simple regime map."""
    nu_star = np.array([0.01, 1.0, 100.0])
    k_i = select_poloidal_rotation_coefficient(nu_star)
    assert np.allclose(k_i, np.array([-1.17, -0.5, 0.5]))


def test_transport_reconstruction_reduces_to_exact_convention_without_transport_terms():
    """The transport reconstruction must reduce to the reduced Vphi formula."""
    common = {
        'n_spec': np.array([5.0e13, 2.0e13]),
        'T_spec': np.array([8.0e-9, 3.0e-9]),
        'dn_spec_ov_ds': np.array([-1.0e13, -2.0e13]),
        'dT_spec_ov_ds': np.array([-1.0e-9, -3.0e-9]),
        'species_tag': np.array([1, 2]),
        'species_tag_vphi': 2,
        'z_spec': np.array([-1.0, 1.0]),
        'Vphi': 2.0e5,
        'aiota': 0.6,
        'sqrtg_bctrvr_phi': 9.0e5,
        'av_nabla_stor': 0.02,
        'bcovar_tht': -2.0e3,
        'bcovar_phi': -4.0e4,
        'row_ind': np.array([1, 1]),
        'col_ind': np.array([0, 1]),
        'D31_AX': np.zeros(2),
        'D32_AX': np.zeros(2),
        'D33_AX': np.zeros(2),
        'avEparB_ov_avb2': 0.0,
    }
    om_recon, er_recon = compute_neo2_omte_from_transport_coefficients(**common)
    om_exact, er_exact = compute_omte_toroidal_rotation_neo2_convention(
        n=common['n_spec'][1],
        T=common['T_spec'][1],
        dn_ds=common['dn_spec_ov_ds'][1],
        dT_ds=common['dT_spec_ov_ds'][1],
        z=common['z_spec'][1],
        aiota=common['aiota'],
        sqrtg_bctrvr_phi=common['sqrtg_bctrvr_phi'],
        av_nabla_stor=common['av_nabla_stor'],
        vphi=common['Vphi'],
        bcovar_tht=common['bcovar_tht'],
        bcovar_phi=common['bcovar_phi'],
    )
    assert np.isclose(er_recon, er_exact, rtol=1e-12)
    assert np.isclose(om_recon, om_exact, rtol=1e-12)


def test_transport_reconstruction_matches_manual_formula():
    """The D31/D32/D33 reconstruction should match the explicit Fortran sums."""
    kwargs = {
        'n_spec': np.array([4.0e13, 2.5e13]),
        'T_spec': np.array([6.0e-9, 2.0e-9]),
        'dn_spec_ov_ds': np.array([-1.5e13, -2.5e13]),
        'dT_spec_ov_ds': np.array([-1.0e-9, -4.0e-9]),
        'species_tag': np.array([1, 2]),
        'species_tag_vphi': 2,
        'z_spec': np.array([-1.0, 1.0]),
        'Vphi': 3.0e5,
        'aiota': 0.5,
        'sqrtg_bctrvr_phi': 8.0e5,
        'av_nabla_stor': 0.015,
        'bcovar_tht': -1.5e3,
        'bcovar_phi': -3.8e4,
        'row_ind': np.array([1, 1, 0, 0]),
        'col_ind': np.array([0, 1, 0, 1]),
        'D31_AX': np.array([1.2e5, -2.5e5, 4.0e4, -6.0e4]),
        'D32_AX': np.array([2.0e5, 3.0e5, -5.0e4, 8.0e4]),
        'D33_AX': np.array([4.0e4, -3.0e4, 2.0e4, -1.0e4]),
        'avEparB_ov_avb2': 2.5e-7,
    }
    er = compute_neo2_er_from_transport_coefficients(**kwargs)

    spec_i = 1
    z_i = kwargs['z_spec'][spec_i]
    T_i = kwargs['T_spec'][spec_i]
    n_i = kwargs['n_spec'][spec_i]
    dp_dr = (
        T_i * kwargs['dn_spec_ov_ds'][spec_i] + n_i * kwargs['dT_spec_ov_ds'][spec_i]
    ) * kwargs['av_nabla_stor']
    pressure = n_i * T_i
    denom = C_CGS * kwargs['bcovar_tht'] / kwargs['sqrtg_bctrvr_phi']
    nom = (
        kwargs['Vphi'] * (kwargs['aiota'] * kwargs['bcovar_tht'] + kwargs['bcovar_phi'])
        + (C_CGS * T_i * kwargs['bcovar_tht'] / (z_i * E_CGS * kwargs['sqrtg_bctrvr_phi']))
        * (dp_dr / pressure)
    )
    for idx, irow in enumerate(kwargs['row_ind']):
        icol = kwargs['col_ind'][idx]
        if irow == spec_i:
            denom += kwargs['D31_AX'][idx] * (kwargs['z_spec'][icol] * E_CGS) / kwargs['T_spec'][icol]
            nom += kwargs['av_nabla_stor'] * kwargs['D31_AX'][idx] * (
                kwargs['dn_spec_ov_ds'][icol] / kwargs['n_spec'][icol]
                + kwargs['dT_spec_ov_ds'][icol] / kwargs['T_spec'][icol]
            )
            nom += kwargs['av_nabla_stor'] * (
                kwargs['dT_spec_ov_ds'][icol] / kwargs['T_spec'][icol]
            ) * (kwargs['D32_AX'][idx] - 2.5 * kwargs['D31_AX'][idx])
            nom += kwargs['D33_AX'][idx] * kwargs['avEparB_ov_avb2'] * (
                kwargs['z_spec'][icol] * E_CGS
            ) / kwargs['T_spec'][icol]

    assert np.isclose(er, nom / denom, rtol=1e-12)


def test_transport_term_decomposition_matches_manual_formula():
    """The transport decomposition should expose the same exact term sums."""
    kwargs = {
        'n_spec': np.array([4.0e13, 2.5e13]),
        'T_spec': np.array([6.0e-9, 2.0e-9]),
        'dn_spec_ov_ds': np.array([-1.5e13, -2.5e13]),
        'dT_spec_ov_ds': np.array([-1.0e-9, -4.0e-9]),
        'species_tag': np.array([1, 2]),
        'species_tag_vphi': 2,
        'z_spec': np.array([-1.0, 1.0]),
        'Vphi': 3.0e5,
        'aiota': 0.5,
        'sqrtg_bctrvr_phi': 8.0e5,
        'av_nabla_stor': 0.015,
        'bcovar_tht': -1.5e3,
        'bcovar_phi': -3.8e4,
        'row_ind': np.array([1, 1, 0, 0]),
        'col_ind': np.array([0, 1, 0, 1]),
        'D31_AX': np.array([1.2e5, -2.5e5, 4.0e4, -6.0e4]),
        'D32_AX': np.array([2.0e5, 3.0e5, -5.0e4, 8.0e4]),
        'D33_AX': np.array([4.0e4, -3.0e4, 2.0e4, -1.0e4]),
        'avEparB_ov_avb2': 2.5e-7,
    }
    terms = decompose_neo2_er_transport_terms(**kwargs)

    assert np.isclose(terms['nom_dia'], 10532547.533628004, rtol=1e-12)
    assert np.isclose(terms['nom_vphi'], -11625000000.0, rtol=1e-12)
    assert np.isclose(terms['nom_d31'], 10275.0, rtol=1e-12)
    assert np.isclose(terms['nom_d32'], -27500.0, rtol=1e-12)
    assert np.isclose(terms['nom_d33'], -0.0026017356462499997, rtol=1e-12)
    assert np.isclose(terms['denom_base'], -56211085.875, rtol=1e-12)
    assert np.isclose(terms['denom_d31'], -69646.46191499999, rtol=1e-12)
    assert np.isclose(terms['er_total'], 206.3669784526051, rtol=1e-12)
    assert np.isclose(
        terms['er_dia']
        + terms['er_vphi']
        + terms['er_d31']
        + terms['er_d32']
        + terms['er_d33']
        + terms['er_denom'],
        terms['er_total'],
        rtol=1e-12,
    )


def test_d31_reference_matches_stored_normalization_ratio():
    """The electron D31 reference should match the stored raw/normalized ratio."""
    with h5py.File(AXISYMMETRIC_OUTPUT_FIXTURE, 'r') as handle:
        d31_ref = compute_d31_reference_electron(
            T_e=float(np.asarray(handle['T_spec'])[0]),
            z_e=float(np.asarray(handle['z_spec'])[0]),
            aiota=float(np.asarray(handle['aiota'])),
            sqrtg_bctrvr_phi=float(np.asarray(handle['sqrtg_bctrvr_phi'])),
            bcovar_phi=float(np.asarray(handle['bcovar_phi'])),
        )
        stored_ratio = float(np.asarray(handle['D31_AX'])[0] / np.asarray(handle['D31_AX_D31ref'])[0])
    assert np.isclose(d31_ref, stored_ratio, rtol=1e-5)


def test_k_cof_transport_model_zero_d31_matches_exact_convention():
    """Zero D31_hat must collapse to the no-transport exact-convention algebra."""
    common = {
        'n_spec': np.array([5.0e13, 2.0e13]),
        'T_spec': np.array([8.0e-9, 3.0e-9]),
        'dn_spec_ov_ds': np.array([-1.0e13, -2.0e13]),
        'dT_spec_ov_ds': np.array([-1.0e-9, -3.0e-9]),
        'species_tag': np.array([1, 2]),
        'species_tag_vphi': 2,
        'z_spec': np.array([-1.0, 1.0]),
        'Vphi': 2.0e5,
        'aiota': 0.6,
        'sqrtg_bctrvr_phi': 9.0e5,
        'av_nabla_stor': 0.02,
        'bcovar_tht': -2.0e3,
        'bcovar_phi': -4.0e4,
    }
    er_model = compute_neo2_er_from_k_cof_transport_model(
        **common,
        d31_hat=0.0,
        k_cof=0.565,
    )
    _, er_exact = compute_omte_toroidal_rotation_neo2_convention(
        n=common['n_spec'][1],
        T=common['T_spec'][1],
        dn_ds=common['dn_spec_ov_ds'][1],
        dT_ds=common['dT_spec_ov_ds'][1],
        z=common['z_spec'][1],
        aiota=common['aiota'],
        sqrtg_bctrvr_phi=common['sqrtg_bctrvr_phi'],
        av_nabla_stor=common['av_nabla_stor'],
        vphi=common['Vphi'],
        bcovar_tht=common['bcovar_tht'],
        bcovar_phi=common['bcovar_phi'],
    )
    assert np.isclose(er_model, er_exact, rtol=1e-12)


def test_k_cof_transport_model_k_only_matters_with_temperature_gradient():
    """Changing k_cof should have no effect when dT/ds vanishes."""
    common = {
        'n_spec': np.array([5.0e13, 2.0e13]),
        'T_spec': np.array([8.0e-9, 3.0e-9]),
        'dn_spec_ov_ds': np.array([-1.0e13, -2.0e13]),
        'dT_spec_ov_ds': np.zeros(2),
        'species_tag': np.array([1, 2]),
        'species_tag_vphi': 2,
        'z_spec': np.array([-1.0, 1.0]),
        'Vphi': 2.0e5,
        'aiota': 0.6,
        'sqrtg_bctrvr_phi': 9.0e5,
        'av_nabla_stor': 0.02,
        'bcovar_tht': -2.0e3,
        'bcovar_phi': -4.0e4,
        'd31_hat': -1.0,
    }
    er_low_k = compute_neo2_er_from_k_cof_transport_model(**common, k_cof=0.2)
    er_high_k = compute_neo2_er_from_k_cof_transport_model(**common, k_cof=1.4)
    assert np.isclose(er_low_k, er_high_k, rtol=1e-12)


def test_k_cof_transport_model_axisymmetric_fixture_regression():
    """The minimal ion-ion transport model should strongly reduce the no-transport overshoot."""
    with h5py.File(AXISYMMETRIC_OUTPUT_FIXTURE, 'r') as handle:
        species_tag = np.asarray(handle['species_tag'])
        species_tag_vphi = int(np.asarray(handle['species_tag_Vphi']).reshape(-1)[0])
        ion_index = np.where(species_tag == species_tag_vphi)[0][0]
        d31_hat = float(np.asarray(handle['D31_AX_D31ref'])[ion_index * species_tag.size + ion_index])
        d32_hat = float(np.asarray(handle['D32_AX_D31ref'])[ion_index * species_tag.size + ion_index])
        k_cof = 2.5 - d32_hat / d31_hat

        er_model = compute_neo2_er_from_k_cof_transport_model(
            n_spec=np.asarray(handle['n_spec']),
            T_spec=np.asarray(handle['T_spec']),
            dn_spec_ov_ds=np.asarray(handle['dn_spec_ov_ds']),
            dT_spec_ov_ds=np.asarray(handle['dT_spec_ov_ds']),
            species_tag=species_tag,
            species_tag_vphi=species_tag_vphi,
            z_spec=np.asarray(handle['z_spec']),
            Vphi=float(np.asarray(handle['Vphi'])),
            aiota=float(np.asarray(handle['aiota'])),
            sqrtg_bctrvr_phi=float(np.asarray(handle['sqrtg_bctrvr_phi'])),
            av_nabla_stor=float(np.asarray(handle['av_nabla_stor'])),
            bcovar_tht=float(np.asarray(handle['bcovar_tht'])),
            bcovar_phi=float(np.asarray(handle['bcovar_phi'])),
            d31_hat=d31_hat,
            k_cof=k_cof,
        )
        er_stored = float(np.asarray(handle['Er']))

    _, er_no_transport = compute_omte_toroidal_rotation_neo2_convention(
        n=float(np.asarray([2.13663262e+13, 2.13663262e+13])[1]),
        T=float(np.asarray([2.49794945e-09, 2.54978413e-09])[1]),
        dn_ds=float(np.asarray([-2.18434868e+13, -2.18434868e+13])[1]),
        dT_ds=float(np.asarray([-3.10182216e-09, -1.81389942e-09])[1]),
        z=1.0,
        aiota=0.46411481338020744,
        sqrtg_bctrvr_phi=898569.3223996271,
        av_nabla_stor=0.02503408702076222,
        vphi=22397.03710787831,
        bcovar_tht=-122680.92474705892,
        bcovar_phi=-2928048.8002521824,
    )
    assert np.isclose(er_model, 0.13877985508811017, rtol=1e-12)
    assert abs(er_model - er_stored) < abs(er_no_transport - er_stored)


def test_compute_omte_from_output_reads_stored_er():
    """Stored-output mode should return the exact stored Er curve."""
    with tempfile.NamedTemporaryFile(suffix='.h5') as tmp:
        with h5py.File(tmp.name, 'w') as handle:
            handle.create_dataset('Er', data=np.array([-0.2, -0.5]))
            handle.create_dataset('aiota', data=np.array([0.4, 0.5]))
            handle.create_dataset('sqrtg_bctrvr_phi', data=np.array([5.0e5, 7.0e5]))
        om_tE, er = compute_omte_from_neo2_output(tmp.name, mode='stored')
    assert np.allclose(er, np.array([-0.2, -0.5]))
    assert np.allclose(
        om_tE,
        C_CGS * np.array([-0.2, -0.5]) / (np.array([0.4, 0.5]) * np.array([5.0e5, 7.0e5])),
    )


def test_compute_omte_from_output_reconstructs_transport_mode():
    """Transport-output mode should reproduce the reconstructed Er curve."""
    species_tag = np.array([1, 2], dtype=np.int32)
    species_tag_vphi = np.array([2], dtype=np.int32)
    species_def = np.array(
        [
            [[-1.0], [1.0]],
            [[9.1e-28], [3.34e-24]],
        ]
    )
    kwargs = {
        'n_spec': np.array([4.0e13, 2.5e13]),
        'T_spec': np.array([6.0e-9, 2.0e-9]),
        'dn_spec_ov_ds': np.array([-1.5e13, -2.5e13]),
        'dT_spec_ov_ds': np.array([-1.0e-9, -4.0e-9]),
        'species_tag': species_tag,
        'species_tag_vphi': 2,
        'z_spec': np.array([-1.0, 1.0]),
        'Vphi': 3.0e5,
        'aiota': 0.5,
        'sqrtg_bctrvr_phi': 8.0e5,
        'av_nabla_stor': 0.015,
        'bcovar_tht': -1.5e3,
        'bcovar_phi': -3.8e4,
        'row_ind': np.array([1, 1, 0, 0]),
        'col_ind': np.array([1, 2, 1, 2]),
        'D31_AX': np.array([1.2e5, -2.5e5, 4.0e4, -6.0e4]),
        'D32_AX': np.array([2.0e5, 3.0e5, -5.0e4, 8.0e4]),
        'D33_AX': np.array([4.0e4, -3.0e4, 2.0e4, -1.0e4]),
        'avEparB_ov_avb2': 2.5e-7,
    }
    er_expected = compute_neo2_er_from_transport_coefficients(**kwargs)

    with tempfile.NamedTemporaryFile(suffix='.h5') as tmp:
        with h5py.File(tmp.name, 'w') as handle:
            handle.create_dataset('species_tag', data=species_tag)
            handle.create_dataset('species_tag_Vphi', data=species_tag_vphi)
            handle.create_dataset('species_def', data=species_def)
            handle.create_dataset('n_spec', data=np.array([kwargs['n_spec']]))
            handle.create_dataset('T_spec', data=np.array([kwargs['T_spec']]))
            handle.create_dataset('dn_ov_ds_prof', data=np.array([kwargs['dn_spec_ov_ds']]))
            handle.create_dataset('dT_ov_ds_prof', data=np.array([kwargs['dT_spec_ov_ds']]))
            handle.create_dataset('aiota', data=np.array([kwargs['aiota']]))
            handle.create_dataset('sqrtg_bctrvr_phi', data=np.array([kwargs['sqrtg_bctrvr_phi']]))
            handle.create_dataset('av_nabla_stor', data=np.array([kwargs['av_nabla_stor']]))
            handle.create_dataset('bcovar_tht', data=np.array([kwargs['bcovar_tht']]))
            handle.create_dataset('bcovar_phi', data=np.array([kwargs['bcovar_phi']]))
            handle.create_dataset('Vphi', data=np.array([kwargs['Vphi']]))
            handle.create_dataset('row_ind_spec', data=np.array([kwargs['row_ind']]))
            handle.create_dataset('col_ind_spec', data=np.array([kwargs['col_ind']]))
            handle.create_dataset('D31_AX', data=np.array([kwargs['D31_AX']]))
            handle.create_dataset('D32_AX', data=np.array([kwargs['D32_AX']]))
            handle.create_dataset('D33_AX', data=np.array([kwargs['D33_AX']]))
            handle.create_dataset('avEparB_ov_avb2', data=np.array([kwargs['avEparB_ov_avb2']]))
            handle.create_dataset('isw_Vphi_loc', data=np.array([0], dtype=np.int32))
        om_tE, er = compute_omte_from_neo2_output(tmp.name, mode='transport')

    assert np.allclose(er, np.array([er_expected]))
    assert np.allclose(
        om_tE,
        np.array([C_CGS * er_expected / (kwargs['aiota'] * kwargs['sqrtg_bctrvr_phi'])]),
    )


def test_compute_omte_from_output_reads_summary_group():
    """Summary files should be readable via the neo2_multispecies_out group."""
    with tempfile.NamedTemporaryFile(suffix='.h5') as tmp:
        with h5py.File(tmp.name, 'w') as handle:
            group = handle.create_group('neo2_multispecies_out')
            group.create_dataset('Er', data=np.array([-0.2, -0.5]))
            group.create_dataset('aiota', data=np.array([0.4, 0.5]))
            group.create_dataset('sqrtg_bctrvr_phi', data=np.array([5.0e5, 7.0e5]))
        om_tE, er = compute_omte_from_neo2_output(tmp.name, mode='stored')
    assert np.allclose(er, np.array([-0.2, -0.5]))
    assert np.allclose(
        om_tE,
        C_CGS * np.array([-0.2, -0.5]) / (np.array([0.4, 0.5]) * np.array([5.0e5, 7.0e5])),
    )


def test_axisymmetric_output_fixture_reconstructs_full_er():
    """Real NEO-2 output should replay the full compute_Er algebra."""
    assert AXISYMMETRIC_OUTPUT_FIXTURE.exists()
    with h5py.File(AXISYMMETRIC_OUTPUT_FIXTURE, 'r') as handle:
        for name in (
            'species_tag_Vphi',
            'isw_Vphi_loc',
            'Vphi',
            'dn_spec_ov_ds',
            'dT_spec_ov_ds',
            'av_nabla_stor',
        ):
            assert name in handle

    om_stored, er_stored = compute_omte_from_neo2_output(
        AXISYMMETRIC_OUTPUT_FIXTURE, mode='stored'
    )
    om_transport, er_transport = compute_omte_from_neo2_output(
        AXISYMMETRIC_OUTPUT_FIXTURE, mode='transport'
    )

    er_stored = np.asarray(er_stored).reshape(-1)
    om_stored = np.asarray(om_stored).reshape(-1)
    er_transport = np.asarray(er_transport).reshape(-1)
    om_transport = np.asarray(om_transport).reshape(-1)

    assert np.allclose(er_stored, np.array([0.11359881049568062]), rtol=0.0, atol=1e-15)
    assert np.allclose(om_stored, np.array([8166.1521825278678]), rtol=0.0, atol=1e-10)
    assert np.allclose(er_transport, np.array([0.11359863968963073]), rtol=0.0, atol=1e-15)
    assert np.allclose(om_transport, np.array([8166.1399039820699]), rtol=0.0, atol=1e-10)

    assert np.allclose(er_transport, er_stored, rtol=2e-6, atol=1e-9)
    assert np.allclose(om_transport, om_stored, rtol=2e-6, atol=1e-5)


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


def test_exact_neo2_convention_zero_pressure_raises():
    """Exact NEO-2 Vphi path should fail on zero pressure."""
    try:
        compute_omte_toroidal_rotation_neo2_convention(
            n=1e13,
            T=0.0,
            dn_ds=-1e13,
            dT_ds=0.0,
            z=1.0,
            aiota=0.5,
            sqrtg_bctrvr_phi=1e6,
            av_nabla_stor=0.01,
            vphi=1e5,
            bcovar_tht=-2e3,
            bcovar_phi=-4e4,
        )
    except ValueError as exc:
        assert str(exc) == 'n*T must be nonzero to compute NEO-2 Vphi model'
    else:
        raise AssertionError('Expected ValueError for zero pressure')


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

    # The rebuilt fixture now stores the exact run-local multispecies inputs
    # from neo2.in, which should be consistent with the scalar state in the
    # multispecies output.
    assert np.allclose(
        ref['n_prof'][:, ion_idx], ref['n_spec'][:, ion_idx], rtol=0.0, atol=0.0
    )
    assert np.allclose(
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
    Er_expected = np.array([-2.93232045, -4.35885986])
    Om_tE_expected = np.array([-215417.3806, -312319.8873])
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

    # Diamagnetic alone gives a fraction of the full neoclassical result
    # because it omits toroidal rotation and the neoclassical transport
    # terms (D31, D32, D33) that enter the full force balance.
    ratio = np.abs(Om_tE_dia / Om_tE_neo2)
    for i in range(len(ratio)):
        assert 1.0 < ratio[i] < 10.0, (
            f'Unexpected ratio at s={ref["boozer_s"][i]:.4f}: '
            f'ratio={ratio[i]:.4f}'
        )


def test_toroidal_rotation_from_neo2_component_pair_regression():
    """Level 1 should accept the NEO-2 contravariant/covariant pair directly."""
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

    om_lvl1, er_lvl1 = compute_omte_toroidal_rotation(
        n=ref['n_prof'][:, ion_idx],
        T=ref['T_prof'][:, ion_idx],
        dn_ds=ref['dn_ov_ds_prof'][:, ion_idx],
        dT_ds=ref['dT_ov_ds_prof'][:, ion_idx],
        z=z_i,
        aiota=ref['aiota'],
        sqrtg_bctrvr_phi=ref['sqrtg_bctrvr_phi'],
        av_nabla_stor=ref['av_nabla_stor'],
        v_phi=ref['Vphi'],
        b_theta=ref['bcovar_tht'],
    )

    er_expected = np.array([-3.14157553, -4.67126557])
    om_expected = np.array([-230789.9101, -334704.2995])
    assert np.allclose(er_lvl1, er_expected, rtol=1e-5), (
        f'Level 1 Er regression: got {er_lvl1}, expected {er_expected}'
    )
    assert np.allclose(om_lvl1, om_expected, rtol=1e-5), (
        f'Level 1 Om_tE regression: got {om_lvl1}, expected {om_expected}'
    )

    assert np.all(np.sign(om_dia) == np.sign(om_lvl1))


def test_aug_reference_reduced_neo2_vphi_convention_regression():
    """The reduced isw_Vphi_loc=0 algebra remains available as a separate mode."""
    ref = np.load(FIXTURE)
    ion_idx = np.where(ref['species_tag'] == ref['species_tag_Vphi'])[0][0]
    z_i = ref['species_def'][0, ion_idx, 0]
    om_exact, er_exact = compute_omte_toroidal_rotation_neo2_convention(
        n=ref['n_prof'][:, ion_idx],
        T=ref['T_prof'][:, ion_idx],
        dn_ds=ref['dn_ov_ds_prof'][:, ion_idx],
        dT_ds=ref['dT_ov_ds_prof'][:, ion_idx],
        z=z_i,
        aiota=ref['aiota'],
        sqrtg_bctrvr_phi=ref['sqrtg_bctrvr_phi'],
        av_nabla_stor=ref['av_nabla_stor'],
        vphi=ref['Vphi'],
        bcovar_tht=ref['bcovar_tht'],
        bcovar_phi=ref['bcovar_phi'],
    )
    er_expected = np.array([53.22434201, 51.62016311])
    om_expected = np.array([3910025.7094, 3698674.4302])
    assert np.allclose(er_exact, er_expected, rtol=1e-5)
    assert np.allclose(om_exact, om_expected, rtol=1e-5)


def test_aug_reference_single_ion_limit_bridges_simple_and_full_neo2():
    """The reduced single-ion NEO-2 limit should explain most of the AUG gap."""
    ref = np.load(FIXTURE)
    ion_idx = np.where(ref['species_tag'] == ref['species_tag_Vphi'])[0][0]
    z_i = ref['species_def'][0, ion_idx, 0]
    d31 = ref['D31_AX'][:, 3]
    d32 = ref['D32_AX'][:, 3]
    k_ii = 2.5 - d32 / d31

    _, er_simple = compute_omte_neoclassical_poloidal(
        n=ref['n_prof'][:, ion_idx],
        T=ref['T_prof'][:, ion_idx],
        dn_ds=ref['dn_ov_ds_prof'][:, ion_idx],
        dT_ds=ref['dT_ov_ds_prof'][:, ion_idx],
        z=z_i,
        aiota=ref['aiota'],
        sqrtg_bctrvr_phi=ref['sqrtg_bctrvr_phi'],
        av_nabla_stor=ref['av_nabla_stor'],
        v_phi=ref['Vphi'],
        b_theta=ref['bcovar_tht'],
        b_phi=ref['bcovar_phi'],
        k_i=k_ii,
    )
    om_reduced, er_reduced = compute_omte_neo2_single_ion_limit(
        n=ref['n_prof'][:, ion_idx],
        T=ref['T_prof'][:, ion_idx],
        dn_ds=ref['dn_ov_ds_prof'][:, ion_idx],
        dT_ds=ref['dT_ov_ds_prof'][:, ion_idx],
        z=z_i,
        aiota=ref['aiota'],
        sqrtg_bctrvr_phi=ref['sqrtg_bctrvr_phi'],
        av_nabla_stor=ref['av_nabla_stor'],
        v_phi=ref['Vphi'],
        b_theta=ref['bcovar_tht'],
        b_phi=ref['bcovar_phi'],
        k_i=k_ii,
    )

    er_reduced_expected = np.array([-0.68487916, -1.69680401])
    om_reduced_expected = np.array([-50313.3532212, -121578.95746675])
    assert np.allclose(er_reduced, er_reduced_expected, rtol=1e-5)
    assert np.allclose(om_reduced, om_reduced_expected, rtol=1e-5)

    er_stored = ref['Er_neo2']
    assert np.all(np.abs(er_reduced - er_stored) < np.abs(er_simple - er_stored))
    assert np.all(np.sign(er_reduced) == np.sign(er_stored))


def test_aug_reference_mars_banana_shortcut_is_closer_by_compensating_errors():
    """The MARS k=1.17 shortcut can look closer for the wrong reason."""
    ref = np.load(FIXTURE)
    ion_idx = np.where(ref['species_tag'] == ref['species_tag_Vphi'])[0][0]
    z_i = ref['species_def'][0, ion_idx, 0]

    _, er_mars = compute_omte_neoclassical_poloidal(
        n=ref['n_prof'][:, ion_idx],
        T=ref['T_prof'][:, ion_idx],
        dn_ds=ref['dn_ov_ds_prof'][:, ion_idx],
        dT_ds=ref['dT_ov_ds_prof'][:, ion_idx],
        z=z_i,
        aiota=ref['aiota'],
        sqrtg_bctrvr_phi=ref['sqrtg_bctrvr_phi'],
        av_nabla_stor=ref['av_nabla_stor'],
        v_phi=ref['Vphi'],
        b_theta=ref['bcovar_tht'],
        b_phi=ref['bcovar_phi'],
        k_i=1.17,
    )
    _, er_reduced_banana = compute_omte_neo2_single_ion_limit(
        n=ref['n_prof'][:, ion_idx],
        T=ref['T_prof'][:, ion_idx],
        dn_ds=ref['dn_ov_ds_prof'][:, ion_idx],
        dT_ds=ref['dT_ov_ds_prof'][:, ion_idx],
        z=z_i,
        aiota=ref['aiota'],
        sqrtg_bctrvr_phi=ref['sqrtg_bctrvr_phi'],
        av_nabla_stor=ref['av_nabla_stor'],
        v_phi=ref['Vphi'],
        b_theta=ref['bcovar_tht'],
        b_phi=ref['bcovar_phi'],
        k_i=1.17,
    )

    assert np.allclose(er_mars, np.array([-0.83550899, -1.24332478]), rtol=1e-5)
    assert np.allclose(er_reduced_banana, np.array([0.37948147, 0.07855097]), rtol=1e-5)
    assert np.all(np.sign(er_mars) == np.sign(ref['Er_neo2']))
    assert np.all(np.sign(er_reduced_banana) != np.sign(ref['Er_neo2']))


def test_aug_reference_transport_replay_matches_stored_er():
    """The rebuilt AUG fixture must replay the stored transport solution."""
    models = get_omte_reference_models(FIXTURE)
    assert np.allclose(models['om_transport_replay'], models['om_neo2'], rtol=0.0, atol=1e-1)
    assert np.allclose(models['er_transport_replay'], models['er_neo2'], rtol=0.0, atol=2e-6)


def test_aug_reference_fixture_contains_full_geometry_fields():
    """The AUG reference fixture should carry the extra geometry from source HDF5 files."""
    ref = np.load(FIXTURE)
    for name in [
        'avbhat2',
        'avb2',
        'av_inv_bhat',
        'av_gphph',
        'D31ref0',
        'avEparB_ov_avb2',
        'D31_AX',
        'D32_AX',
        'D33_AX',
        'D31_AX_D31ref',
        'D32_AX_D31ref',
        'D33_AX_norm',
        'row_ind_spec',
        'col_ind_spec',
        'source_h5_paths',
    ]:
        assert name in ref.files

    assert np.allclose(ref['avb2'], ref['avbhat2'] * ref['Bref'] ** 2, rtol=1e-12)
    assert np.allclose(ref['Er_neo2'], ref['Er'], rtol=0.0, atol=0.0)
    assert np.allclose(ref['sqrtg_bctrvr_phi'] / (1.0 / ref['aiota']), ref['sqrtg_bctrvr_tht'], rtol=1e-12)


def test_transport_plot_decomposition_regression():
    """The transport plot helper should expose the stored exact term breakdown."""
    terms = get_transport_reference_decomposition()
    assert np.isclose(terms['er_dia'], -0.23040107384987163, rtol=1e-12)
    assert np.isclose(terms['er_vphi'], 16.333784941637138, rtol=1e-12)
    assert np.isclose(terms['er_d31'], -11.40113413101951, rtol=1e-12)
    assert np.isclose(terms['er_d32'], 2.537576526707158, rtol=1e-12)
    assert np.isclose(terms['er_d33'], -1.2868034606727117, rtol=1e-12)
    assert np.isclose(terms['er_denom'], -5.839424163112569, rtol=1e-12)
    assert np.isclose(terms['er_total'], 0.11359863968963077, rtol=1e-12)
    assert np.isclose(terms['er_stored'], 0.11359881049568062, rtol=0.0, atol=1e-15)


if __name__ == '__main__':
    test_uniform_profiles_give_zero()
    test_negative_density_gradient_gives_negative_omte()
    test_analytic_linear_profiles()
    test_temperature_gradient_contribution()
    test_force_balance_without_vphi_matches_diamagnetic()
    test_charge_number_scaling()
    test_array_input()
    test_toroidal_rotation_contribution()
    test_toroidal_rotation_neo2_convention_regression()
    test_toroidal_rotation_neo2_convention_without_vphi_matches_diamagnetic()
    test_neoclassical_poloidal_rotation_formula()
    test_neoclassical_poloidal_rotation_preserves_force_balance_product()
    test_select_poloidal_rotation_coefficient()
    test_transport_reconstruction_reduces_to_exact_convention_without_transport_terms()
    test_transport_reconstruction_matches_manual_formula()
    test_compute_omte_from_output_reads_stored_er()
    test_compute_omte_from_output_reconstructs_transport_mode()
    test_invalid_zero_density_raises()
    test_exact_neo2_convention_zero_pressure_raises()
    test_missing_toroidal_rotation_pair_raises()
    test_missing_poloidal_rotation_pair_raises()
    test_diamagnetic_vs_neo2_sign_and_order_of_magnitude()
    test_toroidal_rotation_from_neo2_component_pair_reduces_aug_reference_error()
    test_neoclassical_poloidal_with_neo2_component_pair_regression()
    test_neoclassical_poloidal_auto_k_with_neo2_component_pair_regression()
    test_aug_reference_reduced_neo2_vphi_convention_regression()
    test_aug_reference_single_ion_limit_bridges_simple_and_full_neo2()
    test_aug_reference_mars_banana_shortcut_is_closer_by_compensating_errors()
    test_reference_plot_models_regression()
    print('\nAll tests passed.')
