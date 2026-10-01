"""Python tools on multi-species inputs with Om_tE only (isw_calc_Er=2).

Oracles are analytic: profiles are cubic polynomials in boozer_s, which a
not-a-knot cubic spline reproduces exactly on any new grid, and the E_r to
Om_tE conversion is checked against the large-aspect-ratio relation
Omega_E = E_r / (R0 B_p) in SI units.
"""
import h5py
import numpy as np
import pytest

from neo2_util.hdf5tools import (add_species_to_profile_file,
                                 change_isw_vphi_loc, new_grid,
                                 remove_species_from_profile_file)
from neo2_ql import get_neo2_ql_input_profiles
from neo2_ql import er_to_om_te, om_te_from_neo2_geometry
from neo2_mars.generate_vrot_from_neo2ql import get_vrot_from_neo2ql

S = np.linspace(0.05, 0.95, 7)
ROT_NAMES = ('Vphi', 'species_tag_Vphi', 'isw_Vphi_loc')


def om_te(s):
    return 1.0e4 * (1.0 - 2.0 * s + 0.5 * s**3)


def vphi(s):
    return 3.0e4 * (1.0 - s**2) + 1.0e3 * s**3


def write_input(path, nspec, with_vphi, with_om_te=True):
    nrad = S.size
    scale = np.arange(1, nspec + 1)[:, None]
    with h5py.File(path, 'w') as f:
        f['num_radial_pts'] = np.array([nrad], dtype=np.int32)
        f['num_species'] = np.array([nspec], dtype=np.int32)
        f['species_tag'] = np.arange(1, nspec + 1, dtype=np.int32)
        species_def = np.zeros((2, nspec, nrad))
        species_def[0] = np.array([-1.0] + [1.0] * (nspec - 1))[:, None]
        species_def[1] = np.array([9.1e-28] + [3.3e-24] * (nspec - 1))[:, None]
        f['species_def'] = species_def
        f['boozer_s'] = S
        f['rho_pol'] = 0.1 + S - 0.2 * S**2
        f['rel_stages'] = np.full(nrad, nspec, dtype=np.int32)
        f['T_prof'] = scale * 1.0e-9 * (2.0 - S**2)
        f['dT_ov_ds_prof'] = scale * 1.0e-9 * (-2.0 * S)
        f['n_prof'] = scale * 1.0e13 * (3.0 - S + S**3)
        f['dn_ov_ds_prof'] = scale * 1.0e13 * (-1.0 + 3.0 * S**2)
        f['kappa_prof'] = scale * 1.0e-5 * (1.0 + S**2)
        if with_om_te:
            f['Om_tE'] = om_te(S)
        if with_vphi:
            f['Vphi'] = vphi(S)
            f['species_tag_Vphi'] = np.array([2], dtype=np.int32)
            f['isw_Vphi_loc'] = np.array([0], dtype=np.int32)


def test_new_grid_interpolates_om_te_on_longer_grid(tmp_path):
    src, out = tmp_path / 'in.h5', tmp_path / 'out.h5'
    write_input(src, 2, with_vphi=False)
    s_new = np.linspace(0.1, 0.9, 13)
    new_grid(str(src), str(out), s_new)
    with h5py.File(out, 'r') as f:
        assert all(name not in f for name in ROT_NAMES)
        assert np.allclose(f['Om_tE'][()], om_te(s_new), rtol=1e-12, atol=0)
        assert np.allclose(f['n_prof'][()][1], 2.0e13 * (3 - s_new + s_new**3),
                           rtol=1e-12, atol=0)
        assert np.allclose(f['dn_ov_ds_prof'][()][0],
                           1.0e13 * (-1 + 3 * s_new**2), rtol=1e-12, atol=0)
        assert np.array_equal(f['boozer_s'][()], s_new)
        assert f['num_radial_pts'][()].tolist() == [13]
        assert f['rel_stages'][()].tolist() == [2] * 13
        assert f['species_def'].shape == (2, 2, 13)


def test_new_grid_interpolates_vphi_and_om_te(tmp_path):
    src, out = tmp_path / 'in.h5', tmp_path / 'out.h5'
    write_input(src, 2, with_vphi=True)
    s_new = np.linspace(0.1, 0.9, 5)
    new_grid(str(src), str(out), s_new)
    with h5py.File(out, 'r') as f:
        assert np.allclose(f['Vphi'][()], vphi(s_new), rtol=1e-12, atol=0)
        assert np.allclose(f['Om_tE'][()], om_te(s_new), rtol=1e-12, atol=0)
        assert f['species_tag_Vphi'][()].tolist() == [2]


@pytest.mark.parametrize('with_vphi', [False, True])
def test_add_and_remove_species_keep_rotation(tmp_path, with_vphi):
    src = tmp_path / 'in.h5'
    write_input(src, 2, with_vphi=with_vphi)
    added = tmp_path / 'added.h5'
    add_species_to_profile_file(str(src), str(added), Zeff=1.5, Ztrace=6.0,
                                mtrace=2.0e-23)
    src4 = tmp_path / 'in4.h5'
    write_input(src4, 4, with_vphi=with_vphi)
    removed = tmp_path / 'removed.h5'
    remove_species_from_profile_file(str(src4), str(removed), 0)
    for name in (added, removed):
        with h5py.File(name, 'r') as f:
            assert np.array_equal(f['Om_tE'][()], om_te(S))
            if with_vphi:
                assert np.array_equal(f['Vphi'][()], vphi(S))
            else:
                assert all(n not in f for n in ROT_NAMES)
    with h5py.File(added, 'r') as f:
        assert f['num_species'][()].tolist() == [3]
        if with_vphi:
            assert f['species_tag_Vphi'][()].tolist() == [2]
    with h5py.File(removed, 'r') as f:
        assert f['num_species'][()].tolist() == [3]
        if with_vphi:  # removing species 0 shifts tag 2 to 1
            assert f['species_tag_Vphi'][()].tolist() == [1]


def test_change_isw_vphi_loc(tmp_path):
    src = tmp_path / 'in.h5'
    write_input(src, 2, with_vphi=False)
    with pytest.raises(ValueError):
        change_isw_vphi_loc(str(src), str(tmp_path / 'x.h5'))
    write_input(src, 2, with_vphi=True)
    change_isw_vphi_loc(str(src), str(tmp_path / 'y.h5'))
    with h5py.File(tmp_path / 'y.h5', 'r') as f:
        assert np.array_equal(f['Om_tE'][()], om_te(S))
        assert int(np.ravel(f['isw_Vphi_loc'][()])[0]) == 2


def test_plot_reader_and_mars_export(tmp_path):
    om_only, both = tmp_path / 'om.h5', tmp_path / 'both.h5'
    write_input(om_only, 2, with_vphi=False)
    write_input(both, 2, with_vphi=True, with_om_te=False)
    profiles, _ = get_neo2_ql_input_profiles(str(om_only))
    assert 'Vphi' not in profiles
    assert np.array_equal(profiles['Om_tE']['y'], om_te(S))
    profiles, _ = get_neo2_ql_input_profiles(str(both))
    assert 'Om_tE' not in profiles
    assert np.array_equal(profiles['Vphi']['y'], vphi(S))
    with pytest.raises(ValueError):
        get_vrot_from_neo2ql(str(om_only))
    _, vrot = get_vrot_from_neo2ql(str(both))
    assert np.array_equal(vrot, vphi(S))


def test_er_to_om_te_large_aspect_ratio():
    # Circular large-aspect-ratio tokamak: dpsi_pol/dr = R0 B_p, hence
    # Omega_E = E_r / (R0 B_p) in SI units.
    r0_m = np.array([1.65, 6.2, 3.0])
    bp_t = np.array([0.4, 1.1, 0.25])
    er = np.array([4.0e3, -1.5e4, 2.5e2])
    sqrtg_bctrvr_tht = (r0_m * 100.0) * (bp_t * 1.0e4)  # G cm
    assert np.allclose(er_to_om_te(er, sqrtg_bctrvr_tht), er / (r0_m * bp_t),
                       rtol=1e-14, atol=0)
    with pytest.raises(ValueError):
        er_to_om_te(1.0, 0.0)


def test_om_te_from_neo2_geometry(tmp_path):
    s_geo = np.array([0.2, 0.5, 0.8])
    geom = np.array([2.0e5, 4.0e5, 5.0e5])
    files = []
    for k in (2, 0, 1):  # order of files must not matter
        name = tmp_path / f'out{k}.h5'
        with h5py.File(name, 'w') as f:
            f['boozer_s'] = s_geo[k]
            f['sqrtg_bctrvr_tht'] = geom[k]
        files.append(str(name))
    s_er = np.linspace(0.0, 1.0, 11)

    def er(s):
        return 1.0e3 * (2.0 - 5.0 * s + s**3)

    result = om_te_from_neo2_geometry(er(s_er), s_er, files, s_geo)
    assert np.allclose(result, er(s_geo) * 1.0e6 / geom, rtol=1e-12, atol=0)
    with pytest.raises(ValueError):
        om_te_from_neo2_geometry(er(s_er), s_er, files, [0.3])
