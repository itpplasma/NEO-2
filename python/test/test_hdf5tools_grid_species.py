"""new_grid rel_stages and remove_species_from_profile_file with 2 species left.

Oracles:
- neo2.f90 (prepare_mulitspecies_scan) accepts a surface only if rel_stages
  equals the number of species with n_prof > 0 there; every regridded
  surface must satisfy that, and a surface between two old surfaces must not
  get fewer species than both of them.
- The ion collisionality recomputed after removing a species must equal the
  value generate_multispec_input's formula gives for the same ion density,
  temperature, charge and electron Coulomb logarithm, and the electron row
  must be untouched.
"""
import h5py
import numpy as np
import pytest

from neo2_ql import get_coulomb_logarithm, get_kappa
from neo2_util.hdf5tools import new_grid, remove_species_from_profile_file

E_CGS = 1.60217662e-19 * 2.99792458e8 * 10
S_OLD = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
# Shifted surfaces, same count (new_grid writes datasets in place).
S_NEW = np.array([0.05, 0.3, 0.5, 0.55, 0.7, 0.95])


def write_input(path, n_prof, rel_stages, t_scale=None):
    nspec, nrad = n_prof.shape
    t_scale = np.ones(nspec) if t_scale is None else np.asarray(t_scale)
    t_prof = t_scale[:, None] * 1.0e-9 * (2.0 - S_OLD)
    species_def = np.zeros((2, nspec, nrad))
    species_def[0] = np.array([-1.0, 1.0, 6.0, 2.0][:nspec])[:, None]
    species_def[1] = np.array([9.1e-28, 3.3e-24, 2.0e-23, 6.6e-24][:nspec])[:, None]
    charge = species_def[0] * E_CGS
    log_lambda = get_coulomb_logarithm(n_prof[0], t_prof[0])
    with h5py.File(path, 'w') as f:
        f['num_radial_pts'] = np.array([nrad], dtype=np.int32)
        f['num_species'] = np.array([nspec], dtype=np.int32)
        f['species_tag'] = np.arange(1, nspec + 1, dtype=np.int32)
        f['species_def'] = species_def
        f['boozer_s'] = S_OLD
        f['rho_pol'] = np.sqrt(S_OLD)
        f['rel_stages'] = np.asarray(rel_stages, dtype=np.int32)
        f['n_prof'] = n_prof
        f['dn_ov_ds_prof'] = np.gradient(n_prof, S_OLD, axis=1)
        f['T_prof'] = t_prof
        f['dT_ov_ds_prof'] = np.gradient(t_prof, S_OLD, axis=1)
        with np.errstate(divide='ignore'):
            kappa = get_kappa(n_prof, t_prof, charge, log_lambda)
        f['kappa_prof'] = np.where(n_prof > 0, kappa, 0.0)  # absent species
        f['Vphi'] = 1.0e4 * (1.0 - S_OLD)
        f['species_tag_Vphi'] = np.array([2], dtype=np.int32)
        f['isw_Vphi_loc'] = np.array([0], dtype=np.int32)


def carbon_edge_profiles(s_old):
    # Carbon (species 3) is present inside s = 0.4 and absent outside;
    # deuterium keeps the plasma quasi-neutral (n_e = n_D + 6 n_C).
    n_e = 4.0e13 * (1.2 - s_old**2)
    n_c = 1.0e12 * np.array([2.0, 1.8, 1.0, 0.0, 0.0, 0.0])
    return np.vstack([n_e, n_e - 6.0 * n_c, n_c])


def test_rel_stages_follow_species_with_density(tmp_path):
    n_prof = carbon_edge_profiles(S_OLD)
    rel_old = (n_prof > 0).sum(axis=0)
    assert rel_old.tolist() == [3, 3, 3, 2, 2, 2]
    src, out = tmp_path / 'in.h5', tmp_path / 'out.h5'
    write_input(src, n_prof, rel_old)
    s_new = np.array([0.05, 0.3, 0.5, 0.55, 0.8, 1.0])
    new_grid(str(src), str(out), s_new)
    with h5py.File(out, 'r') as f:
        rel_new = f['rel_stages'][()]
        n_new = f['n_prof'][()]
        z = f['species_def'][0][:, 0]
        assert f['rel_stages'].dtype == np.int32
    # Contract of the Fortran reader on every surface.
    assert rel_new.tolist() == (n_new > 0).sum(axis=0).tolist()
    # Never fewer species than both neighbouring old surfaces.
    for k, s in enumerate(s_new):
        right = min(np.searchsorted(S_OLD, s), S_OLD.size - 1)
        left = right if np.isclose(S_OLD[right], s) else right - 1
        assert rel_new[k] >= min(rel_old[left], rel_old[right]), (s, rel_new)
    assert rel_new.tolist() == [3, 3, 3, 3, 2, 2]
    # Densities are not altered: quasi-neutrality holds on every surface.
    assert np.allclose(z @ n_new, 0.0, rtol=0, atol=1e-12 * n_new[0].max())


def test_overshoot_into_absent_region_is_refused(tmp_path):
    # Between s = 0.8 and 1.0 carbon is absent on both surfaces, but its
    # spline overshoots above zero at 0.95. Existing output stays untouched.
    src, out = tmp_path / 'in.h5', tmp_path / 'out.h5'
    n_prof = carbon_edge_profiles(S_OLD)
    write_input(src, n_prof, (n_prof > 0).sum(axis=0))
    out.write_bytes(b'keep')
    with pytest.raises(ValueError, match=r'neighbouring.*0\.95'):
        new_grid(str(src), str(out),
                 np.array([0.05, 0.3, 0.5, 0.55, 0.8, 0.95]))
    assert out.read_bytes() == b'keep'


def test_negative_density_is_refused(tmp_path):
    # Between s = 0.6 and 0.8 the carbon spline dips below zero (s = 0.7);
    # neo2.f90 would drop carbon there and lose quasi-neutrality.
    src = tmp_path / 'in.h5'
    n_prof = carbon_edge_profiles(S_OLD)
    write_input(src, n_prof, (n_prof > 0).sum(axis=0))
    with pytest.raises(ValueError, match='negative'):
        new_grid(str(src), str(tmp_path / 'out.h5'),
                 np.array([0.05, 0.3, 0.5, 0.55, 0.7, 0.8]))


def test_identity_grid_keeps_zero_density_at_knot(tmp_path):
    # Spline roundoff at a knot with zero density must not add a species.
    n_e = np.full(S_OLD.size, 4.0e13)
    n_c = 1.0e11 * np.array([1.0, 1.0, 1.0, 10.0, 1.0, 0.0])
    n_prof = np.vstack([n_e, n_e - 6.0 * n_c, n_c])
    src, out = tmp_path / 'in.h5', tmp_path / 'out.h5'
    write_input(src, n_prof, (n_prof > 0).sum(axis=0))
    new_grid(str(src), str(out), S_OLD)
    with h5py.File(out, 'r') as f:
        assert f['rel_stages'][()].tolist() == [3, 3, 3, 3, 3, 2]
        assert f['n_prof'][2, -1] == 0.0


def test_extrapolation_compares_with_endpoint(tmp_path):
    # Beyond the last old surface only that surface counts (2 species);
    # the spline revives carbon at s = 0.945, which must be refused.
    s_old = 0.9 * S_OLD
    n_e = np.full(S_OLD.size, 4.0e13)
    n_c = 1.0e11 * np.array([1.0, 1.0, 1.0, 10.0, 1.0, 0.0])
    n_prof = np.vstack([n_e, n_e - 6.0 * n_c, n_c])
    src = tmp_path / 'in.h5'
    write_input(src, n_prof, (n_prof > 0).sum(axis=0))
    with h5py.File(src, 'a') as f:
        f['boozer_s'][...] = s_old
    with pytest.raises(ValueError, match=r'neighbouring.*0\.945'):
        new_grid(str(src), str(tmp_path / 'out.h5'),
                 np.append(s_old[:-1], 0.945))


def test_rel_stages_unchanged_when_constant(tmp_path):
    n_e = 4.0e13 * (1.2 - S_OLD**2)
    n_prof = np.vstack([n_e, 0.9 * n_e, 0.0125 * n_e])
    src, out = tmp_path / 'in.h5', tmp_path / 'out.h5'
    write_input(src, n_prof, [3] * S_OLD.size)
    new_grid(str(src), str(out), S_NEW)
    with h5py.File(src, 'r') as fin, h5py.File(out, 'r') as fout:
        assert np.array_equal(fout['rel_stages'][()], fin['rel_stages'][()])


@pytest.mark.parametrize('index, z_ion', [(2, 1.0), (1, 6.0)])
def test_remove_species_leaving_two_recomputes_ion_kappa(tmp_path, index,
                                                          z_ion):
    # Remove carbon (keep D, Z = 1) or deuterium (keep C, Z = 6); electron
    # and ion temperatures differ so the Coulomb logarithm source matters.
    n_e = 4.0e13 * (1.2 - S_OLD**2)
    n_prof = np.vstack([n_e, 0.9 * n_e, 0.0125 * n_e])
    src, out = tmp_path / 'in.h5', tmp_path / 'out.h5'
    write_input(src, n_prof, [3] * S_OLD.size, t_scale=[1.0, 0.7, 0.5])
    remove_species_from_profile_file(str(src), str(out), index)
    with h5py.File(src, 'r') as fin, h5py.File(out, 'r') as fout:
        assert fout['num_species'][()].tolist() == [2]
        assert fout['rel_stages'][()].tolist() == [2] * S_OLD.size
        kappa = fout['kappa_prof'][()]
        n_out, t_out = fout['n_prof'][()], fout['T_prof'][()]
        assert kappa.shape == (2, S_OLD.size)
        # Quasi-neutrality with the remaining ion: n_i = n_e / Z.
        assert np.allclose(n_out[1], n_e / z_ion, rtol=1e-15, atol=0)
        assert np.array_equal(kappa[0], fin['kappa_prof'][()][0])
        # Independent evaluation: 2 / mean free path with the electron
        # Coulomb logarithm, as written by the input generator.
        log_lambda = 52.43 - 1.15 * np.log10(n_e) + 2.3 * np.log10(t_out[0])
        charge = z_ion * E_CGS
        mfp = (3.0 / (4.0 * np.sqrt(np.pi)) * (t_out[1] / charge)**2
               / (n_out[1] * charge**2 * log_lambda))
        assert np.allclose(kappa[1], 2.0 / mfp, rtol=1e-12, atol=0)


def test_remove_locally_absent_species_keeps_counts(tmp_path):
    # Carbon exists only inside s = 0.4; removing it must leave 2 species on
    # every surface, not 1 where carbon was already absent.
    n_e = 4.0e13 * (1.2 - S_OLD**2)
    n_c = 1.0e12 * np.array([2.0, 1.8, 1.0, 0.0, 0.0, 0.0])
    n_prof = np.vstack([n_e, n_e - 6.0 * n_c, n_c, 1.0e-3 * n_e])
    src, out = tmp_path / 'in.h5', tmp_path / 'out.h5'
    write_input(src, n_prof, (n_prof > 0).sum(axis=0))
    remove_species_from_profile_file(str(src), str(out), 2)
    with h5py.File(out, 'r') as f:
        assert f['rel_stages'][()].tolist() == [3] * S_OLD.size
        assert np.array_equal(f['n_prof'][()][2], 1.0e-3 * n_e)
