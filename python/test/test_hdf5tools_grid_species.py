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

from neo2_ql import get_coulomb_logarithm, get_kappa
from neo2_util.hdf5tools import new_grid, remove_species_from_profile_file

E_CGS = 1.60217662e-19 * 2.99792458e8 * 10
S_OLD = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
# Shifted surfaces, same count (new_grid writes datasets in place).
S_NEW = np.array([0.05, 0.3, 0.5, 0.55, 0.7, 0.95])


def write_input(path, n_prof, rel_stages):
    nspec, nrad = n_prof.shape
    t_prof = np.tile(1.0e-9 * (2.0 - S_OLD), (nspec, 1))
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


def test_rel_stages_follow_species_with_density(tmp_path):
    # Carbon (species 3) is present inside s = 0.4 and absent outside.
    n_e = 4.0e13 * (1.2 - S_OLD**2)
    n_c = 1.0e12 * np.array([2.0, 1.8, 1.0, 0.0, 0.0, 0.0])
    n_prof = np.vstack([n_e, n_e - 6.0 * n_c, n_c])
    rel_old = (n_prof > 0).sum(axis=0)
    assert rel_old.tolist() == [3, 3, 3, 2, 2, 2]
    src, out = tmp_path / 'in.h5', tmp_path / 'out.h5'
    write_input(src, n_prof, rel_old)
    new_grid(str(src), str(out), S_NEW)
    with h5py.File(out, 'r') as f:
        rel_new = f['rel_stages'][()]
        n_new = f['n_prof'][()]
        assert f['rel_stages'].dtype == np.int32
    # Contract of the Fortran reader on every surface.
    assert rel_new.tolist() == (n_new > 0).sum(axis=0).tolist()
    # Never fewer species than both neighbouring old surfaces.
    for k, s in enumerate(S_NEW):
        right = np.searchsorted(S_OLD, s)
        left = max(right - 1, 0)
        assert rel_new[k] >= min(rel_old[left], rel_old[right]), (s, rel_new)
    # Inside the carbon region all three species stay.
    assert rel_new[:2].tolist() == [3, 3]


def test_rel_stages_unchanged_when_constant(tmp_path):
    n_e = 4.0e13 * (1.2 - S_OLD**2)
    n_prof = np.vstack([n_e, 0.9 * n_e, 0.0125 * n_e])
    src, out = tmp_path / 'in.h5', tmp_path / 'out.h5'
    write_input(src, n_prof, [3] * S_OLD.size)
    new_grid(str(src), str(out), S_NEW)
    with h5py.File(src, 'r') as fin, h5py.File(out, 'r') as fout:
        assert np.array_equal(fout['rel_stages'][()], fin['rel_stages'][()])


def test_remove_species_leaving_two_recomputes_ion_kappa(tmp_path):
    n_e = 4.0e13 * (1.2 - S_OLD**2)
    n_prof = np.vstack([n_e, 0.9 * n_e, 0.0125 * n_e])
    src, out = tmp_path / 'in.h5', tmp_path / 'out.h5'
    write_input(src, n_prof, [3] * S_OLD.size)
    remove_species_from_profile_file(str(src), str(out), 2)
    with h5py.File(src, 'r') as fin, h5py.File(out, 'r') as fout:
        assert fout['num_species'][()].tolist() == [2]
        kappa = fout['kappa_prof'][()]
        n_out, t_out = fout['n_prof'][()], fout['T_prof'][()]
        assert kappa.shape == (2, S_OLD.size)
        # Two species left: the ion takes the electron density (Z = 1).
        assert np.array_equal(n_out[1], n_e)
        assert np.array_equal(kappa[0], fin['kappa_prof'][()][0])
        # Independent evaluation: 2 / mean free path with the electron
        # Coulomb logarithm, as written by the input generator.
        log_lambda = 52.43 - 1.15 * np.log10(n_e) + 2.3 * np.log10(t_out[0])
        mfp = (3.0 / (4.0 * np.sqrt(np.pi)) * (t_out[1] / E_CGS)**2
               / (n_out[1] * E_CGS**2 * log_lambda))
        assert np.allclose(kappa[1], 2.0 / mfp, rtol=1e-12, atol=0)
