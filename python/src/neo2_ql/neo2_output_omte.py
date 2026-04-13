"""Reconstruct Om_tE directly from NEO-2 output data."""

from pathlib import Path

import h5py
import numpy as np

from .compute_omte import C_CGS, E_CGS


def _species_index_from_tag(species_tag, species_tag_vphi):
    species_tag = np.asarray(species_tag, dtype=int)
    matches = np.where(species_tag == int(species_tag_vphi))[0]
    if matches.size != 1:
        raise ValueError('species_tag_Vphi must match exactly one species')
    return int(matches[0])


def _remap_output_indices(indices, species_tag):
    indices = np.asarray(indices, dtype=int)
    species_tag = np.asarray(species_tag, dtype=int)
    num_species = species_tag.size

    if np.all((indices >= 0) & (indices < num_species)):
        return indices

    tag_to_index = {int(tag): idx for idx, tag in enumerate(species_tag)}
    try:
        return np.array([tag_to_index[int(value)] for value in indices], dtype=int)
    except KeyError as exc:
        raise ValueError('row/col indices do not map to species_tag') from exc


def compute_neo2_er_from_transport_coefficients(
    n_spec,
    T_spec,
    dn_spec_ov_ds,
    dT_spec_ov_ds,
    species_tag,
    species_tag_vphi,
    z_spec,
    Vphi,
    aiota,
    sqrtg_bctrvr_phi,
    av_nabla_stor,
    bcovar_tht,
    bcovar_phi,
    row_ind,
    col_ind,
    D31_AX,
    D32_AX,
    D33_AX=None,
    avEparB_ov_avb2=0.0,
    isw_Vphi_loc=0,
):
    """Mirror NEO-2 `compute_Er()` for the tested `isw_Vphi_loc=0` branch."""
    if int(isw_Vphi_loc) != 0:
        raise NotImplementedError('Only isw_Vphi_loc=0 is implemented')

    n_spec = np.asarray(n_spec, dtype=float)
    T_spec = np.asarray(T_spec, dtype=float)
    dn_spec_ov_ds = np.asarray(dn_spec_ov_ds, dtype=float)
    dT_spec_ov_ds = np.asarray(dT_spec_ov_ds, dtype=float)
    species_tag = np.asarray(species_tag, dtype=int)
    z_spec = np.asarray(z_spec, dtype=float)
    row_ind = _remap_output_indices(row_ind, species_tag)
    col_ind = _remap_output_indices(col_ind, species_tag)
    D31_AX = np.asarray(D31_AX, dtype=float)
    D32_AX = np.asarray(D32_AX, dtype=float)
    D33_AX = np.zeros_like(D31_AX) if D33_AX is None else np.asarray(D33_AX, dtype=float)

    if np.any(n_spec == 0.0):
        raise ValueError('n_spec must be nonzero to reconstruct NEO-2 E_r')
    if np.any(T_spec == 0.0):
        raise ValueError('T_spec must be nonzero to reconstruct NEO-2 E_r')
    if aiota == 0.0:
        raise ValueError('aiota must be nonzero to reconstruct NEO-2 E_r')
    if sqrtg_bctrvr_phi == 0.0:
        raise ValueError('sqrtg_bctrvr_phi must be nonzero to reconstruct NEO-2 E_r')
    if bcovar_tht == 0.0:
        raise ValueError('bcovar_tht must be nonzero to reconstruct NEO-2 E_r')

    spec_i = _species_index_from_tag(species_tag, species_tag_vphi)

    z_ions = z_spec[spec_i]
    T_ions = T_spec[spec_i]
    n_ions = n_spec[spec_i]
    dT_ions_ov_dr = dT_spec_ov_ds[spec_i] * av_nabla_stor
    dn_ions_ov_dr = dn_spec_ov_ds[spec_i] * av_nabla_stor
    p_ions = n_ions * T_ions
    if p_ions == 0.0:
        raise ValueError('ion pressure must be nonzero to reconstruct NEO-2 E_r')
    dp_ions_ov_dr = T_ions * dn_ions_ov_dr + n_ions * dT_ions_ov_dr

    denom_er = C_CGS * bcovar_tht / sqrtg_bctrvr_phi
    nom_er = (
        Vphi * (aiota * bcovar_tht + bcovar_phi)
        + (C_CGS * T_ions * bcovar_tht / (z_ions * E_CGS * sqrtg_bctrvr_phi))
        * (dp_ions_ov_dr / p_ions)
    )

    for idx in range(D31_AX.size):
        irow_spec = row_ind[idx]
        icol_spec = col_ind[idx]
        if irow_spec == spec_i:
            denom_er += D31_AX[idx] * (z_spec[icol_spec] * E_CGS) / T_spec[icol_spec]
            nom_er += av_nabla_stor * D31_AX[idx] * (
                dn_spec_ov_ds[icol_spec] / n_spec[icol_spec]
                + dT_spec_ov_ds[icol_spec] / T_spec[icol_spec]
            )
            nom_er += av_nabla_stor * (
                dT_spec_ov_ds[icol_spec] / T_spec[icol_spec]
            ) * (D32_AX[idx] - 2.5 * D31_AX[idx])
            nom_er += D33_AX[idx] * avEparB_ov_avb2 * (
                z_spec[icol_spec] * E_CGS
            ) / T_spec[icol_spec]

    return nom_er / denom_er


def compute_neo2_omte_from_transport_coefficients(**kwargs):
    """Reconstruct Om_tE by mirroring NEO-2 `compute_Er()`."""
    er = compute_neo2_er_from_transport_coefficients(**kwargs)
    om_tE = C_CGS * er / (kwargs['aiota'] * kwargs['sqrtg_bctrvr_phi'])
    return om_tE, er


def _dataset_or_none(handle, name):
    return None if name not in handle else np.asarray(handle[name])


def _surface_vector(data, surface_index):
    data = np.asarray(data)
    if data.ndim <= 1:
        return data
    return data[surface_index]


def _surface_z_spec(species_def, species_tag, surface_index, num_surfaces):
    species_def = np.asarray(species_def)
    num_species = np.asarray(species_tag).size

    if species_def.ndim == 3:
        z_data = species_def[0]
    elif species_def.ndim == 2:
        z_data = species_def
    else:
        raise ValueError('species_def must be 2D or 3D')

    if z_data.ndim == 1:
        return z_data
    if z_data.shape == (num_species, num_surfaces):
        return z_data[:, surface_index]
    if z_data.shape == (num_surfaces, num_species):
        return z_data[surface_index]
    raise ValueError('species_def shape is incompatible with species/surface axes')


def _species_state_from_handle(handle, surface_index):
    n_spec = _dataset_or_none(handle, 'n_spec')
    T_spec = _dataset_or_none(handle, 'T_spec')
    dn_spec_ov_ds = _dataset_or_none(handle, 'dn_spec_ov_ds')
    dT_spec_ov_ds = _dataset_or_none(handle, 'dT_spec_ov_ds')

    if n_spec is None:
        n_prof = np.asarray(handle['n_prof'])
        n_spec = _surface_vector(n_prof, surface_index)
    else:
        n_spec = _surface_vector(n_spec, surface_index)

    if T_spec is None:
        T_prof = np.asarray(handle['T_prof'])
        T_spec = _surface_vector(T_prof, surface_index)
    else:
        T_spec = _surface_vector(T_spec, surface_index)

    if dn_spec_ov_ds is None:
        dn_prof = np.asarray(handle['dn_ov_ds_prof'])
        dn_spec_ov_ds = _surface_vector(dn_prof, surface_index)
    else:
        dn_spec_ov_ds = _surface_vector(dn_spec_ov_ds, surface_index)

    if dT_spec_ov_ds is None:
        dT_prof = np.asarray(handle['dT_ov_ds_prof'])
        dT_spec_ov_ds = _surface_vector(dT_prof, surface_index)
    else:
        dT_spec_ov_ds = _surface_vector(dT_spec_ov_ds, surface_index)

    return n_spec, T_spec, dn_spec_ov_ds, dT_spec_ov_ds


def compute_omte_from_neo2_output(path, mode='stored'):
    """Return Om_tE from a NEO-2 output file.

    Parameters
    ----------
    path : str or path-like
        Path to a `neo2_multispecies_out.h5` or combined output file.
    mode : {'stored', 'transport'}
        `stored` reads the `Er` dataset directly.
        `transport` reconstructs `Er` from the stored transport coefficients
        and the stored `avEparB_ov_avb2` term for `isw_Vphi_loc=0`.
    """
    path = Path(path)
    with h5py.File(path, 'r') as handle:
        aiota = np.asarray(handle['aiota'])
        sqrtg_bctrvr_phi = np.asarray(handle['sqrtg_bctrvr_phi'])

        if mode == 'stored':
            er = np.asarray(handle['Er'])
            om_tE = C_CGS * er / (aiota * sqrtg_bctrvr_phi)
            return om_tE, er

        if mode != 'transport':
            raise ValueError("mode must be 'stored' or 'transport'")

        species_tag = np.asarray(handle['species_tag'])
        species_tag_vphi = np.asarray(handle['species_tag_Vphi']).reshape(-1)[0]
        species_def = np.asarray(handle['species_def'])

        row_ind = np.asarray(handle['row_ind_spec'])
        col_ind = np.asarray(handle['col_ind_spec'])
        D31_AX = np.asarray(handle['D31_AX'])
        D32_AX = np.asarray(handle['D32_AX'])
        D33_AX = np.asarray(handle['D33_AX'])
        avEparB = np.asarray(handle['avEparB_ov_avb2'])
        Vphi = np.asarray(handle['Vphi'])
        av_nabla_stor = np.asarray(handle['av_nabla_stor'])
        bcovar_tht = np.asarray(handle['bcovar_tht'])
        bcovar_phi = np.asarray(handle['bcovar_phi'])
        isw_Vphi_loc = int(np.asarray(handle['isw_Vphi_loc']).reshape(-1)[0])

        num_surfaces = int(np.atleast_1d(aiota).shape[0])
        er = np.empty(num_surfaces, dtype=float)
        for surface_index in range(num_surfaces):
            n_spec, T_spec, dn_spec_ov_ds, dT_spec_ov_ds = _species_state_from_handle(
                handle, surface_index
            )
            z_spec = _surface_z_spec(species_def, species_tag, surface_index, num_surfaces)
            er[surface_index] = compute_neo2_er_from_transport_coefficients(
                n_spec=n_spec,
                T_spec=T_spec,
                dn_spec_ov_ds=dn_spec_ov_ds,
                dT_spec_ov_ds=dT_spec_ov_ds,
                species_tag=species_tag,
                species_tag_vphi=species_tag_vphi,
                z_spec=z_spec,
                Vphi=np.atleast_1d(Vphi)[surface_index],
                aiota=np.atleast_1d(aiota)[surface_index],
                sqrtg_bctrvr_phi=np.atleast_1d(sqrtg_bctrvr_phi)[surface_index],
                av_nabla_stor=np.atleast_1d(av_nabla_stor)[surface_index],
                bcovar_tht=np.atleast_1d(bcovar_tht)[surface_index],
                bcovar_phi=np.atleast_1d(bcovar_phi)[surface_index],
                row_ind=_surface_vector(row_ind, surface_index),
                col_ind=_surface_vector(col_ind, surface_index),
                D31_AX=_surface_vector(D31_AX, surface_index),
                D32_AX=_surface_vector(D32_AX, surface_index),
                D33_AX=_surface_vector(D33_AX, surface_index),
                avEparB_ov_avb2=np.atleast_1d(avEparB)[surface_index],
                isw_Vphi_loc=isw_Vphi_loc,
            )

        om_tE = C_CGS * er / np.atleast_1d(aiota) / np.atleast_1d(sqrtg_bctrvr_phi)
        return om_tE, er
