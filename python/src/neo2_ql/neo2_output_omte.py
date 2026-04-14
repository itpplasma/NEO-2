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
    return decompose_neo2_er_transport_terms(
        n_spec=n_spec,
        T_spec=T_spec,
        dn_spec_ov_ds=dn_spec_ov_ds,
        dT_spec_ov_ds=dT_spec_ov_ds,
        species_tag=species_tag,
        species_tag_vphi=species_tag_vphi,
        z_spec=z_spec,
        Vphi=Vphi,
        aiota=aiota,
        sqrtg_bctrvr_phi=sqrtg_bctrvr_phi,
        av_nabla_stor=av_nabla_stor,
        bcovar_tht=bcovar_tht,
        bcovar_phi=bcovar_phi,
        row_ind=row_ind,
        col_ind=col_ind,
        D31_AX=D31_AX,
        D32_AX=D32_AX,
        D33_AX=D33_AX,
        avEparB_ov_avb2=avEparB_ov_avb2,
        isw_Vphi_loc=isw_Vphi_loc,
    )['er_total']


def decompose_neo2_er_transport_terms(
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
    """Return the exact term breakdown of the tested NEO-2 `compute_Er()` branch."""
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

    denom_base = C_CGS * bcovar_tht / sqrtg_bctrvr_phi
    nom_dia = (
        C_CGS * T_ions * bcovar_tht / (z_ions * E_CGS * sqrtg_bctrvr_phi)
    ) * (dp_ions_ov_dr / p_ions)
    nom_vphi = Vphi * (aiota * bcovar_tht + bcovar_phi)
    nom_d31 = 0.0
    nom_d32 = 0.0
    nom_d33 = 0.0
    denom_d31 = 0.0

    for idx in range(D31_AX.size):
        irow_spec = row_ind[idx]
        icol_spec = col_ind[idx]
        if irow_spec == spec_i:
            denom_d31 += D31_AX[idx] * (z_spec[icol_spec] * E_CGS) / T_spec[icol_spec]
            nom_d31 += av_nabla_stor * D31_AX[idx] * (
                dn_spec_ov_ds[icol_spec] / n_spec[icol_spec]
                + dT_spec_ov_ds[icol_spec] / T_spec[icol_spec]
            )
            nom_d32 += av_nabla_stor * (
                dT_spec_ov_ds[icol_spec] / T_spec[icol_spec]
            ) * (D32_AX[idx] - 2.5 * D31_AX[idx])
            nom_d33 += D33_AX[idx] * avEparB_ov_avb2 * (
                z_spec[icol_spec] * E_CGS
            ) / T_spec[icol_spec]

    nom_base = nom_dia + nom_vphi
    nom_total = nom_base + nom_d31 + nom_d32 + nom_d33
    denom_total = denom_base + denom_d31

    er_dia = nom_dia / denom_base
    er_vphi = nom_vphi / denom_base
    er_d31 = nom_d31 / denom_base
    er_d32 = nom_d32 / denom_base
    er_d33 = nom_d33 / denom_base
    er_before_denom = nom_total / denom_base
    er_denom = (nom_total / denom_total) - er_before_denom
    er_total = nom_total / denom_total

    return {
        'species_index': spec_i,
        'dp_ions_ov_dr': dp_ions_ov_dr,
        'p_ions': p_ions,
        'denom_base': denom_base,
        'denom_d31': denom_d31,
        'denom_total': denom_total,
        'nom_dia': nom_dia,
        'nom_vphi': nom_vphi,
        'nom_d31': nom_d31,
        'nom_d32': nom_d32,
        'nom_d33': nom_d33,
        'nom_total': nom_total,
        'er_dia': er_dia,
        'er_vphi': er_vphi,
        'er_d31': er_d31,
        'er_d32': er_d32,
        'er_d33': er_d33,
        'er_denom': er_denom,
        'er_before_denom': er_before_denom,
        'er_total': er_total,
    }


def compute_neo2_omte_from_transport_coefficients(**kwargs):
    """Reconstruct Om_tE by mirroring NEO-2 `compute_Er()`."""
    er = compute_neo2_er_from_transport_coefficients(**kwargs)
    om_tE = C_CGS * er / (kwargs['aiota'] * kwargs['sqrtg_bctrvr_phi'])
    return om_tE, er


def _electron_species_index(z_spec):
    z_spec = np.asarray(z_spec, dtype=float)
    matches = np.where(z_spec < 0.0)[0]
    if matches.size == 0:
        raise ValueError('z_spec must contain at least one electron species')
    return int(matches[0])


def _normalize_z_spec(z_spec, num_species):
    z_spec = np.asarray(z_spec, dtype=float)
    if z_spec.ndim == 1:
        return z_spec
    if z_spec.ndim == 2 and z_spec.shape[1] == num_species:
        if np.allclose(z_spec, z_spec[0], rtol=0.0, atol=0.0):
            return z_spec[0]
    raise ValueError('z_spec must be 1D or surface-invariant along the first axis')


def compute_d31_reference_electron(T_e, z_e, aiota, sqrtg_bctrvr_phi, bcovar_phi):
    """Return the D31 reference coefficient used for normalized transport entries.

    This mirrors the effective reference scale written by the Fortran transport
    pipeline: the stored normalized coefficients `D31_AX_D31ref` and
    `D32_AX_D31ref` can be converted back to raw `D31_AX` and `D32_AX` by
    multiplying with this factor.
    """
    if z_e == 0.0:
        raise ValueError('z_e must be nonzero to compute D31 reference')
    if aiota == 0.0:
        raise ValueError('aiota must be nonzero to compute D31 reference')
    if sqrtg_bctrvr_phi == 0.0:
        raise ValueError('sqrtg_bctrvr_phi must be nonzero to compute D31 reference')
    return C_CGS * T_e * bcovar_phi / (z_e * E_CGS * aiota * sqrtg_bctrvr_phi)


def decompose_neo2_er_k_cof_transport_model(
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
    d31_hat,
    k_cof,
    transport_species_tag=None,
    d33_ax=0.0,
    avEparB_ov_avb2=0.0,
):
    """Evaluate a minimal Level 2.5 transport model on the exact `compute_Er()` algebra.

    This is the smallest transport closure that still reuses the same numerator
    and denominator decomposition as the full NEO-2 replay. It only models the
    measured-ion row of the transport matrix:

    - `D31_AX = d31_hat * D31_ref`
    - `D32_AX = (2.5 - k_cof) * D31_AX`

    The default column species is the measured-ion species itself, which gives a
    compact ion-ion transport model. `k_cof` is taken as a direct user input so
    it can be swept across collisionality regimes without computing it first.
    """
    species_tag = np.asarray(species_tag, dtype=int)
    z_spec = _normalize_z_spec(z_spec, species_tag.size)
    spec_i = _species_index_from_tag(species_tag, species_tag_vphi)
    col_spec = spec_i
    if transport_species_tag is not None:
        col_spec = _species_index_from_tag(species_tag, transport_species_tag)

    electron_index = _electron_species_index(z_spec)
    d31_ref = compute_d31_reference_electron(
        T_e=float(np.asarray(T_spec, dtype=float)[electron_index]),
        z_e=float(z_spec[electron_index]),
        aiota=float(aiota),
        sqrtg_bctrvr_phi=float(sqrtg_bctrvr_phi),
        bcovar_phi=float(bcovar_phi),
    )
    d31_ax = float(d31_hat) * d31_ref
    d32_ax = (2.5 - float(k_cof)) * d31_ax

    return decompose_neo2_er_transport_terms(
        n_spec=n_spec,
        T_spec=T_spec,
        dn_spec_ov_ds=dn_spec_ov_ds,
        dT_spec_ov_ds=dT_spec_ov_ds,
        species_tag=species_tag,
        species_tag_vphi=species_tag_vphi,
        z_spec=z_spec,
        Vphi=Vphi,
        aiota=aiota,
        sqrtg_bctrvr_phi=sqrtg_bctrvr_phi,
        av_nabla_stor=av_nabla_stor,
        bcovar_tht=bcovar_tht,
        bcovar_phi=bcovar_phi,
        row_ind=np.array([spec_i], dtype=int),
        col_ind=np.array([col_spec], dtype=int),
        D31_AX=np.array([d31_ax], dtype=float),
        D32_AX=np.array([d32_ax], dtype=float),
        D33_AX=np.array([float(d33_ax)], dtype=float),
        avEparB_ov_avb2=avEparB_ov_avb2,
        isw_Vphi_loc=0,
    )


def compute_neo2_er_from_k_cof_transport_model(**kwargs):
    """Return `E_r` from the minimal `D31_hat` + `k_cof` transport model."""
    return decompose_neo2_er_k_cof_transport_model(**kwargs)['er_total']


def compute_neo2_omte_from_k_cof_transport_model(**kwargs):
    """Return `Omega_tE` from the minimal `D31_hat` + `k_cof` transport model."""
    er = compute_neo2_er_from_k_cof_transport_model(**kwargs)
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


def _get_output_handle(handle):
    if 'Er' in handle:
        return handle
    required_transport_fields = {'aiota', 'sqrtg_bctrvr_phi', 'row_ind_spec', 'D31_AX'}
    if required_transport_fields.issubset(handle.keys()):
        return handle
    if 'neo2_multispecies_out' in handle:
        return handle['neo2_multispecies_out']
    raise KeyError('Could not find neo2_multispecies_out datasets in file')


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
        handle = _get_output_handle(handle)
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
        if 'z_spec' in handle:
            z_spec = np.asarray(handle['z_spec'])
        else:
            species_def = np.asarray(handle['species_def'])
            num_surfaces = int(np.atleast_1d(aiota).shape[0])
            z_spec = _surface_z_spec(species_def, species_tag, 0, num_surfaces)

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
            if 'z_spec' not in handle:
                z_spec = _surface_z_spec(
                    np.asarray(handle['species_def']),
                    species_tag,
                    surface_index,
                    num_surfaces,
                )
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
