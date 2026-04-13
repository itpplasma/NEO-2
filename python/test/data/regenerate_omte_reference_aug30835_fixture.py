"""Regenerate the AUG 30835 Om_tE reference fixture from real NEO-2 outputs."""

from pathlib import Path
import re

import h5py
import numpy as np


FIXTURE_PATH = Path(__file__).with_name('omte_reference_aug30835.npz')
RUN_DIRS = [
    Path('/home/ert/data/AUG/NEO-2/30835/gorilla_axisymmetric_baseline/runs/ql_two_surfaces/es_0p25271'),
    Path('/home/ert/data/AUG/NEO-2/30835/gorilla_axisymmetric_baseline/runs/ql_two_surfaces/es_0p49841'),
]
SOURCE_H5_PATHS = [run_dir / 'neo2_multispecies_out.h5' for run_dir in RUN_DIRS]
SOURCE_NEO2_INPUT_PATHS = [run_dir / 'neo2.in' for run_dir in RUN_DIRS]

SURFACE_FIELDS = [
    'Bref',
    'Dp0',
    'D31ref0',
    'Er',
    'aiota',
    'avEparB_ov_avb2',
    'av_gphph',
    'av_inv_bhat',
    'av_nabla_stor',
    'avbhat',
    'avbhat2',
    'bcovar_phi',
    'bcovar_tht',
    'boozer_s',
    'dbcovar_phi_ds',
    'dbcovar_theta_ds',
    'diota_ds',
    'eps_M_2',
    'm_phi',
    'psi_pr_hat',
    'sqrtg_bctrvr_phi',
    'sqrtg_bctrvr_tht',
]

INVARIANT_VECTOR_FIELDS = [
    'species_tag',
    'z_spec',
    'm_spec',
]

STACKED_VECTOR_FIELDS = [
    'n_spec',
    'T_spec',
    'collpar_spec',
    'nu_star_spec',
    'MtOvR',
    'VphiB_spec',
    'VphiB_Ware_spec',
    'VthtB_spec',
    'VthtB_Ware_spec',
    'Gamma_AX_spec',
    'Gamma_AX_Ware_spec',
    'Gamma_NA_spec',
    'Gamma_NA_Ware_spec',
    'Qflux_AX_spec',
    'Qflux_AX_Ware_spec',
    'Qflux_NA_spec',
    'Qflux_NA_Ware_spec',
    'ParFlow_AX_spec',
    'ParFlow_AX_Ware_spec',
    'ParFlow_NA_spec',
    'ParFlow_NA_Ware_spec',
    'TphiNA_spec',
    'TphiNA_Ware_spec',
]

STACKED_MATRIX_FIELDS = [
    'row_ind_spec',
    'col_ind_spec',
    'D11_AX',
    'D11_AX_Dpl',
    'D11_NA',
    'D11_NA_Dpl',
    'D12_AX',
    'D12_AX_Dpl',
    'D12_NA',
    'D12_NA_Dpl',
    'D13_AX',
    'D13_AX_D31ref',
    'D13_NA',
    'D13_NA_D31ref',
    'D21_AX',
    'D21_AX_Dpl',
    'D21_NA',
    'D21_NA_Dpl',
    'D22_AX',
    'D22_AX_Dpl',
    'D22_NA',
    'D22_NA_Dpl',
    'D23_AX',
    'D23_AX_D31ref',
    'D23_NA',
    'D23_NA_D31ref',
    'D31_AX',
    'D31_AX_D31ref',
    'D31_NA',
    'D31_NA_D31ref',
    'D32_AX',
    'D32_AX_D31ref',
    'D32_NA',
    'D32_NA_D31ref',
    'D33_AX',
    'D33_AX_norm',
    'D33_NA',
    'D33_NA_norm',
]

PROFILE_VECTOR_FIELDS = {
    'n_prof': 'n_vec',
    'T_prof': 't_vec',
    'dn_ov_ds_prof': 'dn_vec_ov_ds',
    'dT_ov_ds_prof': 'dt_vec_ov_ds',
}


def _extract_namelist_value(text, name):
    match = re.search(rf'^\s*{name}\s*=\s*(.+?)\s*!', text, re.M)
    if not match:
        raise ValueError(f'missing {name} in neo2.in')
    return match.group(1).strip()


def _parse_repeated_list(expr, dtype=float):
    values = []
    for token in [part.strip() for part in expr.split(',')]:
        if not token:
            continue
        if '*' in token:
            count_text, value_text = token.split('*', 1)
            values.extend([dtype(value_text)] * int(count_text))
        else:
            values.append(dtype(token))
    return np.asarray(values)


def _load_multispec_input(path):
    text = path.read_text()
    parsed = {
        target_name: _parse_repeated_list(_extract_namelist_value(text, source_name))
        for target_name, source_name in PROFILE_VECTOR_FIELDS.items()
    }
    parsed['species_tag'] = _parse_repeated_list(
        _extract_namelist_value(text, 'species_tag_vec'),
        dtype=int,
    )
    parsed['species_tag_Vphi'] = int(float(_extract_namelist_value(text, 'species_tag_vphi')))
    parsed['z_spec'] = _parse_repeated_list(_extract_namelist_value(text, 'z_vec'))
    parsed['m_spec'] = _parse_repeated_list(_extract_namelist_value(text, 'm_vec'))
    parsed['Vphi'] = float(_extract_namelist_value(text, 'vphi'))
    return parsed


def _species_def(z_spec, m_spec, num_surfaces):
    species_def = np.empty((2, z_spec.size, num_surfaces), dtype=float)
    species_def[0] = np.repeat(z_spec[:, np.newaxis], num_surfaces, axis=1)
    species_def[1] = np.repeat(m_spec[:, np.newaxis], num_surfaces, axis=1)
    return species_def


def _scalar_field(path, field_name):
    with h5py.File(path, 'r') as handle:
        return float(np.asarray(handle[field_name]).reshape(-1)[0])


def _stack_field(path_list, field_name):
    values = []
    for path in path_list:
        with h5py.File(path, 'r') as handle:
            values.append(np.asarray(handle[field_name]))
    return np.stack(values, axis=0)


def main():
    existing = dict(np.load(FIXTURE_PATH))
    regenerated = dict(existing)

    for path in SOURCE_H5_PATHS:
        if not path.exists():
            raise FileNotFoundError(path)
    for path in SOURCE_NEO2_INPUT_PATHS:
        if not path.exists():
            raise FileNotFoundError(path)

    parsed_inputs = [_load_multispec_input(path) for path in SOURCE_NEO2_INPUT_PATHS]
    num_surfaces = len(parsed_inputs)
    invariant = parsed_inputs[0]

    regenerated['source_h5_paths'] = np.array([str(path) for path in SOURCE_H5_PATHS])
    regenerated['source_neo2in_paths'] = np.array([str(path) for path in SOURCE_NEO2_INPUT_PATHS])

    regenerated['species_tag'] = invariant['species_tag']
    regenerated['species_tag_Vphi'] = invariant['species_tag_Vphi']
    regenerated['z_spec'] = invariant['z_spec']
    regenerated['m_spec'] = invariant['m_spec']
    regenerated['species_def'] = _species_def(
        z_spec=invariant['z_spec'],
        m_spec=invariant['m_spec'],
        num_surfaces=num_surfaces,
    )

    for field_name in PROFILE_VECTOR_FIELDS:
        regenerated[field_name] = np.stack(
            [parsed[field_name] for parsed in parsed_inputs],
            axis=0,
        )
    regenerated['Vphi'] = np.array([parsed['Vphi'] for parsed in parsed_inputs], dtype=float)

    for field_name in SURFACE_FIELDS:
        regenerated[field_name] = np.array([
            _scalar_field(path, field_name) for path in SOURCE_H5_PATHS
        ])

    regenerated['Er_neo2'] = regenerated['Er'].copy()
    regenerated['avb2'] = regenerated['avbhat2'] * regenerated['Bref'] ** 2

    for field_name in STACKED_VECTOR_FIELDS + STACKED_MATRIX_FIELDS:
        regenerated[field_name] = _stack_field(SOURCE_H5_PATHS, field_name)

    if not np.allclose(regenerated['n_prof'], regenerated['n_spec'], rtol=0.0, atol=0.0):
        raise ValueError('n_prof parsed from neo2.in does not match n_spec from HDF5 output')
    if not np.allclose(regenerated['T_prof'], regenerated['T_spec'], rtol=0.0, atol=0.0):
        raise ValueError('T_prof parsed from neo2.in does not match T_spec from HDF5 output')

    np.savez(FIXTURE_PATH, **regenerated)
    print(FIXTURE_PATH)


if __name__ == '__main__':
    main()