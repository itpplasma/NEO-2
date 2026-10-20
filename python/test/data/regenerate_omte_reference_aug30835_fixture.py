"""Regenerate the AUG 30835 Om_tE reference fixture from real NEO-2 outputs."""

import os
from pathlib import Path
import re

import h5py
import numpy as np


FIXTURE_PATH = Path(__file__).with_name('omte_reference_aug30835.npz')
DATA_ROOT = Path(os.environ.get('NEO2_DATA_ROOT', '/home/ert/data'))
RUN_BASE = DATA_ROOT / 'AUG/NEO-2/30835/gorilla_axisymmetric_baseline/runs/ql_two_surfaces'
RUN_DIRS = [
    RUN_BASE / 'es_0p25271',
    RUN_BASE / 'es_0p49841',
]
SOURCE_H5_PATHS = [run_dir / 'neo2_multispecies_out.h5' for run_dir in RUN_DIRS]
SOURCE_NEO2_INPUT_PATHS = [run_dir / 'neo2.in' for run_dir in RUN_DIRS]

SURFACE_FIELDS = [
    'Er', 'aiota', 'avEparB_ov_avb2', 'av_nabla_stor',
    'bcovar_phi', 'bcovar_tht', 'boozer_s', 'sqrtg_bctrvr_tht',
]

STACKED_FIELDS = [
    'MtOvR', 'row_ind_spec', 'col_ind_spec', 'D31_AX', 'D32_AX', 'D33_AX',
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
    regenerated = {}

    for path in SOURCE_H5_PATHS:
        if not path.exists():
            raise FileNotFoundError(path)
    for path in SOURCE_NEO2_INPUT_PATHS:
        if not path.exists():
            raise FileNotFoundError(path)

    parsed_inputs = [_load_multispec_input(path) for path in SOURCE_NEO2_INPUT_PATHS]
    invariant = parsed_inputs[0]

    regenerated['source_h5_paths'] = np.array([str(path) for path in SOURCE_H5_PATHS])
    regenerated['source_neo2in_paths'] = np.array([str(path) for path in SOURCE_NEO2_INPUT_PATHS])

    regenerated['species_tag'] = invariant['species_tag']
    regenerated['species_tag_Vphi'] = invariant['species_tag_Vphi']
    regenerated['z_spec'] = invariant['z_spec']
    regenerated['m_spec'] = invariant['m_spec']
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

    for field_name in STACKED_FIELDS:
        regenerated[field_name] = _stack_field(SOURCE_H5_PATHS, field_name)

    if not np.allclose(regenerated['n_prof'], _stack_field(SOURCE_H5_PATHS, 'n_spec'), rtol=0.0, atol=0.0):
        raise ValueError('n_prof parsed from neo2.in does not match n_spec from HDF5 output')
    if not np.allclose(regenerated['T_prof'], _stack_field(SOURCE_H5_PATHS, 'T_spec'), rtol=0.0, atol=0.0):
        raise ValueError('T_prof parsed from neo2.in does not match T_spec from HDF5 output')

    np.savez(FIXTURE_PATH, **regenerated)
    print(FIXTURE_PATH)


if __name__ == '__main__':
    main()
