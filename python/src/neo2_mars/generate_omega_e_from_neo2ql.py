from pathlib import Path

import h5py
import numpy as np


def _surface_species_array(data, num_surfaces):
    data = np.asarray(data)
    if data.ndim == 0:
        return data.reshape(1, 1)
    if data.ndim == 1:
        if num_surfaces == 1:
            return data.reshape(1, -1)
        if data.size == num_surfaces:
            return data.reshape(num_surfaces, 1)
        raise ValueError('Could not align 1D NEO-2 data with boozer_s surfaces')
    if data.ndim == 2:
        if data.shape[0] == num_surfaces:
            return data
        if data.shape[1] == num_surfaces:
            return data.T
    raise ValueError('Expected NEO-2 data with surface/species axes')


def generate_omega_e_for_mars(neo2ql_output_file, neo2ql_input_file, output_file='PROFWE.IN'):
    species_omega_e, stor = get_omega_e_from_neo2ql(neo2ql_output_file)
    sqrtspol = neo2ql_stor2sqrtspol(neo2ql_input_file, stor)
    omega_e = species_omega_e[:, 0]
    write_omega_e_to_mars_input(-omega_e, sqrtspol, output_file=output_file)


def get_omega_e_from_neo2ql(neo2ql_output_file):
    with h5py.File(neo2ql_output_file, 'r') as neo2ql:
        stor = np.atleast_1d(np.asarray(neo2ql['boozer_s'], dtype=float)).reshape(-1)
        num_surfaces = stor.size

        species_mach_over_major_radius = _surface_species_array(
            np.asarray(neo2ql['MtOvR'], dtype=float),
            num_surfaces,
        )
        species_temperature = _surface_species_array(
            np.asarray(neo2ql['T_spec'], dtype=float),
            num_surfaces,
        )
        species_mass = np.asarray(neo2ql['m_spec'], dtype=float).reshape(-1)
        species_omega_e = species_mach_over_major_radius * np.sqrt(
            2.0 * species_temperature / species_mass[np.newaxis, :]
        )
    return species_omega_e, stor


def neo2ql_stor2sqrtspol(neo2ql_input_file, stor):
    with h5py.File(neo2ql_input_file, 'r') as neo2ql:
        stor_i = np.asarray(neo2ql['boozer_s'], dtype=float).ravel()
        sqrtspol_i = np.asarray(neo2ql['rho_pol'], dtype=float).ravel()
    sqrtspol = np.interp(np.asarray(stor, dtype=float), stor_i, sqrtspol_i)
    return sqrtspol


def write_omega_e_to_mars_input(omega_e, sqrtspol, output_file='PROFWE.IN'):
    output_file = Path(output_file)
    omega_e = np.asarray(omega_e, dtype=float).ravel()
    sqrtspol = np.asarray(sqrtspol, dtype=float).ravel()
    number_of_surfaces = len(omega_e)
    type_of_radial_variable = 1  # 1 means sqrtspol for MARS
    header = f'{number_of_surfaces} {type_of_radial_variable}'
    np.savetxt(output_file, np.column_stack([sqrtspol, omega_e]), header=header, comments='')
