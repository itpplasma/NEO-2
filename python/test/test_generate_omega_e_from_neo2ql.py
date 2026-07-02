from pathlib import Path

import h5py
import numpy as np

from neo2_mars import (
    generate_omega_e_for_mars,
    get_omega_e_from_neo2ql,
    neo2ql_stor2sqrtspol,
    write_omega_e_to_mars_input,
)


FIXTURE_DIR = Path(__file__).with_name('data')
AXISYMMETRIC_OUTPUT_FIXTURE = FIXTURE_DIR / 'neo2_ql_axisymmetric_multispecies_out.h5'


def _write_mapping_input(path, boozer_s, rho_pol):
    with h5py.File(path, 'w') as handle:
        handle.create_dataset('boozer_s', data=np.asarray(boozer_s, dtype=float))
        handle.create_dataset('rho_pol', data=np.asarray(rho_pol, dtype=float))


def test_neo2_stor2sqrtspol():
    input_file = Path('/tmp/test_generate_omega_e_mapping.h5')
    boozer_s = np.array([0.0, 0.25, 1.0])
    rho_pol = np.array([0.0, 0.5, 1.0])
    _write_mapping_input(input_file, boozer_s, rho_pol)

    stor = np.array([0.0, 0.125, 0.25, 0.625, 1.0])
    sqrtspol = neo2ql_stor2sqrtspol(input_file, stor)

    assert np.allclose(sqrtspol, np.array([0.0, 0.25, 0.5, 0.75, 1.0]))


def test_write_omega_e_to_mars_input(tmp_path):
    sqrtspol = np.linspace(0.0, 1.0, 11)
    omega_e = 10.0 * sqrtspol
    output_file = tmp_path / 'PROFWE.IN'

    write_omega_e_to_mars_input(omega_e, sqrtspol, output_file=output_file)

    data_read = np.loadtxt(output_file)
    header = data_read[0]
    sqrtspol_read = data_read[1:, 0]
    omega_e_read = data_read[1:, 1]

    assert header[0] == len(sqrtspol)
    assert header[1] == 1
    assert np.allclose(sqrtspol_read, sqrtspol)
    assert np.allclose(omega_e_read, omega_e)


def test_get_omega_e_from_neo2ql_handles_single_surface_output():
    species_omega_e, stor = get_omega_e_from_neo2ql(AXISYMMETRIC_OUTPUT_FIXTURE)

    assert species_omega_e.shape == (1, 2)
    assert np.allclose(stor, np.array([0.5082424242424243]))
    assert np.allclose(species_omega_e[0], np.array([8166.0852282, 8166.0852282]))


def test_generate_omega_e_for_mars_creates_profwe_from_single_surface_output(tmp_path):
    input_file = tmp_path / 'multi_spec_demo.in'
    output_file = tmp_path / 'PROFWE.IN'

    _write_mapping_input(
        input_file,
        boozer_s=np.array([0.5082424242424243]),
        rho_pol=np.array([0.71291109]),
    )

    generate_omega_e_for_mars(
        AXISYMMETRIC_OUTPUT_FIXTURE,
        input_file,
        output_file=output_file,
    )

    data_read = np.loadtxt(output_file)
    header = data_read[0]
    profile = data_read[1:]

    assert header[0] == 1
    assert header[1] == 1
    assert np.allclose(profile[:, 0], np.array([0.71291109]))
    assert np.allclose(profile[:, 1], np.array([-8166.0852282]))
