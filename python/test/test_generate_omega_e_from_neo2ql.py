#%% Stadart imports
import h5py
import matplotlib.pyplot as plt
import os
import numpy as np
import pytest

# Homebrew imports
from neo2_mars import mars_sqrtspol2stor
from neo2_mars import mars_sqrtspol2sqrtstor

# Modules to test
from neo2_mars import get_omega_e_from_neo2ql
from neo2_mars import neo2ql_stor2sqrtspol
from neo2_mars import write_omega_e_to_mars_input
from neo2_mars import generate_omega_e_for_mars

test_neo2ql_output_file = "/itp/MooseFS/grassl_g/comparison_vary_coilwidth_MARS_NEO2/more_points_equidistant_sqrtspol_run_000/neo2_multispecies_out.h5"
test_neo2ql_input_file = "/itp/MooseFS/grassl_g/comparison_vary_coilwidth_MARS_NEO2/more_points_equidistant_sqrtspol_run_000/multi_spec_demo.in"
test_mars_dir = "/proj/plasma/DATA/DEMO/MARS/MARSQ_OUTPUTS_100kAt_dBkinetic_NTVkinetic_NEO2profs_KEYTORQ_1"
test_neo2ql_input_file = "/temp/grassl_g/comparison_vary_coilwidth_MARS_NEO2/correct_electric_rotation_profile_run_000/multi_spec_demo.in"
test_neo2ql_output_file = "/temp/grassl_g/comparison_vary_coilwidth_MARS_NEO2/correct_electric_rotation_profile_run_000/neo2_multispecies_out.h5"
faulty_mars_dir = "/proj/plasma/DATA/DEMO/MARS/script_get_omega_e_from_NEO_2_result"


def requires(*paths):
    """Skip a test that compares against cluster data when the data is absent."""
    missing = [path for path in paths if not os.path.exists(path)]
    return pytest.mark.skipif(bool(missing), reason=f"missing test data: {missing}")


def _write_h5(path, **datasets):
    with h5py.File(path, "w") as handle:
        for name, value in datasets.items():
            handle.create_dataset(name, data=value)


# Two species chosen so that v_th = sqrt(2 T / m) is exact: electrons
# 2e9 cm/s, ions 1e8 cm/s. With MtOvR = Om_tE / v_th the species Om_tE are
# 2e4 rad/s and 3e4 rad/s; MARS gets the first species with opposite sign.
T_SPEC = np.array([5.0e-9, 2.0e-9])
M_SPEC = np.array([2.5e-27, 4.0e-25])
MTOVR = np.array([1.0e-5, 3.0e-4])
OMEGA_E = np.array([2.0e4, 3.0e4])


def test_get_omega_e_from_neo2ql_single_surface(tmp_path):
    # A single-surface run stores boozer_s as a scalar and per-species data
    # as 1D arrays of length num_species.
    output_file = tmp_path / "neo2_multispecies_out.h5"
    _write_h5(output_file, boozer_s=0.5, MtOvR=MTOVR, T_spec=T_SPEC, m_spec=M_SPEC)
    omega_e, stor = get_omega_e_from_neo2ql(output_file)
    assert omega_e.shape == (1, 2)
    assert np.allclose(stor, [0.5])
    assert np.allclose(omega_e[0], OMEGA_E, rtol=1e-12)


def test_get_omega_e_from_neo2ql_multiple_surfaces(tmp_path):
    # Collected multi-surface output: (num_surfaces, num_species).
    output_file = tmp_path / "neo2_multispecies_out.h5"
    scale = np.array([1.0, 2.0, 3.0])[:, np.newaxis]
    _write_h5(output_file, boozer_s=np.array([0.2, 0.5, 0.8]),
              MtOvR=scale * MTOVR, T_spec=np.tile(T_SPEC, (3, 1)), m_spec=M_SPEC)
    omega_e, stor = get_omega_e_from_neo2ql(output_file)
    assert np.allclose(stor, [0.2, 0.5, 0.8])
    assert np.allclose(omega_e, scale * OMEGA_E, rtol=1e-12)


def test_generate_omega_e_for_mars_single_surface(tmp_path):
    # rho_pol at the surface is interpolated from the multispec input; the
    # PROFWE.IN header is "<number of surfaces> 1" (1: sqrt(s_pol) grid).
    output_file = tmp_path / "neo2_multispecies_out.h5"
    input_file = tmp_path / "multi_spec.in"
    profwe_file = tmp_path / "PROFWE.IN"
    _write_h5(output_file, boozer_s=0.5, MtOvR=MTOVR, T_spec=T_SPEC, m_spec=M_SPEC)
    _write_h5(input_file, boozer_s=np.array([0.0, 0.25, 1.0]),
              rho_pol=np.array([0.0, 0.5, 1.0]))
    generate_omega_e_for_mars(output_file, input_file, output_file=profwe_file)
    data = np.loadtxt(profwe_file)
    assert np.array_equal(data[0], [1, 1])
    assert np.allclose(data[1:], [[2.0 / 3.0, -OMEGA_E[0]]], rtol=1e-12)



@requires(test_neo2ql_input_file, test_mars_dir)
def test_neo2_stor2sqrtspol():
    neo2ql = h5py.File(test_neo2ql_input_file, "r")
    stor = np.array(neo2ql["boozer_s"])
    sqrtspol = neo2ql_stor2sqrtspol(test_neo2ql_input_file, stor)
    assert np.allclose(stor, mars_sqrtspol2stor(test_mars_dir, sqrtspol), atol=1e-6)
    assert np.allclose(sqrtspol, np.array(neo2ql["rho_pol"]))

def test_write_omega_e_to_mars_input(tmp_path):
    sqrtspol = np.linspace(0,10)
    omega_e = 10*sqrtspol
    write_omega_e_to_mars_input(omega_e, sqrtspol, output_file=tmp_path / "PROFWE.IN")
    data_read = np.loadtxt(tmp_path / "PROFWE.IN")
    header = data_read[0]
    sqrtspol_read = data_read[1:,0]
    omega_e_read = data_read[1:,1]
    assert header[0] == len(sqrtspol)
    assert header[1] == 1 # MARS reads the profile in sqrtspol then
    assert np.allclose(sqrtspol_read, sqrtspol)
    assert np.allclose(omega_e_read, omega_e)

@requires(test_neo2ql_output_file, test_neo2ql_input_file, "/proj/plasma/DATA/DEMO/MARS/PROFWE")
def test_generate_omega_e_for_mars_visual_check():
    generate_omega_e_for_mars(test_neo2ql_output_file, test_neo2ql_input_file)
    omega_e_file = "PROFWE.IN"
    mars_dir = "/proj/plasma/DATA/DEMO/MARS/PROFWE"
    omegate_mars, sqrtspol_mars = get_mars_omega_e(mars_dir)
    omega_e = np.loadtxt(omega_e_file, skiprows=1)
    sqrtspol = omega_e[:,0]
    omega_e = omega_e[:,1]
    plt.figure()
    plt.plot(sqrtspol, omega_e, '-ob', label='NEO-2-QL')
    plt.plot(sqrtspol_mars, omegate_mars, '--r', label='MARS')
    plt.xlabel(r"$\rho_\mathrm{{pol}}$ [1]")
    plt.ylabel(r"$\Omega_\mathrm{{E}}$ [1/s]")
    plt.xlim(0,1)
    plt.show()

@requires(faulty_mars_dir)
def test_get_omega_e_from_neo2ql_visual_check():
    faulty_neo2ql_output_file = "/proj/plasma/DATA/DEMO/MARS/script_get_omega_e_from_NEO_2_result/neo2_multispecies_out_22052023.h5"
    omegate_neo2ql, stor_neo2ql = get_omega_e_from_neo2ql(faulty_neo2ql_output_file)
    omegate_mars, sqrtspol_mars = get_mars_omega_e(faulty_mars_dir)
    stor_mars = sqrtspol_mars**2 # the faulty part was a mistranslation stor -> sqrtspol
    plt.figure()
    plt.plot(stor_neo2ql,omegate_neo2ql[:,0],'-b', label='NEO-2-QL species 1')
    plt.plot(stor_neo2ql,omegate_neo2ql[:,1],'--g', label='NEO-2-QL species 2')
    plt.plot(stor_mars,omegate_mars,'--r', label='MARS')
    plt.xlabel(r"$s_\mathrm{{tor}}$ [1]")
    plt.ylabel(r"$\Omega_\mathrm{{E}}$ [1/s]")
    plt.xlim(0,1)
    plt.ylim(-5e3,7.5e3)
    plt.title("Comparison of $\Omega_\mathrm{{E}}$ from NEO-2-QL and MARS \n post-/preprocessing (faulty profiles)")
    plt.legend()
    plt.show()

def get_mars_omega_e(mars_dir):
    omega_e_file = os.path.join(mars_dir, "PROFWE.IN")
    omega_e = np.loadtxt(omega_e_file, skiprows=1)
    sqrtspol = omega_e[10:,0]
    omega_e = np.array(omega_e[10:,1])
    return omega_e, sqrtspol

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
