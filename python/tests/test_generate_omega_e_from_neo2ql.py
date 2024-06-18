#%% Stadart imports
import h5py
import matplotlib.pyplot as plt
import os
import numpy as np

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

def test_neo2_stor2sqrtspol():
    neo2ql = h5py.File(test_neo2ql_input_file, "r")
    stor = np.array(neo2ql["boozer_s"])
    sqrtspol = neo2ql_stor2sqrtspol(test_neo2ql_input_file, stor)
    assert np.allclose(stor, mars_sqrtspol2stor(test_mars_dir, sqrtspol), atol=1e-6)
    assert np.allclose(sqrtspol, np.array(neo2ql["rho_pol"]))

def test_write_omega_e_to_mars_input():
    sqrtspol = np.linspace(0,10)
    omega_e = 10*sqrtspol
    write_omega_e_to_mars_input(omega_e, sqrtspol)
    data_read = np.loadtxt("PROFWE.IN")
    header = data_read[0]
    sqrtspol_read = data_read[1:,0]
    omega_e_read = data_read[1:,1]
    assert header[0] == len(sqrtspol)
    assert header[1] == 1 # MARS reads the profile in sqrtspol then
    assert np.allclose(sqrtspol_read, sqrtspol)
    assert np.allclose(omega_e_read, omega_e)

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

def test_get_omega_e_from_neo2ql_visual_check():
    faulty_neo2ql_output_file = "/proj/plasma/DATA/DEMO/MARS/script_get_omega_e_from_NEO_2_result/neo2_multispecies_out_22052023.h5"
    faulty_mars_dir = "/proj/plasma/DATA/DEMO/MARS/script_get_omega_e_from_NEO_2_result"
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
    test_neo2_stor2sqrtspol()
    test_write_omega_e_to_mars_input()
    test_generate_omega_e_for_mars_visual_check()
    test_get_omega_e_from_neo2ql_visual_check()