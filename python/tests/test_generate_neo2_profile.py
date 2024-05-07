#%% Standard modules
import os
import numpy as np
import h5py

# Custion modules
from neo2_mars import write_neo2ql_inputs_to_hdf5
from neo2_mars import get_coulomb_logarithm, get_kappa, derivative

test_mars_dir = '/proj/plasma/DATA/DEMO/MARS/MARSQ_INPUTS_KNTV21_NEO2profs_RUN/'
test_output_dir = '/tmp/'
test_hdf5_filename = os.path.join(test_output_dir, 'test_neo2ql_profiles.in')
test_src = {
    'sqrtspol': {'filename': os.path.join(test_output_dir, 'sqrtstor.dat'), 'column': 0},
    'sqrtstor': {'filename': os.path.join(test_output_dir, 'sqrtstor.dat'), 'column': 1},
    'ne': {'filename': os.path.join(test_output_dir, 'ne.dat'), 'column': 1},
    'Te': {'filename': os.path.join(test_output_dir, 'Te.dat'), 'column': 1},
    'Ti': {'filename': os.path.join(test_output_dir, 'Ti.dat'), 'column': 1},
    'vrot': {'filename': os.path.join(test_output_dir, 'vrot.dat'), 'column': 1}   
}

test_inputs = {
    '/num_radial_pts': 10,
    '/num_species': 2,
    '/species_tag': np.array([1, 2], dtype=np.int32),
    '/species_def': np.tile(np.array([[-1,1],[1,10]]), 10),
    '/boozer_s': np.linspace(0, 1, 10),
    '/rho_pol': np.sqrt(np.linspace(0, 1, 10)),
    '/Vphi': np.linspace(10,1,10),
    '/species_tag_Vphi': 1,
    '/isw_Vphi_loc': 0,
    '/rel_stages': 1*np.ones((10,)),
    '/T_prof': np.column_stack([np.linspace(20,2,10), np.linspace(30,3,10)]),
    '/dT_ov_ds_prof': np.column_stack([2*np.ones((10,)), 3*np.ones((10,))]),
    '/n_prof': np.column_stack([np.linspace(40,4,10), np.linspace(50,5,10)]),
    '/dn_ov_ds_prof': np.column_stack([4*np.ones((10,)), 5*np.ones((10,))]),
    '/kappa_prof': np.column_stack([np.linspace(60,6,10), np.linspace(70,7,10)]), 
}

test_inputs_types = {
    '/num_radial_pts': 'int32',
    '/num_species': 'int32',
    '/species_tag': 'int32',
    '/species_def': 'float64',
    '/boozer_s': 'float64',
    '/rho_pol': 'float64',
    '/Vphi': 'float64',
    '/species_tag_Vphi': 'int32',
    '/isw_Vphi_loc': 'int32',
    '/rel_stages': 'int32',
    '/T_prof': 'float64',
    '/dT_ov_ds_prof': 'float64',
    '/n_prof': 'float64',
    '/dn_ov_ds_prof': 'float64',
    '/kappa_prof': 'float64'
}

test_inputs_attributes = {
    '/species_def': {'unit': 'e; g'},
    '/boozer_s': {'unit': '1'},
    '/rho_pol': {'unit': '1'},
    '/Vphi': {'unit': 'rad / s'},
    '/T_prof': {'unit': 'erg'},
    '/n_prof': {'unit': '1 / cm^3'}
}

def test_write_neo2ql_inputs_to_hdf5():
    write_neo2ql_inputs_to_hdf5(test_hdf5_filename, test_inputs)
    with h5py.File(test_hdf5_filename, 'r') as hdf5_file:
        assert are_all_quantities_present(hdf5_file, test_inputs)
        assert are_data_types_correct(hdf5_file, test_inputs_types)
        assert are_values_correct(hdf5_file, test_inputs)
        assert are_attributes_correct(hdf5_file, test_inputs_attributes)

def are_all_quantities_present(hdf5_file, inputs):
    for key in inputs.keys():
        if key not in hdf5_file:
            return False
    return True

def are_data_types_correct(hdf5_file, inputs_type):
    for key, dtype in inputs_type.items():
        if hdf5_file[key].dtype != dtype:
            return False
    return True

def are_values_correct(hdf5_file, inputs):
    for key, value in inputs.items():
        if not np.allclose(hdf5_file[key][()], value):
            return False
    return True

def are_attributes_correct(hdf5_file, inputs_attributes):
    for key, attributes in inputs_attributes.items():
        for attr_key, attr_value in attributes.items():
            if attr_key not in hdf5_file[key].attrs:
                return False
            if hdf5_file[key].attrs[attr_key] != attr_value:
                return False
    return True

def test_derivative():
    x = np.linspace(0, 2*np.pi, 60)
    y = np.cos(x)**2
    dy_dx = derivative(x, y)
    control_dy_dx = -2*np.cos(x)*np.sin(x)
    print(np.max(np.abs(dy_dx-control_dy_dx)))
    assert np.allclose(dy_dx, control_dy_dx, atol=1e-2)

def test_derivative_visual_check():
    import matplotlib.pyplot as plt
    x = np.linspace(0, 2*np.pi, 60)
    y = np.cos(x)**2
    dy_dx = derivative(x, y)
    control_dy_dx = -2*np.cos(x)*np.sin(x)
    plt.plot(x, control_dy_dx - dy_dx)
    plt.title('Difference between control and calculated derivative')
    plt.show()

def test_get_kappa():
    n_CGS = 1.22931e+14
    T_CGS = 6.6154e-08
    charge_CGS = -1 * 1.60217662e-19 * 2.99792458e8 * 10
    kappa = get_kappa(n_CGS, T_CGS, charge_CGS)
    control_kappa = 1.39308e-07
    assert np.isclose(kappa, control_kappa, rtol=1e-4)

def test_coulomb_logarithm():
    n_SI = 1e19
    T_eV = 1e6
    log_Lambda_from_SI = get_coulomb_logarithm_from_SI_input(n_SI, T_eV)
    n_CGS = n_SI * 1e-6
    T_CGS = T_eV * 1.60217662e-19 * 1e7
    log_Lambda = get_coulomb_logarithm(n_CGS, T_CGS)
    assert np.isclose(log_Lambda_from_SI, log_Lambda, rtol=1e-4)

def get_coulomb_logarithm_from_SI_input(n_SI, T_eV):
    # Coulomb logarithm (set as species-independent - see E A Belli and J Candy PPCF 54 015015 (2012))
    log_Lambda = 39.1 - 1.15 * np.log10(n_SI) + 2.3 * np.log10(T_eV * 1e-3)
    return log_Lambda

    

if __name__ == '__main__':
    test_write_neo2ql_inputs_to_hdf5()
    test_derivative()
    test_derivative_visual_check()
    test_coulomb_logarithm()
    test_get_kappa()
    print('All tests passed.')