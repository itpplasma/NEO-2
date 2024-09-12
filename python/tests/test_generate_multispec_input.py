#%% Standard modules
import os
import numpy as np
import h5py

# Custom modules
from test_load_profile_data import write_trial_profiles
from test_load_profile_data import profiles_src

# modules to test
from neo2_ql import write_multispec_to_hdf5
from neo2_ql import get_coulomb_logarithm, get_kappa, derivative
from neo2_ql import generate_multispec_input
from neo2_ql import get_species_def_array

test_output_dir = '/tmp/'
test_hdf5_filename = os.path.join(test_output_dir, 'test_multispec.in')
test_config = {
    'hdf5_filename': test_hdf5_filename, 
    'species_tag_of_vrot': 1, 
    'Ze': -1, 
    'Zi': 1, 
    'me': 9.10938356e-28,
    'mi': 3.3436e-24,
    }

test_multispec = {
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

test_multispec_types = {
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

test_multispec_attributes = {
    '/species_def': {'unit': 'e; g'},
    '/boozer_s': {'unit': '1'},
    '/rho_pol': {'unit': '1'},
    '/Vphi': {'unit': 'rad / s'},
    '/T_prof': {'unit': 'erg'},
    '/n_prof': {'unit': '1 / cm^3'}
}

def test_write_multispec_to_hdf5():
    write_multispec_to_hdf5(test_hdf5_filename, test_multispec)
    with h5py.File(test_hdf5_filename, 'r') as hdf5_file:
        assert are_all_quantities_present(hdf5_file, test_multispec)
        assert are_data_types_correct(hdf5_file, test_multispec_types)
        assert are_values_correct(hdf5_file, test_multispec)
        assert are_attributes_correct(hdf5_file, test_multispec_attributes)

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
    assert np.allclose(dy_dx, control_dy_dx, atol=1e-2)

def test_get_kappa():
    n_cgs = 1.22931e+14
    T_cgs = 6.6154e-08
    charge_cgs = -1 * 1.60217662e-19 * 2.99792458e8 * 10
    kappa = get_kappa(n_cgs, T_cgs, charge_cgs)
    control_kappa = 1.39308e-07
    assert np.isclose(kappa, control_kappa, rtol=1e-4)

def test_coulomb_logarithm():
    n_si = 1e19
    T_eV = 1e6
    log_Lambda_from_si = get_coulomb_logarithm_from_si_input(n_si, T_eV)
    n_cgs = n_si * 1e-6
    T_cgs = T_eV * 1.60217662e-19 * 1e7
    log_Lambda = get_coulomb_logarithm(n_cgs, T_cgs)
    assert np.isclose(log_Lambda_from_si, log_Lambda, rtol=1e-4)

def test_generate_multispec_input_call():
    write_trial_profiles(trial_profiles_sqrtspol_non_zero, test_output_dir)
    generate_multispec_input(test_config, profiles_src, profiles_interp_config={})

def trial_profiles_sqrtspol_non_zero(sqrtspol):
    from test_load_profile_data import trial_profiles_sqrtspol
    profiles, sqrtspol, sqrtstor = trial_profiles_sqrtspol(sqrtspol)
    for profile in profiles:
        profiles[profile] += 1
    return profiles, sqrtspol, sqrtstor

def test_get_species_def_array():
    species_def_control = np.array([[[test_config['Ze']], [test_config['Zi']]], 
                                    [[test_config['me']], [test_config['mi']]]])
    species_def = get_species_def_array(test_config, 1)
    assert species_def.shape == species_def_control.shape
    assert np.allclose(species_def, species_def_control)
    species_def = get_species_def_array(test_config, 10)
    for i in range(10):
        assert species_def[:,:,i:i+1].shape == species_def_control.shape
        assert np.allclose(species_def[:,:,i:i+1], species_def_control)

def get_coulomb_logarithm_from_si_input(n_si, T_eV):
    # Coulomb logarithm (set as species-independent - see E A Belli and J Candy PPCF 54 015015 (2012))
    log_Lambda = 39.1 - 1.15 * np.log10(n_si) + 2.3 * np.log10(T_eV * 1e-3)
    return log_Lambda

def test_derivative_visual_check():
    import matplotlib.pyplot as plt
    x = np.linspace(0, 2*np.pi, 60)
    y = np.cos(x)**2
    dy_dx = derivative(x, y)
    control_dy_dx = -2*np.cos(x)*np.sin(x)
    plt.plot(x, control_dy_dx - dy_dx)
    plt.title('Difference between control and calculated derivative')
    plt.show()

    
if __name__ == '__main__':
    test_write_multispec_to_hdf5()
    test_derivative()
    test_coulomb_logarithm()
    test_get_kappa()
    test_generate_multispec_input_call()
    test_get_species_def_array()
    print('All tests passed.')
    test_derivative_visual_check()