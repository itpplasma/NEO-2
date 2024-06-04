#%% Standard modules
import numpy as np
import matplotlib.pyplot as plt
import os
import h5py

# Homebrew modules
from neo2_ql import get_neo2_ql_input_profiles
from neo2_ql import make_figure_neo2_ql_input_profiles, add_profiles_to_axes

# Modules to test
from neo2_mars import get_species_cgs_from_mars
from neo2_mars import generate_multispec_input_from_mars

test_mars_dir = '/proj/plasma/DATA/DEMO/MARS/MARSQ_INPUTS_KNTV21_NEO2profs_RUN/'
test_outpu_dir = '/tmp/'

def test_get_species_cgs_from_mars():
    Ze, Zi, me, mi = get_species_cgs_from_mars(test_mars_dir)
    assert Ze == -1
    assert Zi == 1
    assert np.isclose(me, 9.10938356e-28, rtol=1e-4, atol=0)
    assert np.isclose(mi, 3.343583719e-24, rtol=1e-3, atol=0) # off due MARS using not atomic unit, but proton mass

def test_generate_multispec_input_from_mars_init():
    output_file = os.path.join(test_outpu_dir, 'multi_spec.in')
    generate_multispec_input_from_mars(test_mars_dir, output_file, 
                                       number_of_surfaces=10, bounds=[0.0, 1.0])

def test_generate_multispec_input_from_mars_ion_tag_vrot():
    output_file = os.path.join(test_outpu_dir, 'multi_spec.in')
    generate_multispec_input_from_mars(test_mars_dir, output_file, 
                                       number_of_surfaces=10, bounds=[0.0, 1.0])
    with h5py.File(output_file, 'r') as f:
        Z, m = np.array(f['/species_def'])[:,int(f['/species_tag_Vphi'][()])-1,0]
        assert Z == 1
        assert np.isclose(m, 3.343583719e-24, rtol=1e-3, atol=0) # off due MARS using not atomic unit, but proton mass

def test_generate_multispec_input_from_mars_visual_check():
    output_file = os.path.join(test_outpu_dir, 'multi_spec.in')
    generate_multispec_input_from_mars(test_mars_dir, output_file, 
                                       number_of_surfaces=20, bounds=[0.0, 1.0])
    profiles, _ = get_neo2_ql_input_profiles(output_file)
    fig, axes = make_figure_neo2_ql_input_profiles(profiles, color='b', label='MARS converted')
    eps_2023_inputfile = "/proj/plasma/DATA/DEMO/teams/Equilibrium_DEMO2019_CHEASE/multi_spec_demo_version_1_2_equidistant_s.in"
    eps_2023_profiles, _ = get_neo2_ql_input_profiles(eps_2023_inputfile)
    add_profiles_to_axes(eps_2023_profiles, axes, color='k', label='EPS 2023')
    add_mars_profiles_to_axes(test_mars_dir, axes, color='r', label='MARS')
    plt.show()

def add_mars_profiles_to_axes(mars_dir, axes, color: str='r', label: str='MARS'):
    from neo2_mars import get_profiles_mars
    profiles_mars = get_profiles_mars(mars_dir)
    profiles = {}
    EV2ERG = 1.60217662e-12
    M3_TO_CM3 = 1e+6
    profiles['T_prof'] = {'x': profiles_mars['sqrtstor'][:, 1]**2, 
                          'y': np.row_stack([profiles_mars['Te'][:, 1] * EV2ERG, 
                                             profiles_mars['Ti'][:, 1] * EV2ERG])}
    profiles['dT_ov_ds_prof'] = {'x': np.array([]), 'y': np.array([])}
    profiles['n_prof'] = {'x': profiles_mars['sqrtstor'][:, 1]**2,
                          'y': np.row_stack([profiles_mars['ne'][:, 1] / M3_TO_CM3, 
                                             profiles_mars['ne'][:, 1] / M3_TO_CM3])}
    profiles['dn_ov_ds_prof'] = {'x': np.array([]), 'y': np.array([])}
    profiles['kappa_prof'] = {'x': np.array([]), 'y': np.array([])}
    profiles['Vphi'] = {'x': profiles_mars['sqrtstor'][:, 1]**2,
                        'y': profiles_mars['vrot'][:, 1]}
    profiles['rho_pol'] = {'x': profiles_mars['sqrtstor'][:, 1]**2,
                           'y': profiles_mars['sqrtstor'][:, 0]}
    add_profiles_to_axes(profiles, axes, color=color, label=label)
    
if __name__ == '__main__':
    test_get_species_cgs_from_mars()
    test_generate_multispec_input_from_mars_init()
    test_generate_multispec_input_from_mars_ion_tag_vrot()
    print('All tests passed.')
    #test_generate_multispec_input_from_mars_visual_check()