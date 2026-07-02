#%% Standard modules
import numpy as np
import matplotlib.pyplot as plt
import os
import h5py

# Homebrew modules
from neo2_ql import get_neo2_ql_input_profiles
from neo2_ql import make_figure_neo2_ql_input_profiles, add_profiles_to_axes
from mars_test_utils import create_test_mars_dir

# Modules to test
from neo2_mars import get_species_cgs_from_mars
from neo2_mars import generate_multispec_input_from_mars

test_outpu_dir = '/tmp/'

def test_get_species_cgs_from_mars(tmp_path):
    test_mars_dir = create_test_mars_dir(tmp_path / 'mars')['mars_dir']
    Z, m = get_species_cgs_from_mars(test_mars_dir)
    assert Z[0] == -1
    assert Z[1] == 1
    assert np.isclose(m[0], 9.10938356e-28, rtol=1e-4, atol=0)
    assert np.isclose(m[1], 3.34524384738e-24, rtol=1e-12, atol=0)

def test_generate_multispec_input_from_mars_init(tmp_path):
    test_mars_dir = create_test_mars_dir(tmp_path / 'mars')['mars_dir']
    output_file = os.path.join(test_outpu_dir, 'multi_spec.in')
    generate_multispec_input_from_mars(test_mars_dir, output_file, 
                                       number_of_surfaces=10, bounds=[0.0, 1.0])

def test_generate_multispec_input_from_mars_ion_tag_vrot(tmp_path):
    test_mars_dir = create_test_mars_dir(tmp_path / 'mars')['mars_dir']
    output_file = os.path.join(test_outpu_dir, 'multi_spec.in')
    generate_multispec_input_from_mars(test_mars_dir, output_file, 
                                       number_of_surfaces=10, bounds=[0.0, 1.0])
    with h5py.File(output_file, 'r') as f:
        species_tag_vphi = int(np.asarray(f['/species_tag_Vphi']).reshape(-1)[0])
        Z, m = np.array(f['/species_def'])[:, species_tag_vphi - 1, 0]
        assert Z == 1
        assert np.isclose(m, 3.34524384738e-24, rtol=1e-12, atol=0)

def test_generate_multispec_input_from_mars_visual_check(tmp_path):
    fixture = create_test_mars_dir(tmp_path / 'mars')
    test_mars_dir = fixture['mars_dir']
    output_file = os.path.join(test_outpu_dir, 'multi_spec.in')
    generate_multispec_input_from_mars(test_mars_dir, output_file, 
                                       number_of_surfaces=20, bounds=[0.0, 1.0])
    profiles, _ = get_neo2_ql_input_profiles(output_file)
    fig, axes = make_figure_neo2_ql_input_profiles(profiles, color='b', label='MARS converted')
    add_mars_profiles_to_axes(test_mars_dir, axes, color='r', label='MARS')
    figure_file = os.path.join(test_outpu_dir, 'test_generate_multispec_input_from_mars_visual_check.png')
    fig.savefig(figure_file)
    plt.close(fig)
    assert os.path.exists(figure_file)

def add_mars_profiles_to_axes(mars_dir, axes, color: str='r', label: str='MARS'):
    from neo2_mars import get_profiles_mars
    profiles_mars = get_profiles_mars(mars_dir)
    profiles = {}
    EV2ERG = 1.60217662e-12
    M3_TO_CM3 = 1e+6
    profiles['T_prof'] = {'x': profiles_mars['sqrtstor'][:, 1]**2, 
                          'y': profiles_mars['T'][:, 1:].T * EV2ERG}
    profiles['dT_ov_ds_prof'] = {'x': np.array([]), 'y': np.array([])}
    profiles['n_prof'] = {'x': profiles_mars['sqrtstor'][:, 1]**2,
                          'y': profiles_mars['n'][:, 1:].T / M3_TO_CM3}
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
    test_generate_multispec_input_from_mars_visual_check()
