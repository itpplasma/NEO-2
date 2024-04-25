#%% Standard modules
import os

# Custion modules
from neo2_mars import write_neo2_input_profile_from_mars

# Modules to test
from neo2_mars import get_profiles_over_equidist_grid

test_mars_dir = '/proj/plasma/DATA/DEMO/MARS/MARSQ_INPUTS_KNTV21_NEO2profs_RUN/'
test_output_dir = '/tmp/'
test_src = {
    'sqrtspol': {'filename': os.path.join(test_output_dir, 'sqrtstor.dat'), 'column': 0},
    'sqrtstor': {'filename': os.path.join(test_output_dir, 'sqrtstor.dat'), 'column': 1},
    'ne': {'filename': os.path.join(test_output_dir, 'ne.dat'), 'column': 1},
    'Te': {'filename': os.path.join(test_output_dir, 'Te.dat'), 'column': 1},
    'Ti': {'filename': os.path.join(test_output_dir, 'Ti.dat'), 'column': 1},
    'vrot': {'filename': os.path.join(test_output_dir, 'vrot.dat'), 'column': 1}   
}

def test_equidistant_grid():
    write_neo2_input_profile_from_mars(test_mars_dir, test_output_dir)
    profiles = get_profiles_over_equidist_grid(test_src, 'sqrtspol')
    assert is_equidistant(profiles['sqrtspol'])
    profiles = get_profiles_over_equidist_grid(test_src, 'sqrtstor')
    assert is_equidistant(profiles['sqrtstor'])
    profiles = get_profiles_over_equidist_grid(test_src, 'spol')
    assert is_equidistant(profiles['sqrtspol']**2)

def is_equidistant(x):
    return np.allclose(np.diff(x), np.diff(x)[0])

if __name__ == '__main__':
    test_equidistant_grid()