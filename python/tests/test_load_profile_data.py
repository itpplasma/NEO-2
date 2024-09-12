#%% Standard modules
import os
import numpy as np
import copy

# Custom modules
from neo2_mars import write_neo2_input_profiles_from_mars

# Modules to test
from neo2_ql import load_cgs_profiles_and_interp
from neo2_ql import convert_units_from_norm_to_si, convert_units_from_si_to_cgs

mars_dir = '/proj/plasma/DATA/DEMO/MARS/MARSQ_INPUTS_KNTV21_NEO2profs_RUN/'
output_dir = '/tmp/'
profiles_src = {
    'sqrtspol': {'filename': os.path.join(output_dir, 'sqrtstor.dat'), 'column': 0},
    'sqrtstor': {'filename': os.path.join(output_dir, 'sqrtstor.dat'), 'column': 1},
    'n': {'filename': os.path.join(output_dir, 'n.dat'), 'column': 1},
    'T': {'filename': os.path.join(output_dir, 'T.dat'), 'column': 1},
    'vrot': {'filename': os.path.join(output_dir, 'vrot.dat'), 'column': 1}   
}
multispec_profiles_src = copy.deepcopy(profiles_src)
multispec_profiles_src['n']['column'] = [1,2]
multispec_profiles_src['T']['column'] = [1,2]
mars_profiles_src = copy.deepcopy(multispec_profiles_src)

def test_write_and_read_compatiblity():
    write_neo2_input_profiles_from_mars(mars_dir, output_dir)
    profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(mars_profiles_src)

def test_output_type():
    write_neo2_input_profiles_from_mars(mars_dir, output_dir)
    profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(mars_profiles_src, interp_config={'grid':'sqrtspol'})
    assert is_profile_types_correct(profiles, sqrtspol, sqrtstor)
    profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(mars_profiles_src, interp_config={'grid':'spol'})
    assert is_profile_types_correct(profiles, sqrtspol, sqrtstor)
    profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(mars_profiles_src, interp_config={'grid':'sqrtstor'})
    assert is_profile_types_correct(profiles, sqrtspol, sqrtstor)

def is_profile_types_correct(profiles, sqrtspol, sqrtstor):
    bool = True
    bool = bool and type(sqrtspol) == np.ndarray
    bool = bool and type(sqrtstor) == np.ndarray
    for profile in profiles:
        bool = bool and type(profiles[profile]) == np.ndarray
    return bool

def test_output_shape():
    write_neo2_input_profiles_from_mars(mars_dir, output_dir)
    profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(mars_profiles_src, interp_config={'grid':'sqrtspol'})
    assert is_profile_shapes_correct(profiles, sqrtspol, sqrtstor)
    profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(mars_profiles_src, interp_config={'grid':'spol'})
    assert is_profile_shapes_correct(profiles, sqrtspol, sqrtstor)
    profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(mars_profiles_src, interp_config={'grid':'sqrtstor'})
    assert is_profile_shapes_correct(profiles, sqrtspol, sqrtstor)

def is_profile_shapes_correct(profiles, sqrtspol, sqrtstor):
    bool = True
    n_s = len(sqrtspol)
    bool = bool and sqrtspol.shape == (n_s,)
    bool = bool and sqrtstor.shape == (n_s,)
    for profile in profiles:
        if profile == 'T' or profile == 'n':
            bool = bool and profiles[profile].shape == (2,n_s)
        else:
            bool = bool and profiles[profile].shape == (n_s,)
    return bool

def test_equidistant_grid():
    write_neo2_input_profiles_from_mars(mars_dir, output_dir)
    _, sqrtspol, _ = load_cgs_profiles_and_interp(mars_profiles_src, interp_config={'grid':'sqrtspol'})
    assert is_equidistant(sqrtspol)
    _, _, sqrtstor = load_cgs_profiles_and_interp(mars_profiles_src, interp_config={'grid':'sqrtstor'})
    assert is_equidistant(sqrtstor)
    _, sqrtspol, _ = load_cgs_profiles_and_interp(mars_profiles_src, interp_config={'grid':'spol'})
    assert is_equidistant(sqrtspol**2)
    _, _, sqrtstor = load_cgs_profiles_and_interp(mars_profiles_src, interp_config={'grid':'stor'})
    assert is_equidistant(sqrtstor**2)

def test_interpolation():
    test_equidist_sqrtspol_interpolation()
    test_equidist_spol_interpolation()
    test_equidist_sqrtstor_interpolation()
    test_equidist_stor_interpolation()

def test_equidist_sqrtspol_interpolation():
    write_trial_profiles(trial_profiles_sqrtspol, output_dir)
    trial_profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(profiles_src, interp_config={'grid':'sqrtspol'})
    assert is_equidistant(sqrtspol)
    assert not is_equidistant(sqrtstor)
    trial_profiles = backconvert_trial_profiles_from_cgs(trial_profiles)
    assert is_equidistant(trial_profiles['n']/1)
    assert is_equidistant(trial_profiles['T']/3)
    assert is_equidistant(trial_profiles['vrot']/5)

def test_equidist_spol_interpolation():
    write_trial_profiles(trial_profiles_spol, output_dir)
    trial_profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(profiles_src, interp_config={'grid':'spol'})
    assert is_equidistant(sqrtspol**2)
    assert not is_equidistant(sqrtstor)
    trial_profiles = backconvert_trial_profiles_from_cgs(trial_profiles)
    assert is_equidistant((trial_profiles['n']/1)**2)
    assert is_equidistant((trial_profiles['T']/3)**2)
    assert is_equidistant((trial_profiles['vrot']/5)**2)

def test_equidist_sqrtstor_interpolation():
    write_trial_profiles(trial_profiles_sqrtstor, output_dir)
    trial_profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(profiles_src, interp_config={'grid':'sqrtstor'})
    assert not is_equidistant(sqrtspol)
    assert is_equidistant(sqrtstor)
    trial_profiles = backconvert_trial_profiles_from_cgs(trial_profiles)
    assert is_equidistant(trial_sqrtspol2sqrtstor(trial_profiles['n']/1))
    assert is_equidistant(trial_sqrtspol2sqrtstor(trial_profiles['T']/3))
    assert is_equidistant(trial_sqrtspol2sqrtstor(trial_profiles['vrot']/5))

def test_equidist_stor_interpolation():

    def trial_profiles_stor(stor): # Note that this is not consistent with the
        sqrtspol = stor**(3/1)     # conversion sqrtspol-sqrtstor used usually
        profiles, sqrtspol, _ = trial_profiles_sqrtspol(sqrtspol)
        return profiles, sqrtspol, np.sqrt(stor)

    def trial_sqrtspol2stor(sqrtspol): # We do take this inconsistency to keep
        return sqrtspol**(1/3)         # stor-sqrtspol exact cubic interpolatable

    write_trial_profiles(trial_profiles_stor, output_dir)
    trial_profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(profiles_src, interp_config={'grid':'stor'})
    assert not is_equidistant(sqrtspol)
    assert is_equidistant(sqrtstor**2)
    trial_profiles = backconvert_trial_profiles_from_cgs(trial_profiles)
    assert is_equidistant(trial_sqrtspol2stor(trial_profiles['n']/1))
    assert is_equidistant(trial_sqrtspol2stor(trial_profiles['T']/3))
    assert is_equidistant(trial_sqrtspol2stor(trial_profiles['vrot']/5))

def test_multispec_interpolation():
    write_trial_profiles(trial_multispec_profiles_sqrtspol, output_dir)
    profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(multispec_profiles_src, 
                                                      interp_config={'grid':'sqrtstor'})
    assert not is_equidistant(sqrtspol)
    assert is_equidistant(sqrtstor)
    profiles = backconvert_trial_profiles_from_cgs(profiles)
    assert is_equidistant(trial_sqrtspol2sqrtstor(profiles['n'][0]/1))
    assert is_equidistant(trial_sqrtspol2sqrtstor(profiles['n'][1]/2))
    assert is_equidistant(trial_sqrtspol2sqrtstor(profiles['T'][0]/3))
    assert is_equidistant(trial_sqrtspol2sqrtstor(profiles['T'][1]/4))
    assert is_equidistant(trial_sqrtspol2sqrtstor(profiles['vrot']/5))

def trial_multispec_profiles_sqrtspol(sqrtspol):
    profiles = {
        'n': np.array([sqrtspol*1, sqrtspol*2]),
        'T': np.array([sqrtspol*3, sqrtspol*4]),
        'vrot': sqrtspol*5
    }
    sqrtstor = trial_sqrtspol2sqrtstor(sqrtspol)
    return profiles, sqrtspol, sqrtstor

def is_equidistant(x):
    return np.allclose(np.diff(x), np.diff(x)[0])

def test_unit_conversion():
    write_trial_profiles(trial_profiles_sqrtspol, output_dir)
    trial_profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(profiles_src, interp_config={'grid':'sqrtspol'})
    trial_profiles_cgs = trial_profiles_sqrtspol_cgs(sqrtspol)
    assert are_same_profiles(trial_profiles, trial_profiles_cgs)

def write_trial_profiles(trial_profiles, output_dir):
    profiles, sqrtspol, sqrtstor = trial_profiles(np.linspace(0, 1, 11))
    for profile in profiles:
        np.savetxt(os.path.join(output_dir, profile+'.dat'), np.vstack([sqrtspol, profiles[profile]]).T)
    np.savetxt(os.path.join(output_dir, 'sqrtstor.dat'), np.vstack([sqrtspol, sqrtstor]).T)

def trial_profiles_spol(spol):
    sqrtspol = np.sqrt(spol)
    return trial_profiles_sqrtspol(sqrtspol)

def trial_profiles_sqrtstor(sqrtstor):
    sqrtspol = trial_sqrtstor2sqrtspol(sqrtstor)
    return trial_profiles_sqrtspol(sqrtspol)

def trial_profiles_stor(stor):
    sqrtstor = np.sqrt(stor)
    return trial_profiles_sqrtstor(sqrtstor)

def trial_sqrtstor2sqrtspol(sqrtstor):
    sqrtspol = 0.5*(sqrtstor + sqrtstor**(2/1))
    return sqrtspol

def trial_sqrtspol2sqrtstor(sqrtspol):
    sqrtstor = -0.5 + np.sqrt(0.25 + 2*sqrtspol)
    return sqrtstor

def trial_profiles_sqrtspol_cgs(sqrtspol):
    profiles, _, _ = trial_profiles_sqrtspol(sqrtspol)
    profiles = convert_trial_profiles_to_cgs(profiles)
    return profiles

def trial_profiles_sqrtspol(sqrtspol):
    profiles = {
        'n': sqrtspol*1,
        'T': sqrtspol*3,
        'vrot': sqrtspol*5,
    }
    sqrtstor = trial_sqrtspol2sqrtstor(sqrtspol)
    return profiles, sqrtspol, sqrtstor

def convert_trial_profiles_to_cgs(profiles):
    profiles['n'] *= 1e-6
    profiles['T'] *= 1.60217662e-19 * 1e7
    profiles['vrot'] *= 1.0
    return profiles

def backconvert_trial_profiles_from_cgs(profiles):
    profiles['n'] /= 1e-6
    profiles['T'] /= 1.60217662e-19 * 1e7
    profiles['vrot'] /= 1.0
    return profiles

def are_same_profiles(profiles1, profiles2):
    bool = True
    for profile in profiles1:
        bool = bool and np.allclose(profiles1[profile], profiles2[profile])
    return bool

def test_interpolation_visual_check():
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(4, 4, figsize=(12, 6))
    write_trial_profiles(trial_profiles_sqrtspol, output_dir)
    trial_profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(profiles_src, interp_config={'grid':'sqrtspol'})
    plot_profiles(ax[:,0], trial_profiles, sqrtspol)
    write_trial_profiles(trial_profiles_spol, output_dir)
    trial_profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(profiles_src, interp_config={'grid':'spol'})
    plot_profiles(ax[:,1], trial_profiles, sqrtspol)
    write_trial_profiles(trial_profiles_sqrtstor, output_dir)
    trial_profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(profiles_src, interp_config={'grid':'sqrtstor'})
    plot_profiles(ax[:,2], trial_profiles, sqrtspol)
    write_trial_profiles(trial_profiles_stor, output_dir)
    trial_profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(profiles_src, interp_config={'grid':'stor'})
    plot_profiles(ax[:,3], trial_profiles, sqrtspol)
    plt.subplots_adjust(hspace=0.8, wspace=0.3)
    ax[0,0].set_title('equidist sqrtspol grid')
    ax[0,1].set_title('equidist spol grid')
    ax[0,2].set_title('equidist sqrtstor grid')
    ax[0,3].set_title('equidist stor grid')
    plt.show()

def plot_profiles(ax, profiles, sqrtspol):
    profiles = backconvert_trial_profiles_from_cgs(profiles)
    xs = [sqrtspol, sqrtspol**2, trial_sqrtspol2sqrtstor(sqrtspol), trial_sqrtspol2sqrtstor(sqrtspol)**2]
    xnames = ['sqrtspol', 'spol', 'sqrtstor', 'stor']
    for i,x in enumerate(xs):
        for profile in profiles:
            ax[i].plot(x, profiles[profile],'o', label=profile)
        ax[i].set_xlabel(xnames[i] + ' [1]')
        ax[i].set_ylabel('profiles [1]')
    

if __name__ == '__main__':
    test_write_and_read_compatiblity()
    test_output_type()
    test_output_shape()
    test_equidistant_grid()
    test_interpolation()
    test_multispec_interpolation()
    test_unit_conversion()
    print('All tests passed')
    test_interpolation_visual_check()