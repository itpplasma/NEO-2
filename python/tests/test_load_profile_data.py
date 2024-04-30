#%% Standard modules
import os
import numpy as np

# Custion modules
from neo2_mars import write_neo2_input_profile_from_mars

# Modules to test
from neo2_mars import load_profiles_and_interp

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

def test_write_and_read_compatiblity():
    write_neo2_input_profile_from_mars(test_mars_dir, test_output_dir)
    profiles, sqrtspol, sqrtstor = load_profiles_and_interp(test_src, 'sqrtspol')

def test_equidistant_grid():
    write_neo2_input_profile_from_mars(test_mars_dir, test_output_dir)
    _, sqrtspol, _ = load_profiles_and_interp(test_src, 'sqrtspol')
    assert is_equidistant(sqrtspol)
    _, _, sqrtstor = load_profiles_and_interp(test_src, 'sqrtstor')
    assert is_equidistant(sqrtstor)
    _, sqrtspol, _ = load_profiles_and_interp(test_src, 'spol')
    assert is_equidistant(sqrtspol**2)

def is_equidistant(x):
    return np.allclose(np.diff(x), np.diff(x)[0])

def test_interpolation():
    test_equidist_sqrtspol_interpolation()
    test_equidist_spol_interpolation()
    test_equidist_sqrtstor_interpolation()

def test_equidist_sqrtspol_interpolation():
    write_trial_profiles(trial_profiles_sqrtspol, test_output_dir)
    trial_profiles, sqrtspol, sqrtstor = load_profiles_and_interp(test_src, 'sqrtspol')
    assert is_equidistant(sqrtspol)
    assert not is_equidistant(sqrtstor)
    assert is_equidistant(trial_profiles['ne'])
    assert is_equidistant(trial_profiles['Te'])
    assert is_equidistant(trial_profiles['Ti'])
    assert is_equidistant(trial_profiles['vrot'])

def test_equidist_spol_interpolation():
    write_trial_profiles(trial_profiles_spol, test_output_dir)
    trial_profiles, sqrtspol, sqrtstor = load_profiles_and_interp(test_src, 'spol')
    assert is_equidistant(sqrtspol**2)
    assert not is_equidistant(sqrtstor)
    assert is_equidistant((trial_profiles['ne']/1)**2)
    assert is_equidistant((trial_profiles['Te']/2)**2)
    assert is_equidistant((trial_profiles['Ti']/3)**2)
    assert is_equidistant((trial_profiles['vrot']/4)**2)

def test_equidist_sqrtstor_interpolation():
    write_trial_profiles(trial_profiles_sqrtstor, test_output_dir)
    trial_profiles, sqrtspol, sqrtstor = load_profiles_and_interp(test_src, 'sqrtstor')
    assert not is_equidistant(sqrtspol)
    assert is_equidistant(sqrtstor)
    assert is_equidistant(trial_sqrtspol2sqrtstor(trial_profiles['ne']/1))
    assert is_equidistant(trial_sqrtspol2sqrtstor(trial_profiles['Te']/2))
    assert is_equidistant(trial_sqrtspol2sqrtstor(trial_profiles['Ti']/3))
    assert is_equidistant(trial_sqrtspol2sqrtstor(trial_profiles['vrot']/4))

def write_trial_profiles(trial_profiles, output_dir):
    profiles, sqrtspol, sqrtstor = trial_profiles(np.linspace(0, 1, 11))
    for profile in profiles:
        np.savetxt(os.path.join(output_dir, profile+'.dat'), np.array([sqrtspol, profiles[profile]]).T)
    np.savetxt(os.path.join(output_dir, 'sqrtstor.dat'), np.array([sqrtspol, sqrtstor]).T)

def trial_profiles_sqrtspol(sqrtspol):
    profiles = {
        'ne': sqrtspol*1,
        'Te': sqrtspol*2,
        'Ti': sqrtspol*3,
        'vrot': sqrtspol*4,
    }
    sqrtstor = trial_sqrtspol2sqrtstor(sqrtspol)
    return profiles, sqrtspol, sqrtstor

def trial_profiles_spol(spol):
    sqrtspol = np.sqrt(spol)
    return trial_profiles_sqrtspol(sqrtspol)

def trial_profiles_sqrtstor(sqrtstor):
    sqrtspol = trial_sqrtstor2sqrtspol(sqrtstor)
    return trial_profiles_sqrtspol(sqrtspol)

def trial_sqrtstor2sqrtspol(sqrtstor):
    sqrtspol = sqrtstor**(3/1)
    return sqrtspol

def trial_sqrtspol2sqrtstor(sqrtspol):
    sqrtstor = sqrtspol**(1/3)
    return sqrtstor

def test_interpolation_visual_check():
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(3, 3, figsize=(12, 4))
    write_trial_profiles(trial_profiles_sqrtspol, test_output_dir)
    trial_profiles, sqrtspol, sqrtstor = load_profiles_and_interp(test_src, 'sqrtspol')
    plot_profiles(ax[:,0], trial_profiles, sqrtspol)
    write_trial_profiles(trial_profiles_spol, test_output_dir)
    trial_profiles, sqrtspol, sqrtstor = load_profiles_and_interp(test_src, 'spol')
    plot_profiles(ax[:,1], trial_profiles, sqrtspol)
    write_trial_profiles(trial_profiles_sqrtstor, test_output_dir)
    trial_profiles, sqrtspol, sqrtstor = load_profiles_and_interp(test_src, 'sqrtstor')
    plot_profiles(ax[:,2], trial_profiles, sqrtspol)
    plt.subplots_adjust(hspace=0.8, wspace=0.3)
    ax[0,0].set_title('equidist sqrtspol grid')
    ax[0,1].set_title('equidist spol grid')
    ax[0,2].set_title('equidist sqrtstor grid')
    plt.show()

def plot_profiles(ax, profiles, sqrtspol):
    xs = [sqrtspol, sqrtspol**2, trial_sqrtspol2sqrtstor(sqrtspol)]
    xnames = ['sqrtspol', 'spol', 'sqrtstor']
    for i,x in enumerate(xs):
        for profile in profiles:
            ax[i].plot(x, profiles[profile],'o', label=profile)
        ax[i].set_xlabel(xnames[i] + ' [1]')
        ax[i].set_ylabel('profiles [1]')
    

if __name__ == '__main__':
    test_write_and_read_compatiblity()
    test_equidistant_grid()
    test_interpolation()
    test_interpolation_visual_check()