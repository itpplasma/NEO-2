# %%Standard Python modules
import numpy as np
from numpy.testing import assert_allclose
import matplotlib.pyplot as plt

# Custom modules
from libneo import FluxConverter, read_eqdsk

#Module to test
from neo2_ql import interpolate_profiles_to_same_grid

sqrtspol = np.linspace(0, 1, 10)

def test_interpolate_profiles_to_same_grid():
    test_no_change_to_profiles_on_same_grid()
    test_profiles_on_different_grids()

def test_no_change_to_profiles_on_same_grid():
    profiles = trial_profiles_on_same_grid(sqrtspol)
    interp_profiles = interpolate_profiles_to_same_grid(profiles)
    assert is_same_grid_for_all(interp_profiles)
    assert is_same_for_all_entries(profiles, interp_profiles)

def test_profiles_on_different_grids():
    profiles = trial_profiles_on_different_grids(sqrtspol)
    interp_profiles = interpolate_profiles_to_same_grid(profiles)
    assert is_same_grid_for_all(interp_profiles)
    control_profiles = trial_profiles_on_same_grid(sqrtspol)
    assert is_same_for_all_entries(control_profiles, interp_profiles)

def is_same_grid_for_all(profiles):
    for key, profile in profiles.items():
        if not np.allclose(profile[:,0], profiles['sqrtstor'][:,0]):
            return False
    return True

def is_same_for_all_entries(profiles1, profiles2):
    for key, profile in profiles1.items():
        if not np.allclose(profiles1[key][:,1:], profiles2[key][:,1:], atol=1e-4):
            print(f"Profiles differ for {key}")
            print(profiles1[key][:,1:])
            print(profiles2[key][:,1:])
            return False
    return True

def trial_profiles_on_same_grid(sqrtspol):
    profiles = {}
    profiles['sqrtstor'] = np.array([sqrtspol, sqrtspol**2]).T
    profiles['n'] = np.array([sqrtspol, 2*sqrtspol, 3*sqrtspol]).T
    profiles['T'] = np.array([sqrtspol, 4*sqrtspol, 5*sqrtspol]).T
    profiles['vrot'] = np.array([sqrtspol, 6*sqrtspol]).T
    return profiles

def trial_profiles_on_different_grids(sqrtspol):
    profiles = {}
    profiles['sqrtstor'] = np.array([sqrtspol, sqrtspol**2]).T
    profiles['n'] = np.array([sqrtspol, 2*sqrtspol, 3*sqrtspol]).T
    shift_sqrtspol = sqrtspol[:-1] + np.diff(sqrtspol)/2
    profiles['T'] = np.array([shift_sqrtspol, 4*shift_sqrtspol, 5*shift_sqrtspol]).T
    profiles['vrot'] = np.array([shift_sqrtspol, 6*shift_sqrtspol]).T
    return profiles

if __name__ == "__main__":
    test_interpolate_profiles_to_same_grid()
    print("All tests passed!")