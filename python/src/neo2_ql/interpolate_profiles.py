#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

def interpolate_profiles_to_same_grid(profiles: dict) -> dict:
    sqrtspol = profiles['sqrtstor'][:, 0]
    interp_profiles = {}
    for kind in profiles.keys():
        interp_profiles[kind] = np.zeros((sqrtspol.shape[0], profiles[kind].shape[1]))
        if profiles[kind].shape[1] > 2:
            interp_profiles[kind][:,1:] = matrix_interpolate(sqrtspol, profiles[kind][:,0], 
                                                             profiles[kind][:,1:])
        else:
            interp_profiles[kind][:,1] = vector_interpolate(sqrtspol, profiles[kind][:,0], 
                                                            profiles[kind][:,1])
        interp_profiles[kind][:, 0] = sqrtspol.copy()
    return interp_profiles

def matrix_interpolate(x_new: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    y_new = np.zeros((x_new.shape[0], y.shape[1]))
    for i in range(y.shape[1]):
        y_new[:, i] = vector_interpolate(x_new, x, y[:, i])
    return y_new

def vector_interpolate(x_new: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    from scipy.interpolate import interp1d
    f = interp1d(x, y, fill_value='extrapolate')
    return f(x_new)