#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# %% Standard libraries
import numpy as np
import os

# Custom libraries
from neo2_ql import interpolate_profiles_to_same_grid
from neo2_ql import write_profiles_to_dat_files

########################################################################################

def write_neo2_input_profiles_from_mars(mars_folder: str, output_dir: str):
    profiles_mars = get_profiles_mars(mars_folder)
    profiles_mars = interpolate_profiles_to_same_grid(profiles_mars)
    write_profiles_to_dat_files(profiles_mars, output_dir)

########################################################################################

def get_profiles_mars(mars_dir: str) -> dict:
    profiles_mars = {}
    profiles_mars['n'] = np.loadtxt(os.path.join(mars_dir,'PROFDEN.IN'), skiprows=1)
    ni = profiles_mars['n']
    profiles_mars['n'] = np.hstack((profiles_mars['n'], ni[:,1:]))
    profiles_mars['T'] = np.loadtxt(os.path.join(mars_dir,'PROFTE.IN'), skiprows=1)
    Ti = np.loadtxt(os.path.join(mars_dir,'PROFTI.IN'), skiprows=1)
    profiles_mars['T'] = np.hstack((profiles_mars['T'], Ti[:,1:]))
    profiles_mars['vrot'] = np.loadtxt(os.path.join(mars_dir,'PROFROT.IN'), skiprows=1)
    profiles_mars['sqrtstor'] = get_sqrtstor_profile(mars_dir)
    return profiles_mars

def get_sqrtstor_profile(mars_dir: str) -> np.ndarray:
    sqrtspol = get_mars_sqrtspol(mars_dir)
    sqrtstor = mars_sqrtspol2sqrtstor(mars_dir,sqrtspol)
    return np.array([sqrtspol, sqrtstor]).T

def get_mars_sqrtspol(mars_dir: str) -> np.ndarray:
    path_to_file = os.path.join(mars_dir, 'PROFDEN.IN')
    sqrtspol = np.loadtxt(path_to_file, skiprows=1)[:, 0]
    return sqrtspol

def mars_sqrtspol2sqrtstor(mars_dir,sqrtspol):
    return np.sqrt(mars_sqrtspol2stor(mars_dir,sqrtspol))

def mars_sqrtspol2stor(mars_dir,sqrtspol):
    from libneo import FluxConverter
    q_over_equidist_spol = get_mars_q_over_equidist_spol(mars_dir)
    stor = FluxConverter(q_over_equidist_spol).spol2stor(sqrtspol ** 2)
    return stor

def get_mars_q_over_equidist_spol(mars_dir: str) -> np.ndarray:
    q, sqrtspol = get_mars_q(mars_dir)
    q_over_equidist_spol = interp_to_equidist_grid(q, sqrtspol**2)
    return q_over_equidist_spol

def get_mars_q(mars_dir: str) -> dict:
    from omfit_classes.omfit_mars import OMFITmars
    data = OMFITmars(mars_dir)
    q = data['PROFEQ']['q_prof'].values
    sqrtspol = data['PROFEQ']['q_prof'].coords['s_eq'].values
    return q, sqrtspol

def interp_to_equidist_grid(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    equidist_x = np.linspace(np.min(x), np.max(x), x.shape[0])
    equidist_y = np.interp(equidist_x, x, y)
    return equidist_y