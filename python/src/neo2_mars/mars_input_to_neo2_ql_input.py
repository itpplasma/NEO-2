#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# %% Standard libraries
import numpy as np
import os

########################################################################################

def write_neo2_input_profile_from_mars(mars_folder: str, output_dir: str):
    profiles_mars = get_profiles_mars(mars_folder)
    write_profiles_to_dat_files(profiles_mars, output_dir)

########################################################################################

def get_profiles_mars(mars_dir: str) -> dict:
    profiles_mars = {}
    profiles_mars['ne'] = np.loadtxt(os.path.join(mars_dir,'PROFDEN.IN'), skiprows=1)
    profiles_mars['Te'] = np.loadtxt(os.path.join(mars_dir,'PROFTE.IN'), skiprows=1)
    profiles_mars['Ti'] = np.loadtxt(os.path.join(mars_dir,'PROFTI.IN'), skiprows=1)
    profiles_mars['vrot'] = np.loadtxt(os.path.join(mars_dir,'PROFROT.IN'), skiprows=1)
    profiles_mars['sqrtstor'] = get_sqrtstor_profile(mars_dir)
    return profiles_mars

def get_sqrtstor_profile(mars_dir: str) -> np.ndarray:
    sqrtspol = get_mars_sqrtspol(mars_dir)
    sqrtstor = mars_sqrtspol2sqrtstor(mars_dir,sqrtspol)
    return np.array([sqrtspol, sqrtstor]).T

def get_mars_sqrtspol (mars_dir: str) -> np.ndarray:
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

########################################################################################

def write_profiles_to_dat_files(profiles_mars: dict, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    for profile_name, profile_data in profiles_mars.items():
        filename = f"{output_dir}/{profile_name}.dat"
        np.savetxt(filename, profile_data)

########################################################################################