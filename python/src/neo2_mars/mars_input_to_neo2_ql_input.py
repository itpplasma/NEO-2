#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# %% Standard libraries
import numpy as np
import os
from omfit_classes.omfit_mars import OMFITmars

# Homebrew libraries
from libneo import FluxConverter

########################################################################################

def write_input_for_generate_neo2_profile_from_mars(mars_folder: str, output_dir: str):
    input_profiles_mars = get_input_profiles_mars(mars_folder)
    write_profiles_to_dat_files(input_profiles_mars, output_dir)

########################################################################################

def get_input_profiles_mars(path_to_mars_folder: str) -> dict:
    input_profiles_mars = get_input_profiles_mars_except_sqrt_stor(path_to_mars_folder)
    input_profiles_mars['PROFSQRTSTOR'] = get_sqrt_stor_profile(path_to_mars_folder)
    return input_profiles_mars

def get_input_profiles_mars_except_sqrt_stor(path_to_mars_folder: str) -> dict:
    filenames = ['PROFDEN', 'PROFTE', 'PROFTI', 'PROFROT']
    input_profiles_mars = {}
    for filename in filenames:
        path_to_file = os.path.join(path_to_mars_folder, filename + '.IN')
        input_profiles_mars[filename] = np.loadtxt(path_to_file, skiprows=1)
    return input_profiles_mars

def get_sqrt_stor_profile(path_to_mars_folder: str) -> np.ndarray:
    sqrt_spol = get_sqrt_spol_of_input_profiles_mars(path_to_mars_folder)
    sqrt_stor = get_sqrt_stor_of_input_profiles_mars(path_to_mars_folder)
    sqrt_stor_profile = np.array([sqrt_spol, sqrt_stor]).T
    return sqrt_stor_profile

def get_sqrt_stor_of_input_profiles_mars(path_to_mars_folder: str) -> np.ndarray:
    sqrt_spol_mars = get_sqrt_spol_of_input_profiles_mars(path_to_mars_folder)
    q_over_equidist_spol_mars = get_q_over_equidist_spol_mars(path_to_mars_folder)
    sqrt_stor_mars = convert_sqrt_spol_to_sqrt_stor(q_over_equidist_spol_mars,sqrt_spol_mars)
    return sqrt_stor_mars

def get_sqrt_spol_of_input_profiles_mars (path_to_mars_folder: str) -> np.ndarray:
    path_to_file = os.path.join(path_to_mars_folder, 'PROFDEN.IN')
    sqrt_spol = np.loadtxt(path_to_file, skiprows=1)[:, 0]
    return sqrt_spol

def mars_sqrtspol2stor(path_to_mars_folder,sqrt_spol):
    q_over_equidist_spol = get_q_over_equidist_spol_mars(path_to_mars_folder)
    converter = FluxConverter(q_over_equidist_spol)
    stor = converter.spol2stor(sqrt_spol ** 2)
    return stor

def get_q_over_equidist_spol_mars(path_to_mars_folder: str) -> np.ndarray:
    q_prof = get_q_prof_mars(path_to_mars_folder)
    spol = q_prof['sqrt_spol']**2
    equidist_spol = np.linspace(np.min(spol), np.max(spol), spol.shape[0])
    q_over_equidist_spol = np.interp(equidist_spol, spol, q_prof['values'])
    return q_over_equidist_spol

def get_q_prof_mars(path_to_mars_folder: str) -> dict:
    data = OMFITmars(path_to_mars_folder)
    data.load()
    q_prof = {}
    q_prof['values'] = data['PROFEQ']['q_prof'].values
    q_prof['sqrt_spol'] = data['PROFEQ']['q_prof'].coords['s_eq'].values
    return q_prof

def convert_sqrt_spol_to_sqrt_stor(q_over_equidist_spol, sqrt_spol_to_evaluate):
    return np.sqrt(convert_sqrt_spol_to_stor(q_over_equidist_spol, sqrt_spol_to_evaluate))

def convert_sqrt_spol_to_stor(q_prof_over_equidist_spol, sqrt_spol_to_evaluate):
    converter = FluxConverter(q_prof_over_equidist_spol)
    spol_to_evaluate = sqrt_spol_to_evaluate ** 2
    stor = converter.spol2stor(spol_to_evaluate)
    return stor

########################################################################################

def write_profiles_to_dat_files(input_profiles_mars: dict, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    for profile_name, profile_data in input_profiles_mars.items():
        filename = f"{output_dir}/{profile_name}.dat"
        np.savetxt(filename, profile_data)

########################################################################################