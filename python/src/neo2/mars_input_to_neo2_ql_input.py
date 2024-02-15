# %% Standard libraries
import numpy as np
import os
from omfit_classes.omfit_mars import OMFITmars

# Homebrew libraries
from libneo import FluxConverter

########################################################################################

def write_input_for_generate_neo2_profile_from_mars(mars_folder: str, output_dir: str):
    input_profiles_mars = get_input_profiles_mars(mars_folder)
    write_dat_files_from_profiles(input_profiles_mars, output_dir)

########################################################################################

def get_input_profiles_mars(path_to_mars_folder: str) -> dict:
    filenames = ['PROFDEN', 'PROFTE', 'PROFTI', 'PROFROT']
    input_profiles_mars = {}
    for filename in filenames:
        path_to_file = os.path.join(path_to_mars_folder, filename + '.IN')
        input_profiles_mars[filename] = np.loadtxt(path_to_file, skiprows=1)
    input_profiles_mars['PROFSQRTSTOR'] = get_sqrt_stor_profile(path_to_mars_folder)
    return input_profiles_mars

def get_sqrt_stor_profile(path_to_mars_folder: str) -> np.ndarray:
    sqrt_spol = get_sqrt_spol_of_input_profiles_mars(path_to_mars_folder)
    sqrt_stor = get_sqrt_stor_of_input_profiles_mars(path_to_mars_folder)
    sqrt_stor_profile = np.array([sqrt_spol, sqrt_stor]).T
    return sqrt_stor_profile

########################################################################################

def get_sqrt_stor_of_input_profiles_mars(path_to_mars_folder: str) -> np.ndarray:
    sqrt_spol_mars = get_sqrt_spol_of_input_profiles_mars(path_to_mars_folder)
    q_prof = get_q_prof_mars(path_to_mars_folder)
    sqrt_stor_mars = convert_sqrt_spol_to_sqrt_stor(q_prof['sqrt_spol'], q_prof['values'],
                                                    sqrt_spol_mars)
    return sqrt_stor_mars

def get_sqrt_spol_of_input_profiles_mars (path_to_mars_folder: str) -> np.ndarray:
    path_to_file = os.path.join(path_to_mars_folder, 'PROFDEN.IN')
    sqrt_spol = np.loadtxt(path_to_file, skiprows=1)[:, 0]
    return sqrt_spol

def get_q_prof_mars(path_to_mars_folder: str) -> dict:
    data = OMFITmars(path_to_mars_folder)
    data.load()
    q_prof = {}
    q_prof['values'] = data['PROFEQ']['q_prof'].values
    q_prof['sqrt_spol'] = data['PROFEQ']['q_prof'].coords['s_eq'].values
    return q_prof

def convert_sqrt_spol_to_sqrt_stor(q_prof_sqrt_spol, q_prof_values, sqrt_spol_to_evaluate):
    _, q_prof_over_equidist_spol = get_q_prof_over_equidist_spol(q_prof_sqrt_spol, q_prof_values)
    stor = convert_sqrt_spol_to_stor(q_prof_over_equidist_spol, sqrt_spol_to_evaluate)
    return np.sqrt(stor)

def get_q_prof_over_equidist_spol(q_prof_sqrt_spol, q_prof_values):
    spol_prof = q_prof_sqrt_spol**2
    equidist_spol = np.linspace(np.min(spol_prof), np.max(spol_prof), spol_prof.shape[0])
    q_prof_over_equidist_spol = np.interp(equidist_spol, spol_prof, q_prof_values)
    return equidist_spol, q_prof_over_equidist_spol

def convert_sqrt_spol_to_stor(q_prof_over_equidist_spol, sqrt_spol_to_evaluate):
    converter = FluxConverter(q_prof_over_equidist_spol)
    spol_to_evaluate = sqrt_spol_to_evaluate ** 2
    stor = converter.spol2stor(spol_to_evaluate)
    return stor
########################################################################################

def write_dat_files_from_profiles(input_profiles_mars: dict, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    for profile_name, profile_data in input_profiles_mars.items():
        filename = f"{output_dir}/{profile_name}.dat"
        np.savetxt(filename, profile_data)

########################################################################################

# def add_variable_to_Dataset(Dataset,variable_data,variable_coordinate:str, variable_name:str):
#     existing_variable_name = list(Dataset.variables)[0]
#     Dataset = Dataset.assign(**{variable_name: ((variable_coordinate), variable_data)})
#     return Dataset

########################################################################################

# def write_dat_files_from_dataset(Dataset, variables_to_write: list, output_dir: str):
#     os.makedirs(output_dir, exist_ok=True)
#     for variable_name in  variables_to_write:
#         variable = Dataset[variable_name]
#         filename = f"{output_dir}/{variable_name}.dat"
#         write_dat_file_from_dataset_variable(variable, filename)

# def write_dat_file_from_dataset_variable(variable, filename):

#     dataframe = variable.to_dataframe()
#     for coord in variable.coords:
#         dataframe.insert(0, coord, variable.coords[coord])
#     dataframe.to_csv(filename, sep=' ', index=False, header=False)

########################################################################################

def get_mars_q_prof_over_stor(path_to_mars_folder: str):
    input_profiles_mars = get_input_profiles_mars(path_to_mars_folder)
    stor = convert_sqrt_spol_to_stor(input_profiles_mars['rho_pol'],input_profiles_mars['q_prof_values'])
    return stor, input_profiles_mars['q_prof_values']

########################################################################################

if __name__ == "__main__":

    path_to_mars_folder = "/proj/plasma/DATA/DEMO/MARS/MARSQ_INPUTS_KNTV21_NEO2profs_RUN/"
    output_dir = "/temp/grassl_g/TEST_NTV_DEMO/input_files_for_generate_neo2_profile"
    #dataset_mars = get_input_profiles_mars(path_to_mars_folder)
    write_input_for_generate_neo2_profile_from_mars(path_to_mars_folder, output_dir)
