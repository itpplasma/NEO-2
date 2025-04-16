import numpy as np
import h5py

def generate_omega_e_for_mars(neo2ql_output_file, neo2ql_input_file):
    species_omega_e, stor = get_omega_e_from_neo2ql(neo2ql_output_file)
    sqrtspol = neo2ql_stor2sqrtspol(neo2ql_input_file, stor)
    omega_e = species_omega_e[:,0].T
    write_omega_e_to_mars_input(-omega_e, sqrtspol) # MARS has opposite phi direction

def get_omega_e_from_neo2ql(neo2ql_output_file):
    neo2ql = h5py.File(neo2ql_output_file, "r")
    species_mach_over_major_radius = np.array(neo2ql["MtOvR"])
    species_temperature = np.array(neo2ql["T_spec"])
    species_mass = np.array(neo2ql["m_spec"])
    species_omega_e = species_mach_over_major_radius * np.sqrt(2*species_temperature/species_mass)
    return species_omega_e, np.array(neo2ql["boozer_s"])

def neo2ql_stor2sqrtspol(neo2ql_input_file, stor):
    neo2ql = h5py.File(neo2ql_input_file, "r")
    stor_i = np.array(neo2ql["boozer_s"]).ravel()
    sqrtspol_i = np.array(neo2ql["rho_pol"]).ravel()
    sqrtspol = np.interp(stor, stor_i, sqrtspol_i)
    return sqrtspol

def write_omega_e_to_mars_input(omega_e, sqrtspol):
    omega_e_file = "PROFWE.IN"
    number_of_surfaces = len(omega_e)
    type_of_radial_variable = 1 # 1 means sqrtspol for MARS
    header = f"{number_of_surfaces} {type_of_radial_variable}"
    np.savetxt(omega_e_file, np.array([sqrtspol, omega_e]).T, header=header, comments="")
