import numpy as np
import h5py

def generate_nu_for_mars(neo2ql_input_file):
    sqrtspol, nu = get_nu_from_neo2ql(neo2ql_input_file)
    write_nu_to_mars_input(nu, sqrtspol)

def get_nu_from_neo2ql(neo2ql_input_file):
    neo2ql = h5py.File(neo2ql_input_file, "r")
    sqrtspol = np.array(neo2ql["rho_pol"])
    species_kappa = np.array(neo2ql["kappa_prof"])
    species_temperature = np.array(neo2ql["T_prof"])
    species_mass = np.array(neo2ql["species_def"])[1]
    species_thermal_velocity = np.sqrt(2*species_temperature/species_mass)
    species_nu = species_kappa/2 * species_thermal_velocity
    return sqrtspol, species_nu

def write_nu_to_mars_input(nu, sqrtspol):
    nue_file = "PROFNUE.IN"
    nui_file = "PROFNUI.IN"
    number_of_surfaces = len(nu[0])
    type_of_radial_variable = 1 # 1 means sqrtspol for MARS
    header = f"{number_of_surfaces} {type_of_radial_variable}"
    np.savetxt(nue_file, np.array([sqrtspol, nu[0]]).T, header=header, comments="")
    np.savetxt(nui_file, np.array([sqrtspol, nu[1]]).T, header=header, comments="")