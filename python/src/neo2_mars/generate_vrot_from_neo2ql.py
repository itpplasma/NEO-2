import numpy as np
import h5py

def generate_vrot_for_mars(neo2ql_input_file):
    sqrtspol, ion_vrot = get_vrot_from_neo2ql(neo2ql_input_file)
    write_vrot_to_mars_input(ion_vrot, sqrtspol)

def get_vrot_from_neo2ql(neo2ql_input_file):
    neo2ql = h5py.File(neo2ql_input_file, "r")
    sqrtspol = np.array(neo2ql["rho_pol"])
    ion_vrot = np.array(neo2ql["Vphi"])
    return sqrtspol, ion_vrot

def write_vrot_to_mars_input(vrot, sqrtspol):
    vrot_file = "PROFROT.IN"
    number_of_surfaces = len(vrot)
    type_of_radial_variable = 1 # 1 means sqrtspol for MARS
    header = f"{number_of_surfaces} {type_of_radial_variable}"
    np.savetxt(vrot_file, np.array([sqrtspol, vrot]).T, header=header, comments="")