import numpy as np
import h5py

def generate_omegate_from_neo2ql(neo2ql_output_file):
    omega_te = get_omegate_from_neo2ql(neo2ql_output_file)
    pass

def get_omegate_from_neo2ql(neo2ql_output_file):
    neo2ql = h5py.File(neo2ql_output_file, "r")
    species_mach_over_major_radius = np.array(neo2ql["MtOvR"])
    species_temperature = np.array(neo2ql["T_spec"])
    species_mass = np.array(neo2ql["m_spec"])
    species_omegate = species_mach_over_major_radius * np.sqrt(2*species_temperature/species_mass)
    return species_omegate, np.array(neo2ql["boozer_s"])
    

if __name__ == "__main__":
    pass