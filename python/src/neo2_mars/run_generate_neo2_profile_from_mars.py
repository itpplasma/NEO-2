# %% Standard imports
import numpy as np

# Homebrew imports
from .generate_neo2_profile import generate_neo2_profile

def run_generate_neo2_profile(path_to_input_files, number_of_surfaces, bounds):
    hdf5FileName = 'multi_spec_demo.in'
    path_to_input_files = './input_files_for_generate_neo2_profile/'
    
    data_source = {}
    data_source['rhopoloidal'] = {'filename': 'PROFDEN.dat', 'column': 0}
    data_source['electron_density'] = {'filename': 'PROFDEN.dat', 'column': 1}
    data_source['rhotoroidal'] = {'filename': 'PROFSQRTSTOR.dat', 'column': 1}
    data_source['electron_temperature'] = {'filename': 'PROFTE.dat', 'column': 1}
    data_source['ion_temperature'] = {'filename': 'PROFTI.dat', 'column': 1}
    data_source['rotation_velocity'] = {'filename': 'PROFROT.dat', 'column': 1}
    
    Zi = 1  # ion charge number (deuterium)
    mi = 3.3436e-24  # ion mass (deuterium)
    Ze = -1  # electron charge number
    me = 9.1094e-28  # electron mass
    species_definition = np.array([[[Ze, me], [Zi, mi]]] * number_of_surfaces)
    isw_Vphi_loc = 0
    species_tag_Vphi = 2
    input_unit_type = 2
    switch_grid = 2
    generate_neo2_profile(hdf5FileName, path_to_input_files, data_source, species_definition, isw_Vphi_loc,
                          species_tag_Vphi, input_unit_type, bounds, switch_grid)