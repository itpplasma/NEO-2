# %% Standard imports
import numpy as np

# Homebrew imports
from neo2_ql import generate_multispec_input
from .mars_profiles_to_neo2_ql_profiles import write_neo2_input_profiles_from_mars

def generate_multispec_input_from_mars(mars_dir: str, output_file: str='multi_spec.in',
                                       number_of_surfaces: int=10, bounds: list=[0, 1]):
    write_neo2_input_profiles_from_mars(mars_dir, '/tmp/')

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
    species_definition = species_definition.transpose(2, 1, 0)
    isw_Vphi_loc = 0
    species_tag_Vphi = 2
    input_unit_type = 2
    switch_grid = 2
    generate_multispec_input(hdf5FileName, path_to_input_files, data_source, species_definition, isw_Vphi_loc,
                             species_tag_Vphi, input_unit_type, bounds, switch_grid)

def get_species_cgs_form_mars(mars_dir):
    import f90nml
    PROTON_MASS_CGS = 1.6726219e-24
    run_config = f90nml.read(os.path.join(mars_dir, 'RUN.IN'))
    Z = run_config['KINETIC']['ESPECIES_Z']
    Ze, Zi = Z[Z < 0], Z[Z > 0]
    mass_numbers = run_config['KINETIC']['ESPECIES_M']
    me, mi = (mass_numbers[Z.index(Ze)], mass_numbers[Z.index(Zi)]) * PROTON_MASS_CGS
    return Ze, Zi, me, mi