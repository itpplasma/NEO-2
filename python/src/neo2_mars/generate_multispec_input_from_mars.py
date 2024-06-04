# %% Standard imports
import numpy as np
import os

# Homebrew imports
from neo2_ql import generate_multispec_input
from .mars_profiles_to_neo2_ql_profiles import write_neo2_input_profiles_from_mars

def generate_multispec_input_from_mars(mars_dir: str, output_file: str='multi_spec.in',
                                       number_of_surfaces: int=10, bounds: list=[0, 1], grid: str='spol'):
    working_dir = '/tmp/'
    write_neo2_input_profiles_from_mars(mars_dir, working_dir)
    config = {
            'hdf5_filename': output_file, 
            'species_tag_of_vrot': 2, 
             }
    config['Ze'], config['Zi'], config['me'], config['mi'] = get_species_cgs_from_mars(mars_dir)
    profiles_src = {
    'sqrtspol': {'filename': os.path.join(working_dir, 'sqrtstor.dat'), 'column': 0},
    'sqrtstor': {'filename': os.path.join(working_dir, 'sqrtstor.dat'), 'column': 1},
    'ne': {'filename': os.path.join(working_dir, 'ne.dat'), 'column': 1},
    'ni': {'filename': os.path.join(working_dir, 'ne.dat'), 'column': 1}, 
    'Te': {'filename': os.path.join(working_dir, 'Te.dat'), 'column': 1},
    'Ti': {'filename': os.path.join(working_dir, 'Ti.dat'), 'column': 1},
    'vrot': {'filename': os.path.join(working_dir, 'vrot.dat'), 'column': 1}
    }
    profiles_interp_config = {'n_s': number_of_surfaces, 'min_s': bounds[0], 'max_s': bounds[1],
                              'grid': grid}
    generate_multispec_input(config, profiles_src, profiles_interp_config)

def get_species_cgs_from_mars(mars_dir):
    import f90nml
    import os
    PROTON_MASS_CGS = 1.67262192369e-24
    run_config = f90nml.read(os.path.join(mars_dir, 'RUN.IN'))
    Zs = run_config['KINETIC']['ESPECIES_Z']
    for Z in Zs:
        if Z < 0:
            Ze = Z
        if Z > 0:
            Zi = Z
    mass_numbers = run_config['KINETIC']['ESPECIES_M']
    me, mi = mass_numbers[Zs.index(Ze)] * PROTON_MASS_CGS, mass_numbers[Zs.index(Zi)] * PROTON_MASS_CGS
    return Ze, Zi, me, mi