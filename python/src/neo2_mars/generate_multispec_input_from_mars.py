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
    config['Z'], config['m'] = get_species_cgs_from_mars(mars_dir)
    profiles_src = {
    'sqrtspol': {'filename': os.path.join(working_dir, 'sqrtstor.dat'), 'column': 0},
    'sqrtstor': {'filename': os.path.join(working_dir, 'sqrtstor.dat'), 'column': 1},
    'n': {'filename': os.path.join(working_dir, 'n.dat'), 'column': [1,2]},
    'T': {'filename': os.path.join(working_dir, 'T.dat'), 'column': [1,2]},
    'vrot': {'filename': os.path.join(working_dir, 'vrot.dat'), 'column': 1}
    }
    profiles_interp_config = {'n_s': number_of_surfaces, 'min_s': bounds[0], 'max_s': bounds[1],
                              'grid': grid}
    generate_multispec_input(config, profiles_src, profiles_interp_config)

def get_species_cgs_from_mars(mars_dir):
    PROTON_MASS_CGS = 1.67262192369e-24

    try:
        import f90nml
    except ModuleNotFoundError:
        f90nml = None

    run_in = os.path.join(mars_dir, 'RUN.IN')
    if f90nml is not None:
        run_config = f90nml.read(run_in)
        Zs = run_config['KINETIC']['ESPECIES_Z']
        mass_numbers = run_config['KINETIC']['ESPECIES_M']
    else:
        Zs, mass_numbers = _parse_species_from_run_in(run_in)

    for Z in Zs:
        if Z < 0:
            Ze = Z
        if Z > 0:
            Zi = Z

    me, mi = mass_numbers[Zs.index(Ze)] * PROTON_MASS_CGS, mass_numbers[Zs.index(Zi)] * PROTON_MASS_CGS
    Z = [Ze, Zi]
    m = [me, mi]
    return Z, m


def _parse_species_from_run_in(path):
    with open(path, 'r', encoding='utf-8') as handle:
        content = handle.read()

    def _values(name):
        for line in content.splitlines():
            if name in line:
                _, values = line.split('=', 1)
                return [float(value.strip().strip(',')) for value in values.split(',') if value.strip()]
        raise ValueError(f'Could not find {name} in {path}')

    z_values = [int(value) for value in _values('ESPECIES_Z')]
    m_values = _values('ESPECIES_M')
    return z_values, m_values
