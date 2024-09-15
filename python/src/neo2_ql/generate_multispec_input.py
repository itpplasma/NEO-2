import h5py
import numpy as np

from .load_profile_data import load_cgs_profiles_and_interp

def generate_multispec_input(config: dict, profiles_src: dict, profiles_interp_config: dict={}):
    multispec = {}
    nspecies = len(config['Z'])
    multispec['/num_species'] = nspecies
    multispec['/species_tag'] = np.array(range(nspecies), dtype=np.int32) + 1
    multispec['/species_tag_Vphi'] = config['species_tag_of_vrot']
    multispec['/isw_Vphi_loc'] = 0

    profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(profiles_src, profiles_interp_config)

    multispec['/num_radial_pts'] = len(sqrtspol)
    multispec['/species_def'] = get_species_def_array(config, multispec['/num_radial_pts'])
    multispec['/rel_stages'] = np.full(multispec['/num_radial_pts'], nspecies)
    multispec['/rho_pol'] = sqrtspol
    stor = sqrtstor**2
    multispec['/boozer_s'] = stor
    multispec['/T_prof'] = profiles['T']
    multispec['/dT_ov_ds_prof'] = derivative(stor, profiles['T'])
    multispec['/n_prof'] = profiles['n']
    multispec['/dn_ov_ds_prof'] = derivative(stor, profiles['n'])
    multispec['/Vphi'] = profiles['vrot']
    ELEMENTARY_CHARGE_CGS = 1.60217662e-19 * 2.99792458e8 * 10
    charge = np.array(config['Z'])*ELEMENTARY_CHARGE_CGS
    charge = charge[:, np.newaxis]
    multispec['/kappa_prof'] = get_kappa(profiles['n'], profiles['T'], charge)
    
    write_multispec_to_hdf5(config['hdf5_filename'], multispec)

def get_species_def_array(config, num_radial_pts):
    species_def = np.vstack([config['Z'], config['m']]).T
    #species_def = [[config['Ze'], config['me']], [config['Zi'], config['mi']]]
    return np.array([species_def]*num_radial_pts).T

def derivative(x,y):
    if y.ndim > 1:
        return matrix_derivative(x,y)
    elif y.ndim == 1:
        return vector_derivative(x,y)
    else:
        raise ValueError('y must be a vector or 2D matrix')

def matrix_derivative(x,y):
    dy = np.zeros_like(y)
    for i in range(y.shape[0]):
        dy[i] = vector_derivative(x, y[i])
    return dy

def vector_derivative(x,y):
    left_border = (-1/2*y[2] + 2*y[1] -3/2*y[0]) / (x[1] - x[0])
    middle = (y[2:] - y[:-2]) / (x[2:] - x[:-2])
    right_border = (3/2*y[-1] - 2*y[-2] + 1/2*y[-3]) / (x[-1] - x[-2])
    return np.concatenate([[left_border], middle, [right_border]])

def get_kappa(n_cgs, T_cgs, charge_cgs):
    log_Lambda = get_coulomb_logarithm(n_cgs, T_cgs)
    mean_free_path = (3 / (4 * np.sqrt(np.pi)) * (T_cgs/charge_cgs) ** 2 / 
                                    (n_cgs * (charge_cgs ** 2) * log_Lambda))
    kappa = 2 / mean_free_path
    return kappa

def get_coulomb_logarithm(n_cgs, T_cgs):
    # Coulomb logarithm (set as species-independent - see E A Belli and J Candy PPCF 54 015015 (2012))
    log_Lambda = 52.43 - 1.15 * np.log10(n_cgs) + 2.3 * np.log10(T_cgs)
    return log_Lambda

def write_multispec_to_hdf5(hdf5_filename, multispec):
    with h5py.File(hdf5_filename, 'w') as f:
        f.create_dataset('/num_radial_pts', data=[multispec['/num_radial_pts']], dtype='int32')
        f.create_dataset('/num_species', data=[multispec['/num_species']], dtype='int32')
        f.create_dataset('/species_tag', data=multispec['/species_tag'], dtype='int32')
        f.create_dataset('/species_def', data=multispec['/species_def'], dtype='float64')
        f['species_def'].attrs['unit'] = 'e; g'
        f.create_dataset('/boozer_s', data=multispec['/boozer_s'], dtype='float64')
        f['boozer_s'].attrs['unit'] = '1'
        f.create_dataset('/rho_pol', data=multispec['/rho_pol'], dtype='float64')
        f['rho_pol'].attrs['unit'] = '1'
        f.create_dataset('/Vphi', data=multispec['/Vphi'], dtype='float64')
        f['Vphi'].attrs['unit'] = 'rad / s'
        f.create_dataset('/species_tag_Vphi', data=[multispec['/species_tag_Vphi']], dtype='int32')
        f.create_dataset('/isw_Vphi_loc', data=[multispec['/isw_Vphi_loc']], dtype='int32')
        f.create_dataset('/rel_stages', data=multispec['/rel_stages'], dtype='int32')
        f.create_dataset('/T_prof', data=multispec['/T_prof'], dtype='float64')
        f['T_prof'].attrs['unit'] = 'erg'
        f.create_dataset('/dT_ov_ds_prof', data=multispec['/dT_ov_ds_prof'], dtype='float64')
        f.create_dataset('/n_prof', data=multispec['/n_prof'], dtype='float64')
        f['n_prof'].attrs['unit'] = '1 / cm^3'
        f.create_dataset('/dn_ov_ds_prof', data=multispec['/dn_ov_ds_prof'], dtype='float64')
        f.create_dataset('/kappa_prof', data=multispec['/kappa_prof'], dtype='float64')