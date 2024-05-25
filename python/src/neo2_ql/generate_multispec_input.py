import h5py
import numpy as np

from .load_profile_data import load_cgs_profiles_and_interp

def generate_multispec_input(config: dict, profiles_src: dict, profiles_interp_config: dict={}):
    multispec = {}
    multispec['/num_species'] = 2
    multispec['/species_tag'] = np.array([1, 2], dtype=np.int32)
    multispec['/species_tag_Vphi'] = config['species_tag_Vphi']
    multispec['/isw_Vphi_loc'] = config['isw_Vphi_loc']

    profiles, sqrtspol, sqrtstor = load_cgs_profiles_and_interp(profiles_src, profiles_interp_config)

    multispec['/num_radial_pts'] = len(sqrtspol)
    multispec['/species_def'] = get_species_def_array(config, multispec['/num_radial_pts'])
    multispec['/rel_stages'] = np.full(multispec['/num_radial_pts'], multispec['/species_tag_Vphi'])
    multispec['/rho_pol'] = sqrtspol
    stor = sqrtstor**2
    multispec['/boozer_s'] = stor
    multispec['/T_prof'] = np.column_stack([profiles['Te'], profiles['Ti']])
    multispec['/dT_ov_ds_prof'] = np.column_stack([derivative(stor, profiles['Te']), 
                                                derivative(stor, profiles['Ti'])])
    multispec['/n_prof'] = np.column_stack([profiles['ne'], profiles['ni']])
    multispec['/dn_ov_ds_prof'] = np.column_stack([derivative(stor, profiles['ne']),
                                                derivative(stor, profiles['ni'])])
    multispec['/Vphi'] = profiles['vrot']
    ELEMENTARY_CHARGE_CGS = 1.60217662e-19 * 2.99792458e8 * 10
    charge_e = config['Ze']*ELEMENTARY_CHARGE_CGS
    charge_i = config['Zi']*ELEMENTARY_CHARGE_CGS
    multispec['/kappa_prof'] = np.column_stack([get_kappa(profiles['ne'], profiles['Te'], charge_e), 
                                             get_kappa(profiles['ne'], profiles['Ti'], charge_i)])
    
    write_multispec_to_hdf5(config['hdf5_filename'], multispec)

def get_species_def_array(config, num_radial_pts):
    species_def = np.array([[config['Ze'], config['me']], [config['Zi'], config['mi']]])
    return np.tile(species_def, num_radial_pts)

def derivative(x,y):
    left_border = (-1/2*y[2] + 2*y[1] -3/2*y[0]) / (x[1] - x[0])
    middle = (y[2:] - y[:-2]) / (x[2:] - x[:-2])
    right_border = (3/2*y[-1] - 2*y[-2] + 1/2*y[-3]) / (x[-1] - x[-2])
    return np.concatenate([[left_border], middle, [right_border]])

def get_coulomb_logarithm(n_cgs, T_cgs):
    # Coulomb logarithm (set as species-independent - see E A Belli and J Candy PPCF 54 015015 (2012))
    log_Lambda = 52.43 - 1.15 * np.log10(n_cgs) + 2.3 * np.log10(T_cgs)
    return log_Lambda

def get_kappa(n_cgs, T_cgs, charge_cgs):
    log_Lambda = get_coulomb_logarithm(n_cgs, T_cgs)
    mean_free_path = (3 / (4 * np.sqrt(np.pi)) * (T_cgs/charge_cgs) ** 2 / 
                                    (n_cgs * (charge_cgs ** 2) * log_Lambda))
    kappa = 2 / mean_free_path
    return kappa

def write_multispec_to_hdf5(hdf5_filename, multispec):
    with h5py.File(hdf5_filename, 'w') as f:
        f.create_dataset('/num_radial_pts', data=multispec['/num_radial_pts'], dtype='int32')
        f.create_dataset('/num_species', data=multispec['/num_species'], dtype='int32')
        f.create_dataset('/species_tag', data=multispec['/species_tag'], dtype='int32')
        f.create_dataset('/species_def', data=multispec['/species_def'], dtype='float64')
        f['species_def'].attrs['unit'] = 'e; g'
        f.create_dataset('/boozer_s', data=multispec['/boozer_s'], dtype='float64')
        f['boozer_s'].attrs['unit'] = '1'
        f.create_dataset('/rho_pol', data=multispec['/rho_pol'], dtype='float64')
        f['rho_pol'].attrs['unit'] = '1'
        f.create_dataset('/Vphi', data=multispec['/Vphi'], dtype='float64')
        f['Vphi'].attrs['unit'] = 'rad / s'
        f.create_dataset('/species_tag_Vphi', data=multispec['/species_tag_Vphi'], dtype='int32')
        f.create_dataset('/isw_Vphi_loc', data=multispec['/isw_Vphi_loc'], dtype='int32')
        f.create_dataset('/rel_stages', data=multispec['/rel_stages'], dtype='int32')
        f.create_dataset('/T_prof', data=multispec['/T_prof'], dtype='float64')
        f['T_prof'].attrs['unit'] = 'erg'
        f.create_dataset('/dT_ov_ds_prof', data=multispec['/dT_ov_ds_prof'], dtype='float64')
        f.create_dataset('/n_prof', data=multispec['/n_prof'], dtype='float64')
        f['n_prof'].attrs['unit'] = '1 / cm^3'
        f.create_dataset('/dn_ov_ds_prof', data=multispec['/dn_ov_ds_prof'], dtype='float64')
        f.create_dataset('/kappa_prof', data=multispec['/kappa_prof'], dtype='float64')

    # if hdf5_file_name is None or hdf5_file_name == '':
    #     hdf5_file_name = 'profiles.in'

    # if species_definition is None:
    #     Zi = 1  # ion charge number (deuterium)
    #     mi = 3.3436e-24  # ion mass (deuterium)
    #     Ze = -1  # electron charge number
    #     me = 9.1094e-28  # electron mass
    #     species_definition = np.zeros((gridpoints, 2, 2))
    #     species_definition[0, 0, 0] = Ze
    #     species_definition[0, 0, 1] = me
    #     species_definition[0, 1, 0] = Zi
    #     species_definition[0, 1, 1] = mi
    #     species_definition = np.tile(species_definition, (gridpoints, 1, 1))

    # if isw_Vphi_loc is None:
    #     isw_Vphi_loc = 0

    # if species_tag_Vphi is None:
    #     species_tag_Vphi = 2

    # if input_unit_type is None:
    #     input_unit_type = 1

    # if bounds is None:
    #     gridpoints = species_definition.shape[0]
    # else:
    #     gridpoints = [species_definition.shape[2], bounds[0], bounds[1]]

    # if switch_grid is None:
    #     switch_grid = []

    # num_species = species_definition.shape[0]

    # # UNIT CONVERSION CONSTANTS
    # ELEMENTARY_CHARGE_SI = 1.60217662e-19
    # SPEED_OF_LIGHT_SI = 2.99792458e8
    # EV_TO_SI = ELEMENTARY_CHARGE_SI
    # DENSITY_SI_TO_CGS = 1e-6
    # ENERGY_TO_CGS = 1e7
    # EV_TO_CGS = EV_TO_SI * ENERGY_TO_CGS
    # CHARGE_SI_TO_CGS = 10 * SPEED_OF_LIGHT_SI  # Conversion factor is not speed of light, but 10c_si.

    # # Define profile data
    # rho_pol, rho_tor, ne_si, Ti_eV, Te_eV, vrot = get_profiles_over_equidist_grid(path_to_shot, src, gridpoints, 0, input_unit_type, switch_grid)

    # # Calculate boozer_s
    # boozer_s = rho_tor ** 2

    # # Convert profiles to appropriate units
    # ne = ne_si * DENSITY_SI_TO_CGS
    # ni = ne
    # Te = Te_eV * EV_TO_CGS
    # Ti = Ti_eV * EV_TO_CGS

    # # Calculate profile derivatives
    # def derivate(s, t):
    #     return np.concatenate([[(s[2] - s[0]) / (t[2] - t[0])], (s[2:] - s[:-2]) / (t[2:] - t[:-2]), [(s[-1] - s[-3]) / (t[-1] - t[-3])]])

    # dTe_ov_ds = derivate(Te, boozer_s)
    # dTi_ov_ds = derivate(Ti, boozer_s)
    # dne_ov_ds = derivate(ne, boozer_s)
    # dni_ov_ds = derivate(ni, boozer_s)

    # # Initialize species tags
    # num_radial_pts = len(rho_pol)
    # species_tag = np.arange(1, num_species + 1)

    # rel_stages = np.full(num_radial_pts, species_tag_Vphi)

    # # Coulomb logarithm
    # log_Lambda = 39.1 - 1.15 * np.log10(ne_si) + 2.3 * np.log10(Te_eV * 1e-3)

    # # Electron and ion mean free path
    # e_e = -ELEMENTARY_CHARGE_SI * CHARGE_SI_TO_CGS
    # e_i = +ELEMENTARY_CHARGE_SI * CHARGE_SI_TO_CGS
    # Te_ov_ee = -Te / e_e
    # Ti_ov_ei = Ti / e_i
    # lc_e = (3 / (4 * np.sqrt(np.pi))) * (Te_ov_ee ** 2) / (ne * (e_e ** 2) * log_Lambda)
    # lc_i = (3 / (4 * np.sqrt(np.pi))) * (Ti_ov_ei ** 2) / (ni * (e_i ** 2) * log_Lambda)

    # # Compute kappa
    # kappa_e = 2 / lc_e
    # kappa_i = 2 / lc_i

    # # Write to HDF5
    # with h5py.File(hdf5_file_name, 'w') as f:
    #     f.create_dataset('/num_radial_pts', data=np.int32([num_radial_pts]))
    #     f.create_dataset('/num_species', data=np.int32([num_species]))
    #     f.create_dataset('/species_tag', data=np.int32(species_tag))
    #     f.create_dataset('/species_def', data=species_definition)
    #     f.create_dataset('/boozer_s', data=boozer_s)
    #     f.create_dataset('/rho_pol', data=rho_pol)
    #     f.create_dataset('/Vphi', data=vrot)
    #     f.create_dataset('/species_tag_Vphi', data=np.int32([species_tag_Vphi]))
    #     f.create_dataset('/isw_Vphi_loc', data=np.int32([isw_Vphi_loc]))
    #     f.create_dataset('/rel_stages', data=np.int32(rel_stages))
    #     f.create_dataset('/T_prof', data=np.vstack((Te, Ti)))
    #     f.create_dataset('/dT_ov_ds_prof', data=np.vstack((dTe_ov_ds, dTi_ov_ds)))
    #     f.create_dataset('/n_prof', data=np.vstack((ne, ni)))
    #     f.create_dataset('/dn_ov_ds_prof', data=np.vstack((dne_ov_ds, dni_ov_ds)))
    #     f.create_dataset('/kappa_prof', data=np.vstack((kappa_e, kappa_i)))
