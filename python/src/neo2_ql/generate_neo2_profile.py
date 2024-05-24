import h5py
import numpy as np

from .load_profile_data import load_profiles_and_interp

def generate_neo2_profile(hdf5_file_name: str, src: dict, options: dict):
    
    return

def derivative(x,y):
    left_border = (-1/2*y[2] + 2*y[1] -3/2*y[0]) / (x[1] - x[0])
    middle = (y[2:] - y[:-2]) / (x[2:] - x[:-2])
    right_border = (3/2*y[-1] - 2*y[-2] + 1/2*y[-3]) / (x[-1] - x[-2])
    return np.concatenate([[left_border], middle, [right_border]])

def get_coulomb_logarithm(n_CGS, T_CGS):
    # Coulomb logarithm (set as species-independent - see E A Belli and J Candy PPCF 54 015015 (2012))
    log_Lambda = 52.43 - 1.15 * np.log10(n_CGS) + 2.3 * np.log10(T_CGS)
    return log_Lambda

def get_kappa(n_CGS, T_CGS, charge_CGS):
    log_Lambda = get_coulomb_logarithm(n_CGS, T_CGS)
    mean_free_path = (3 / (4 * np.sqrt(np.pi)) * (T_CGS/charge_CGS) ** 2 / 
                                    (n_CGS * (charge_CGS ** 2) * log_Lambda))
    kappa = 2 / mean_free_path
    return kappa

def write_neo2ql_inputs_to_hdf5(hdf5_filename, input):
    with h5py.File(hdf5_filename, 'w') as f:
        f.create_dataset('/num_radial_pts', data=input['/num_radial_pts'], dtype='int32')
        f.create_dataset('/num_species', data=input['/num_species'], dtype='int32')
        f.create_dataset('/species_tag', data=input['/species_tag'], dtype='int32')
        f.create_dataset('/species_def', data=input['/species_def'], dtype='float64')
        f['species_def'].attrs['unit'] = 'e; g'
        f.create_dataset('/boozer_s', data=input['/boozer_s'], dtype='float64')
        f['boozer_s'].attrs['unit'] = '1'
        f.create_dataset('/rho_pol', data=input['/rho_pol'], dtype='float64')
        f['rho_pol'].attrs['unit'] = '1'
        f.create_dataset('/Vphi', data=input['/Vphi'], dtype='float64')
        f['Vphi'].attrs['unit'] = 'rad / s'
        f.create_dataset('/species_tag_Vphi', data=input['/species_tag_Vphi'], dtype='int32')
        f.create_dataset('/isw_Vphi_loc', data=input['/isw_Vphi_loc'], dtype='int32')
        f.create_dataset('/rel_stages', data=input['/rel_stages'], dtype='int32')
        f.create_dataset('/T_prof', data=input['/T_prof'], dtype='float64')
        f['T_prof'].attrs['unit'] = 'erg'
        f.create_dataset('/dT_ov_ds_prof', data=input['/dT_ov_ds_prof'], dtype='float64')
        f.create_dataset('/n_prof', data=input['/n_prof'], dtype='float64')
        f['n_prof'].attrs['unit'] = '1 / cm^3'
        f.create_dataset('/dn_ov_ds_prof', data=input['/dn_ov_ds_prof'], dtype='float64')
        f.create_dataset('/kappa_prof', data=input['/kappa_prof'], dtype='float64')

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
