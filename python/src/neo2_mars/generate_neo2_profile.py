import h5py
import numpy as np

from .load_profile_data import load_profile_data

def generate_neo2_profile(hdf5FileName=None, path_to_shot=None, data_source=None, species_definition=None, isw_Vphi_loc=None, species_tag_Vphi=None, input_unit_type=None, bounds=None, switch_grid=None):
    if hdf5FileName is None or hdf5FileName == '':
        hdf5FileName = 'multi_spec_Valentin.in'
    if path_to_shot is None or path_to_shot == '':
        path_to_shot = 'SHOT35568/'

    if data_source is None:
        shot_designation = '35568_t2.6881'
        data_source = {
            'rhopoloidal': {'filename': f'ne_ida_{shot_designation}_rhopol.dat', 'column': 1},
            'rhotoroidal': {'filename': f'ne_ida_{shot_designation}_rhotor.dat', 'column': 1},
            'electron_density': {'filename': f'ne_ida_{shot_designation}_rhopol.dat', 'column': 2},
            'electron_temperature': {'filename': f'Te_ida_{shot_designation}_rhopol.dat', 'column': 2},
            'ion_temperature': {'filename': f'Ti_cez_{shot_designation}_rhopol.dat', 'column': 2},
            'rotation_velocity': {'filename': f'vrot_cez_{shot_designation}_rhopol.dat', 'column': 2},
            'major_radius': {'filename': f'vrot_cez_{shot_designation}_R.dat', 'column': 1}
        }

    if species_definition is None:
        Zi = 1  # ion charge number (deuterium)
        mi = 3.3436e-24  # ion mass (deuterium)
        Ze = -1  # electron charge number
        me = 9.1094e-28  # electron mass
        species_definition = np.zeros((gridpoints, 2, 2))
        species_definition[0, 0, 0] = Ze
        species_definition[0, 0, 1] = me
        species_definition[0, 1, 0] = Zi
        species_definition[0, 1, 1] = mi
        species_definition = np.tile(species_definition, (gridpoints, 1, 1))

    if isw_Vphi_loc is None:
        isw_Vphi_loc = 0

    if species_tag_Vphi is None:
        species_tag_Vphi = 2

    if input_unit_type is None:
        input_unit_type = 1

    if bounds is None:
        gridpoints = species_definition.shape[0]
    else:
        gridpoints = [species_definition.shape[0], bounds[0], bounds[1]]

    if switch_grid is None:
        switch_grid = []

    num_species = species_definition.shape[1]

    # UNIT CONVERSION CONSTANTS
    ELEMENTARY_CHARGE_SI = 1.60217662e-19
    SPEED_OF_LIGHT_SI = 2.99792458e8
    EV_TO_SI = ELEMENTARY_CHARGE_SI
    DENSITY_SI_TO_CGS = 1e-6
    ENERGY_TO_CGS = 1e7
    EV_TO_CGS = EV_TO_SI * ENERGY_TO_CGS
    CHARGE_SI_TO_CGS = 10 * SPEED_OF_LIGHT_SI  # Conversion factor is not speed of light, but 10c_si.

    # Define profile data
    rho_pol, rho_tor, ne_si, Ti_eV, Te_eV, vrot = load_profile_data(path_to_shot, data_source, gridpoints, 0, input_unit_type, switch_grid)

    # Calculate boozer_s
    boozer_s = rho_tor ** 2

    # Convert profiles to appropriate units
    ne = ne_si * DENSITY_SI_TO_CGS
    ni = ne
    Te = Te_eV * EV_TO_CGS
    Ti = Ti_eV * EV_TO_CGS

    # Calculate profile derivatives
    def derivate(s, t):
        return np.concatenate([[(s[2] - s[0]) / (t[2] - t[0])], (s[2:] - s[:-2]) / (t[2:] - t[:-2]), [(s[-1] - s[-3]) / (t[-1] - t[-3])]])

    dTe_ov_ds = derivate(Te, boozer_s)
    dTi_ov_ds = derivate(Ti, boozer_s)
    dne_ov_ds = derivate(ne, boozer_s)
    dni_ov_ds = derivate(ni, boozer_s)

    # Initialize species tags
    num_radial_pts = len(rho_pol)
    species_tag = np.arange(1, num_species + 1)

    rel_stages = np.full(num_radial_pts, species_tag_Vphi)

    # Coulomb logarithm
    log_Lambda = 39.1 - 1.15 * np.log10(ne_si) + 2.3 * np.log10(Te_eV * 1e-3)

    # Electron and ion mean free path
    e_e = -ELEMENTARY_CHARGE_SI * CHARGE_SI_TO_CGS
    e_i = +ELEMENTARY_CHARGE_SI * CHARGE_SI_TO_CGS
    Te_ov_ee = -Te / e_e
    Ti_ov_ei = Ti / e_i
    lc_e = (3 / (4 * np.sqrt(np.pi))) * (Te_ov_ee ** 2) / (ne * (e_e ** 2) * log_Lambda)
    lc_i = (3 / (4 * np.sqrt(np.pi))) * (Ti_ov_ei ** 2) / (ni * (e_i ** 2) * log_Lambda)

    # Compute kappa
    kappa_e = 2 / lc_e
    kappa_i = 2 / lc_i

    # Write to HDF5
    with h5py.File(hdf5FileName, 'w') as f:
        f.create_dataset('/num_radial_pts', data=np.int32([num_radial_pts]))
        f.create_dataset('/num_species', data=np.int32([num_species]))
        f.create_dataset('/species_tag', data=np.int32(species_tag))
        f.create_dataset('/species_def', data=species_definition)
        f.create_dataset('/boozer_s', data=boozer_s)
        f.create_dataset('/rho_pol', data=rho_pol)
        f.create_dataset('/Vphi', data=vrot)
        f.create_dataset('/species_tag_Vphi', data=np.int32([species_tag_Vphi]))
        f.create_dataset('/isw_Vphi_loc', data=np.int32([isw_Vphi_loc]))
        f.create_dataset('/rel_stages', data=np.int32(rel_stages))
        f.create_dataset('/T_prof', data=np.vstack((Te, Ti)))
        f.create_dataset('/dT_ov_ds_prof', data=np.vstack((dTe_ov_ds, dTi_ov_ds)))
        f.create_dataset('/n_prof', data=np.vstack((ne, ni)))
        f.create_dataset('/dn_ov_ds_prof', data=np.vstack((dne_ov_ds, dni_ov_ds)))
        f.create_dataset('/kappa_prof', data=np.vstack((kappa_e, kappa_i)))
