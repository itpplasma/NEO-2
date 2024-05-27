# %% Standard imports
import numpy as np
import matplotlib.pyplot as plt
import h5py

########################################################################################

def get_neo2_ql_multispec_profiles(multispec_filename: str):
    multispec = h5py.File(multispec_filename,"r",)
    profiles = read_profiles(multispec)
    profiles['species'] = read_species(multispec)
    return profiles

def read_profiles(multispec):
    profiles = {}
    profiles_names = ['T_prof', 'dT_ov_ds_prof', 'n_prof', 'dn_ov_ds_prof', 'kappa_prof', 'Vphi', 'rho_pol']
    for name in profiles_names:
        profiles[name] = {}
        profiles[name]['x'] = np.array(multispec['boozer_s'])
        profiles[name]['y'] = np.array(multispec[name])
    return profiles

def read_species(multispec):
    species = {}
    species['tag'] = np.array(multispec['species_tag'])
    species['charge'] = np.array(multispec['species_def'])[0][:,0]
    species['mass'] = np.array(multispec['species_def'])[1][:,0]
    return species

########################################################################################

def make_figure_neo2_ql_multispec_profiles(profiles, color: str='b', label: str=''):
    fig, axes = plt.subplots(4, 2, figsize=(10, 5 * 4))
    add_profiles_to_axes(profiles, axes, color=color, label=label)
    return fig, axes

def add_profiles_to_axes(profiles, axes, color, label):
    row = 0
    col = 0
    for profile in profiles.keys():
        if profile not in ['boozer_s', 'species']:
            if profiles[profile]['y'].ndim == 1:
                add_one_species_profile_to_ax(profiles[profile], axes[row, col], color)
            else:
                add_multiple_species_profiles_to_ax(profiles[profile], axes[row, col], color)
            set_axes_properties(axes[row, col], profile, color=color, label=label)
            col += 1
            if col == 2:
                col = 0
                row += 1
    return axes

def add_multiple_species_profiles_to_ax(profiles, ax, color):
    species_marker = ['--', ':', '-', '-.']
    for species in range(profiles['y'].shape[0]):
        ax.plot(profiles['x'], profiles['y'][species], species_marker[species], color=color)

def add_one_species_profile_to_ax(profile, ax, color):
    ax.plot(profile['x'], profile['y'], color=color)

def set_axes_properties(ax, profile, color, label):
    ax.set_xlabel(r'$s_\mathrm{tor}$')
    ax.set_ylabel(profile)
    ax.plot([], [], color=color, label=label)
    ax.legend()  

########################################################################################

eps_2023_multispec = "/proj/plasma/DATA/DEMO/teams/Equilibrium_DEMO2019_CHEASE/multi_spec_demo_version_1_2_equidistant_s.in"

eps_2023_profiles = get_neo2_ql_multispec_profiles(eps_2023_multispec)
fig, axes = make_figure_neo2_ql_multispec_profiles(eps_2023_profiles, label='EPS 2023')
plt.show()
