# %% Standard imports
import numpy as np
import matplotlib.pyplot as plt
import h5py

########################################################################################

def get_input_profiles_neo2_ql(input_profiles_h5filename: str):
    input_neo2 = h5py.File(input_profiles_h5filename,"r",)
    input_profiles = read_input_profiles(input_neo2)
    input_profiles['species'] = read_input_species(input_neo2)
    return input_profiles

def read_input_profiles(input_neo2):
    input_profiles = {}
    input_profiles_names = ['boozer_s', 'T_prof', 'dT_ov_ds_prof', 'n_prof', 'dn_ov_ds_prof', 'kappa_prof', 'Vphi', 'rho_pol']
    for name in input_profiles_names:
        input_profiles[name] = np.array(input_neo2[name])
    return input_profiles

def read_input_species(input_neo2):
    input_species = {}
    input_species['tag'] = np.array(input_neo2['species_tag'])
    input_species['charge'] = np.array(input_neo2['species_def'])[0][:,0]
    input_species['mass'] = np.array(input_neo2['species_def'])[1][:,0]
    return input_species

########################################################################################

def make_figure_of_input_profiles_neo2_ql(input_profiles):
    figure, axes = plt.subplots(4, 2, figsize=(10, 5 * 4))
    add_profiles_to_figure(input_profiles, axes,color='b')
    return figure, axes

def add_profiles_to_figure(input_profiles, figure_axes,color):
    row = 0
    col = 0
    for profile_name in input_profiles:
        if profile_name not in ['boozer_s', 'species']:
            if input_profiles[profile_name].ndim == 1:
                add_one_species_profile_to_axes(input_profiles, profile_name, figure_axes[row, col],color)
            else:
                add_multiple_species_profiles_to_axes(input_profiles, profile_name, figure_axes[row, col],color)
            set_axes_properties(figure_axes[row, col], profile_name)
            row, col = update_row_col(row, col)
    return figure_axes

def add_multiple_species_profiles_to_axes(input_profiles, profile_name, axes,color):
    for species in input_profiles['species']['tag']:
        axes.plot(input_profiles['boozer_s'], input_profiles[profile_name][species-1],color=color)

def add_one_species_profile_to_axes(input_profiles, profile_name, axes,color):
    axes.plot(input_profiles['boozer_s'], input_profiles[profile_name],color=color)

def set_axes_properties(axes, profile_name):
    axes.set_xlabel('s_tor')
    axes.set_ylabel(profile_name)
    #if profile_name in ['dT_ov_ds_prof', 'dn_ov_ds_prof']:
        #axes.set_xscale('log')
    #if profile_name in ['T_prof', 'n_prof', 'kappa_prof']:
        #axes.set_yscale('log')
    #axes.legend()

def update_row_col(row, col):
    col += 1
    if col == 2:
        col = 0
        row += 1
    return row, col

########################################################################################

mars_input_profiles_h5filename = "/temp/grassl_g/TEST_NTV_DEMO/multi_spec_demo_MARS.in"
control_input_profiles_h5filename = "/temp/grassl_g/TEST_NTV_DEMO/control_profile.in"
astra_input_profiles_h5filename = "/temp/grassl_g/TEST_NTV_DEMO/multi_spec_demo_ASTRA.in"
astra2_input_profiles_h5filename = "/temp/grassl_g/TEST_NTV_DEMO/multi_spec_demo_ASTRA2.in"
chease_input_profiles_h5filename = "/temp/grassl_g/TEST_NTV_DEMO/multi_spec_demo_CHEASE.in"
temp_buchholz_paper_input_profiles_h5filename = "/proj/plasma/DATA/DEMO/teams/Equilibrium_DEMO2019_CHEASE/multi_spec_demo_version_1_2_equidistant_s.in"

mars_input_profiles = get_input_profiles_neo2_ql(mars_input_profiles_h5filename)
control_input_profiles = get_input_profiles_neo2_ql(control_input_profiles_h5filename)
astra_input_profiles = get_input_profiles_neo2_ql(astra_input_profiles_h5filename)
astra2_input_profiles = get_input_profiles_neo2_ql(astra2_input_profiles_h5filename)
chease_input_profiles = get_input_profiles_neo2_ql(chease_input_profiles_h5filename)
temp_buchholz_paper_input_profiles = get_input_profiles_neo2_ql(temp_buchholz_paper_input_profiles_h5filename)
#control_input_profiles['boozer_s'] = control_input_profiles['boozer_s']**2

profile_figure, profile_axes = make_figure_of_input_profiles_neo2_ql(mars_input_profiles)
#profile_axes = add_profiles_to_figure(control_input_profiles, profile_axes,color='r')
#profile_axes = add_profiles_to_figure(astra_input_profiles, profile_axes,color='g')
#profile_axes = add_profiles_to_figure(astra2_input_profiles, profile_axes,color='y')
#profile_axes = add_profiles_to_figure(chease_input_profiles, profile_axes,color='c')
profile_axes = add_profiles_to_figure(temp_buchholz_paper_input_profiles, profile_axes,color='m')
plt.show()