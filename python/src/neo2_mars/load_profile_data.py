import numpy as np

def get_profiles_over_equidist_grid(src:dict, equidist_grid:str, options:dict={}):

    options = set_default_options_where_missing(options)
    profile = load_profile(src)
    if options['is_not_SI']:
        profile = convert_units_to_SI(profile)

    if equidist_grid == 'sqrtspol':
        profile['sqrtspol'], profile['sqrtstor'] = interp_to_equidist_s(profile['sqrtspol'], profile['sqrtstor'], options)
    elif equidist_grid == 'sqrtstor':
        profile['sqrtstor'], profile['sqrtspol'] = interp_to_equidist_s(profile['sqrtstor'], profile['sqrtspol'], options)
    elif equidist_grid == 'spol':
        spol, profile['sqrtstor'] = interp_to_equidist_s(profile['sqrtspol']**2, profile['sqrtstor'], options)
        profile['sqrtspol'] = np.sqrt(spol)

    profile = interp_profile(profile, profile['sqrtspol']**2)

    return profile

def set_default_options_where_missing(options):
    default_options = {'n_s': 30, 'min_s': 0.0, 'max_s': 1.0, 'is_not_SI': False}
    for key in default_options:
        if key not in options:
            options[key] = default_options[key]
    options['max_s'] = min(options['max_s'], 1.0)
    options['min_s'] = max(options['min_s'], 0.0)
    return options

def load_profile(source):
    profiles = {}
    s = {}
    for key in source:
        data = np.loadtxt(source[key]['filename'])
        if key == 'sqrtspol' or key == 'sqrtstor':
            s[key] = data[:,source[key]['column']]
        else:
            profiles[key] = data[:,source[key]['column']]
    return profiles, s

def convert_units_to_SI(profiles):
    keV2eV = 1.0e3
    krads2rads = 1.0e3
    norm_density2si_density = 1.0e+19
    profiles['electron_temperature'] *= keV2eV
    profiles['ion_temperature'] *= keV2eV
    profiles['rotation_velocity'] /= profiles['major_radius']
    profiles['rotation_velocity'] *= krads2rads
    profiles['electron_density'] *= norm_density2si_density

def interp_to_equidist_s(s, y, options):
    equidist_s = np.linspace(options['min_s'], options['max_s'], options['n_s'])
    equidist_y = np.interp(equidist_s, s, y)
    if options['max_s'] == 1.0:
        equidist_y[-1] = 1.0
    if options['min_s'] == 0.0:
        equidist_y[0] = 0.0
    return equidist_s, equidist_y

def interp_profile(profiles, spol):
    for key in profiles:
        if key != 'sqrtspol' or key != 'sqrtstor':
            profiles[key] = np.interp(spol, profiles['sqrtspol']**2, profiles[key])
    return profiles