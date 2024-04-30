import numpy as np

def load_profiles_and_interp(src:dict, equidist_grid:str, options:dict={}):

    options = set_default_options_where_missing(options)

    profiles, sqrtspol, sqrtstor = load_profiles(src)
    if options['is_not_SI']:
        profiles = convert_units_to_SI(profiles)

    if equidist_grid == 'sqrtspol':
        sqrtspol, sqrtstor = interp_y_to_equidist_s(s=sqrtspol, y=sqrtstor, options=options)
    elif equidist_grid == 'sqrtstor':
        sqrtstor, sqrtspol = interp_y_to_equidist_s(s=sqrtstor, y=sqrtspol, options=options)
    elif equidist_grid == 'spol':
        spol, sqrtstor = interp_y_to_equidist_s(s=sqrtspol**2, y=sqrtstor, options=options)
        sqrtspol = np.sqrt(spol)

    profiles = interp_profiles(profiles, sqrtspol)

    return profiles, sqrtspol, sqrtstor

def set_default_options_where_missing(options):
    default_options = {'n_s': 10, 'min_s': 0.0, 'max_s': 1.0, 'is_not_SI': False}
    for key in default_options:
        if key not in options:
            options[key] = default_options[key]
    options['max_s'] = min(options['max_s'], 1.0)
    options['min_s'] = max(options['min_s'], 0.0)
    return options

def load_profiles(src):
    profiles = {}
    for profile in src:
        data = np.loadtxt(src[profile]['filename'])
        if profile != 'sqrtspol' and profile != 'sqrtstor':
            profiles[profile] = data[:,src[profile]['column']]
    data = np.loadtxt(src['sqrtspol']['filename'])
    sqrtspol = data[:,src['sqrtspol']['column']]
    sqrtstor = data[:,src['sqrtstor']['column']]
    profiles['sqrtspol'] = sqrtspol.copy()
    return profiles, sqrtspol, sqrtstor

def convert_units_to_SI(profiles):
    keV2eV = 1.0e3
    krads2rads = 1.0e3
    norm_density2si_density = 1.0e+19
    profiles['Te'] *= keV2eV
    profiles['Ti'] *= keV2eV
    profiles['vrot'] /= profiles['major_radius']
    profiles['vrot'] *= krads2rads
    profiles['ne'] *= norm_density2si_density

def interp_y_to_equidist_s(s, y, options):
    equidist_s = np.linspace(options['min_s'], options['max_s'], options['n_s'])
    equidist_y = interp_cubic(equidist_s, s, y)
    if options['max_s'] == 1.0:
        equidist_y[-1] = 1.0
    if options['min_s'] == 0.0:
        equidist_y[0] = 0.0
    return equidist_s, equidist_y

def interp_profiles(profiles, sqrtspol):
    for profile in profiles:
        if profile != 'sqrtspol':
            profiles[profile] = interp_cubic(sqrtspol, profiles['sqrtspol'], profiles[profile])
    del profiles['sqrtspol']
    return profiles

def interp_cubic(x_new, x , y):
    from scipy.interpolate import interp1d
    f = interp1d(x, y, kind='cubic')
    return f(x_new)