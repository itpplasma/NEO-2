import numpy as np

def load_profiles_and_interp(src:dict, equidist_grid:str, options:dict={}):

    options = set_default_options_where_missing(options)

    profiles, sqrtspol, sqrtstor = load_profiles(src)
    if options['is_not_SI']:
        profiles = convert_units_from_norm_to_SI(profiles)
    profiles = convert_units_from_SI_to_CGS(profiles)

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

def convert_units_from_norm_to_SI(profiles):
    KEV2EV = 1.0e3
    KRADS2RADS = 1.0e3
    NORM_DENSITY2SI_DENSITY = 1.0e+19
    profiles['Te'] *= KEV2EV
    profiles['Ti'] *= KEV2EV
    profiles['vrot'] /= profiles['major_radius']
    profiles['vrot'] *= KRADS2RADS
    profiles['ne'] *= NORM_DENSITY2SI_DENSITY
    return profiles

def convert_units_from_SI_to_CGS(profiles):
    ELEMENTARY_CHARGE_SI = 1.60217662e-19
    ENERGY_TO_CGS = 1e7
    DENSITY_TO_CGS = 1e-6
    profiles['Te'] *= ELEMENTARY_CHARGE_SI * ENERGY_TO_CGS
    profiles['Ti'] *= ELEMENTARY_CHARGE_SI * ENERGY_TO_CGS
    profiles['ne'] *= DENSITY_TO_CGS
    return profiles


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