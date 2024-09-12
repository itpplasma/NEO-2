import numpy as np

def load_cgs_profiles_and_interp(src:dict, interp_config:dict={}):

    interp_config = set_default_config_where_missing(interp_config)

    profiles, sqrtspol, sqrtstor = load_profiles(src)
    if interp_config['is_not_si']:
        profiles = convert_units_from_norm_to_si(profiles)
    profiles = convert_units_from_si_to_cgs(profiles)

    if interp_config['grid'] == 'sqrtspol':
        sqrtspol, sqrtstor = interp_y_to_equidist_s(s=sqrtspol, y=sqrtstor, options=interp_config)
    elif interp_config['grid'] == 'sqrtstor':
        sqrtstor, sqrtspol = interp_y_to_equidist_s(s=sqrtstor, y=sqrtspol, options=interp_config)
    elif interp_config['grid'] == 'spol':
        spol, sqrtstor = interp_y_to_equidist_s(s=sqrtspol**2, y=sqrtstor, options=interp_config)
        sqrtspol = np.sqrt(spol)
    elif interp_config['grid'] == 'stor':
        stor, sqrtspol = interp_y_to_equidist_s(s=sqrtstor**2, y=sqrtspol, options=interp_config)
        sqrtstor = np.sqrt(stor)
    else:
        raise ValueError('Unknown grid type: ' + interp_config['grid'])

    profiles = interp_profiles(profiles, sqrtspol)
    return profiles, sqrtspol, sqrtstor

def set_default_config_where_missing(options):
    default_options = {'n_s': 10, 'min_s': 0.0, 'max_s': 1.0, 
                       'grid': 'sqrtspol', 'is_not_si': False}
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
            profiles[profile] = data[:,src[profile]['column']].T
    data = np.loadtxt(src['sqrtspol']['filename'])
    sqrtspol = data[:,src['sqrtspol']['column']]
    sqrtstor = data[:,src['sqrtstor']['column']]
    profiles['sqrtspol'] = sqrtspol.copy()
    return profiles, sqrtspol, sqrtstor

def convert_units_from_norm_to_si(profiles):
    KEV2EV = 1.0e3
    KRADS2RADS = 1.0e3
    NORM_DENSITY2SI_DENSITY = 1.0e+19
    profiles['T'] *= KEV2EV
    profiles['vrot'] /= profiles['major_radius']
    profiles['vrot'] *= KRADS2RADS
    profiles['n'] *= NORM_DENSITY2SI_DENSITY
    return profiles

def convert_units_from_si_to_cgs(profiles):
    ELEMENTARY_CHARGE_SI = 1.60217662e-19
    ENERGY_TO_CGS = 1e7
    DENSITY_TO_CGS = 1e-6
    profiles['T'] *= ELEMENTARY_CHARGE_SI * ENERGY_TO_CGS
    profiles['n'] *= DENSITY_TO_CGS
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
    if len(y.shape) == 2:
        f_eval = matrix_interp_cubic(x_new, x, y)
    elif len(y.shape) == 1:
        f_eval = vector_interp_cubic(x_new, x, y)
    else:
        raise ValueError('Unknown shape of y: ' + y.shape)
    return f_eval

def matrix_interp_cubic(x_new, x , y):
    size_f = np.array(y.shape)
    size_f[1] = len(x_new)
    f_eval = np.zeros(size_f)
    for i in range(size_f[0]):
        f_eval[i,:] = vector_interp_cubic(x_new, x, y[i,:])
    return f_eval

def vector_interp_cubic(x_new, x , y):
    from scipy.interpolate import interp1d
    f = interp1d(x, y, kind='cubic')
    f_eval = f(x_new)
    return f_eval