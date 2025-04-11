import numpy as np
import h5py

def generate_den_ti_te_for_mars(neo2ql_input_file):
    sqrtspol, ne, te, ti = get_profiles_from_neo2ql(neo2ql_input_file)
    CM3_to_M3 = 1e-6
    ERG_to_EV = 1.0/(1.6 * 1e-19 * 1e7)
    ne /= CM3_to_M3 # MARS is SI while NEO-2 is cgs
    te *= ERG_to_EV
    ti *= ERG_to_EV
    write_profiles_to_mars_input(ne, te, ti, sqrtspol)

def get_profiles_from_neo2ql(neo2ql_input_file):
    neo2ql = h5py.File(neo2ql_input_file, "r")
    sqrtspol = np.array(neo2ql["rho_pol"])
    temperatures = np.array(neo2ql["T_prof"])
    densities = np.array(neo2ql["n_prof"])
    ne = densities[0,:]
    te = temperatures[0,:]
    ti = temperatures[1,:] # MARS only takes one ion temperature
    return sqrtspol, ne, te, ti

def write_profiles_to_mars_input(ne, te, ti, sqrtspol):
    number_of_surfaces = len(sqrtspol)
    type_of_radial_variable = 1 # 1 means sqrtspol for MARS
    header = f"{number_of_surfaces} {type_of_radial_variable}"
    np.savetxt("PROFDEN.IN", np.array([sqrtspol, ne]).T, header=header, comments="") # MARS assumes ni = ne
    np.savetxt("PROFTE.IN", np.array([sqrtspol, te]).T, header=header, comments="")
    np.savetxt("PROFTI.IN", np.array([sqrtspol, ti]).T, header=header, comments="")