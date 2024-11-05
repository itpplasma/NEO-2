# %%
import numpy as np
from scipy.integrate import cumulative_trapezoid

def get_integral_ntv_torque_neo2ql(neo2_ql_hdf5):
    dV, _ = get_volume_element(neo2_ql_hdf5)
    torque_density, boozer_s = get_torque_density(neo2_ql_hdf5)
    integral_torque = cumulative_trapezoid(torque_density * dV, x=boozer_s,
                                                                initial=0)
    integral_torque *= 1e-7  # CGS to SI
    return integral_torque, boozer_s

def get_volume_element(neo2_ql_hdf5):
    from hdf5tools import get_hdf5file
    with get_hdf5file(neo2_ql_hdf5) as output:
        s = output['boozer_s'][:]
        volume_element = ((2*np.pi)**2 * output['psi_pr_hat'][-1] * output['Bref'][-1] 
                    * (output['bcovar_tht'][:]*output['aiota'][:] + output['bcovar_phi'][:])
                    / (output['avbhat2'][:] * output['Bref'][:]**2))
    return volume_element, s

def get_torque_density(neo2_ql_hdf5):
    from hdf5tools import get_hdf5file
    with get_hdf5file(neo2_ql_hdf5) as output:
        s = output['boozer_s'][:]
        ntv_density = (output['TphiNA_spec'][:].T)
    return ntv_density, s
