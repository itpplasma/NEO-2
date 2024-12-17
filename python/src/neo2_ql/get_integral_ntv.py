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
    from neo2_util import get_hdf5file, get_hdf5dataset_value
    with get_hdf5file(neo2_ql_hdf5) as output:
        s = get_hdf5dataset_value(output, 'boozer_s')
        Bref = get_hdf5dataset_value(output, 'Bref')
        psi_tor_edge = get_hdf5dataset_value(output, 'psi_pr_hat') * Bref
        bcovar_tht = get_hdf5dataset_value(output, 'bcovar_tht')
        aiota = get_hdf5dataset_value(output, 'aiota')
        bcovar_phi = get_hdf5dataset_value(output, 'bcovar_phi')
        average_b2 = get_hdf5dataset_value(output, 'avbhat2') * Bref**2

        dV_ds = ((2*np.pi)**2 * psi_tor_edge * (bcovar_tht * aiota + bcovar_phi)
                / average_b2)
    return dV_ds, s

def get_torque_density(neo2_ql_hdf5):
    from neo2_util import get_hdf5file
    with get_hdf5file(neo2_ql_hdf5) as output:
        s = output['boozer_s'][:]
        ntv_density = (output['TphiNA_spec'][:].T)
    return ntv_density, s
