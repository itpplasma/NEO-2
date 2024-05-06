# %% Standard imports
import numpy as np

def get_integral_ntv_neo2ql(neo2_ql_hdf5):
    volume_element, _ = get_volume_element(neo2_ql_hdf5)
    ntv_density, boozer_s = get_ntv_density(neo2_ql_hdf5)
    integral_ntv = np.trapz(ntv_density * volume_element, boozer_s)
    integral_ntv *= 1e-7  # CGS to SI
    return integral_ntv

def get_volume_element(neo2_ql_hdf5):
    from hdf5tools import get_hdf5file
    with get_hdf5file(neo2_ql_hdf5) as output:
        s = output['boozer_s'][:]
        volume_element = ((2*np.pi)**2 * output['psi_pr_hat'][-1] * output['Bref'][-1] 
                    * (output['aiota'][:] * output['bcovar_tht'][:] + output['bcovar_phi'][:])
                    / (output['avbhat2'][:] * output['Bref'][:]**2))
    return volume_element, s

def get_ntv_density(neo2_ql_hdf5):
    from hdf5tools import get_hdf5file
    with get_hdf5file(neo2_ql_hdf5) as output:
        s = output['boozer_s'][:]
        ntv_density = (output['TphiNA_spec'][:].T)
    return ntv_density, s

if __name__ == "__main__":
    neo2_ql_hdf5 = '/proj/plasma/DATA/DEMO/NEO-2/mars_100kAt_hamada_fluid-n_minus1/neo-2_1_00_rotation/neo2_multispecies_out.h5'
    integral_ntv = get_integral_ntv_neo2ql(neo2_ql_hdf5)
    print('Integral_ntv neo2ql')
    print('------------')
    print('electrons | ions | total')
    print(integral_ntv, np.sum(integral_ntv))