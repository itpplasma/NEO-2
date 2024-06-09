#%% Stadart imports
import matplotlib.pyplot as plt
import os

# Homebrew imports
from neo2_mars import mars_sqrtspol2stor

# Modules to test
from neo2_mars import get_omegate_from_neo2ql

test_neo2ql_output_file = "/itp/MooseFS/grassl_g/comparison_vary_coilwidth_MARS_NEO2/lagrange_faulty_profiles_run_001/neo2_multispecies_out.h5"
test_mars_dir = "/proj/plasma/DATA/DEMO/MARS/MARSQ_OUTPUTS_100kAt_dBkinetic_NTVkinetic_NEO2profs_KEYTORQ_1"

def test_get_omegate_from_neo2ql_visual_check():
    omegate_neo2ql, stor_neo2ql = get_omegate_from_neo2ql(test_neo2ql_output_file)
    omegate_mars, sqrtspol_mars = get_mars_omegate(test_mars_dir)
    stor_mars = mars_sqrtspol2stor(test_mars_dir, sqrtspol_mars)
    plt.figure()
    plt.plot(stor_neo2ql,omegate_neo2ql[:,0],'-b', label='NEO-2-QL species 1')
    plt.plot(stor_neo2ql,omegate_neo2ql[:,1],'--g', label='NEO-2-QL species 2')
    plt.plot(stor_mars,omegate_mars,'--r', label='MARS')
    plt.xlabel(f"$s_\mathrm{{tor}}$ [1]")
    plt.ylabel(f"$\Omega_\mathrm{{tE}}$ [1/s]")
    plt.ylim(0,1e4)
    plt.show()

def get_mars_omegate(mars_dir):
    omegate_file = os.path.join(mars_dir, "PROFWE.IN")
    omegate = np.loadtxt(omegate_file, skiprows=1)
    sqrtspol = omegate[:,0]
    omegate = omegate[:,1:]*100
    return omegate, sqrtspol

if __name__ == "__main__":
    test_get_omegate_from_neo2ql_visual_check()