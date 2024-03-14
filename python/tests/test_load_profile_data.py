# %%Standard imports
from neo2_mars import load_profile_data

MARS_DIR = '/proj/plasma/DATA/DEMO/teams/MARSQ_OUTPUTS/DATA_equil_qmod/MARSQ_OUTPUTS_100kAt_dBkinetic_NTVkinetic_NEO2profs'
MARS_DIR = '/proj/plasma/DATA/DEMO/MARS/MARSQ_INPUTS_KNTV21_NEO2profs_RUN/'
OUTPUT_DIR = './processed_mars_inputs/'

from omfit_classes.omfit_mars import OMFITmars
data = OMFITmars(MARS_DIR)
data.load()
q_prof = {}
q_prof['values'] = data['PROFEQ']['q_prof'].values
q_prof['sqrt_spol'] = data['PROFEQ']['q_prof'].coords['s_eq'].values

from libneo import read_eqdsk
eqdsk_ASTRA_1 = read_eqdsk('/proj/plasma/DATA/DEMO/teams/ASTRA/BASELINE_2019/eqdska.gsef')
eqdsk_ASTRA_2 = read_eqdsk('/proj/plasma/DATA/DEMO/teams/ASTRA/BASELINE_2019/Correct_ASTRA_eqdsk_2PCR4X_v1_0.dat')
eqdsk_CHEASE = read_eqdsk('/proj/plasma/DATA/DEMO/teams/Equilibrium_DEMO2019_CHEASE/MOD_Qprof_Test/EQDSK_DEMO2019_q1_COCOS_02.OUT')

import matplotlib.pyplot as plt
plt.figure()
plt.plot(eqdsk_ASTRA_1['rho_poloidal'], eqdsk_ASTRA_1['qprof'], label='eqdska.gsef')
plt.plot(eqdsk_ASTRA_2['rho_poloidal'], eqdsk_ASTRA_2['qprof'], label='Correct_ASTRA_eqdsk_2PCR4X_v1_0.dat')
plt.plot(q_prof['sqrt_spol'], q_prof['values'], label='MARS')
plt.plot(np.sqrt(eqdsk_CHEASE['rho_poloidal']), eqdsk_CHEASE['qprof'], '--',label='EQDSK_DEMO2019_q1_COCOS_02.OUT')
plt.legend()

#%%
from neo2_mars import write_input_for_generate_neo2_profile_from_mars
write_input_for_generate_neo2_profile_from_mars(MARS_DIR, OUTPUT_DIR)
