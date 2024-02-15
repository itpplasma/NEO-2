# %%
import numpy as np
import matplotlib.pyplot as plt
from omfit_classes.omfit_mars import OMFITmars

from libneo import FluxConverter, read_eqdsk

########################################################################################

def get_mars_q_profile_over_stor(path_to_mars_folder: str):
    mars_input_profiles = get_mars_input_profiles(path_to_mars_folder)
    stor = convert_sqrt_spol_to_stor(mars_input_profiles['s_eq'],mars_input_profiles['q_prof'])
    return stor, mars_input_profiles['q_prof']

def get_mars_input_profiles(path_to_mars_folder: str) -> dict:
    data = OMFITmars(path_to_mars_folder)
    data.load()
    return data['PROFEQ']

def convert_sqrt_spol_to_stor(sqrt_spol,q_profil):
    spol = sqrt_spol**2
    _, q_profile_over_equidist_spol = get_q_profile_over_equidist_spol(spol, q_profil)
    converter = FluxConverter(q_profile_over_equidist_spol)
    stor = converter.spol2stor(spol)
    return stor

def get_q_profile_over_equidist_spol(spol, q_profil):
    equidist_spol = np.linspace(np.min(spol), np.max(spol), spol.shape[0])
    q_profile_over_equidist_spol = np.interp(equidist_spol, spol, q_profil)
    return equidist_spol, q_profile_over_equidist_spol

########################################################################################

def get_eqdsk_q_profile_over_stor(filename_eqdsk: str):
    eqdsk = read_eqdsk(filename_eqdsk)
    q_profile = eqdsk['qprof']
    converter = FluxConverter(q_profile)
    stor = converter.spol2stor(eqdsk['rho_poloidal'])
    return stor, q_profile

def get_eqdsk_q_profile_over_sqrt_spol(filename_eqdsk: str) -> dict:
    eqdsk = read_eqdsk(filename_eqdsk)
    q_profile = eqdsk['qprof']
    converter = FluxConverter(q_profile)
    sqrt_spol = np.sqrt(eqdsk['rho_poloidal'])
    return sqrt_spol, q_profile

########################################################################################

def convert_sqrt_spol_to_sqrt_stor(sqrt_spol,q_profil):
    return np.sqrt(convert_sqrt_spol_to_stor(sqrt_spol,q_profil))

########################################################################################

path_to_mars_folder = "/proj/plasma/DATA/DEMO/MARS/MARSQ_INPUTS_KNTV21_NEO2profs_RUN/"
mars_input_profiles = get_mars_input_profiles(path_to_mars_folder)
stor_mars = convert_sqrt_spol_to_stor(mars_input_profiles['s_eq'],mars_input_profiles['q_prof'])

filename_eqdsk = '/proj/plasma/DATA/DEMO/teams/Equilibrium_DEMO2019_CHEASE/MOD_Qprof_Test/EQDSK_DEMO2019_q1_COCOS_02.OUT'
stor_eqdsk, q_profile_eqdsk = get_eqdsk_q_profile_over_stor(filename_eqdsk)
sqrt_spol_eqdsk, _ = get_eqdsk_q_profile_over_sqrt_spol(filename_eqdsk)
