#
# Om_tE = c E_r / psi_pol'(r) = c E_r / (sqrtg B^th)
# Vtor = Om_tE + cT/e n'(r)/(n(r) psi_pol'(r))
#
#
# If radial variable is psi_pol = psi then
#
# Om_tE = c E_psi
# Vtor = Om_tE + cT/e n'(psi)/n(psi)
#

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline
import matplotlib.pyplot as plt

# Constants in CGS units
c = 2.997925e10            # speed of light
qe = 4.803204e-10          # elementary charge
eV = 1.602176e-12          # 1 electron volt
mi = 2*1.66054e-24         # deuterium mass

#%%
filename = 'eqdsk_35568_2.68800.dat'
# Read in the file.
with open(filename) as f:
    lines = f.readlines()

# Split the first line, which contains some usefull (e.g. nw, nh) and
# some not so usefull information (e.g. unknownvaluestr and signtodiscard).
unknownvaluestr, datestr, signtodiscard, shootnr, idum,nw,nh = lines[0].split()

# Variables are so far all strings, but change to integer for some of them.
shootnr = int(shootnr)
idum = int(idum)# Some
nw = int(nw)# Number of horizontal R grid points
nh = int(nh)# Number of vertical z grid points

d = 15
# First few (four) lines contain specific values.
rdim, zdim, rcentr, rleft, zmid = [lines[1][i:i+16] for i in range(0,5*16,16)]
rmaxis, zmaxis, simag, sibry, bcentr = [lines[2][i:i+16] for i in range(0,5*16,16)]
current, simag2, dummy1, rmaxis, dummy2 = [lines[3][i:i+16] for i in range(0,5*16,16)]
zmaxis, dummy3, sibry2, dummy4, dummy5 = [lines[4][i:i+16] for i in range(0,5*16,16)]

# Create array of psi points.
psi_a = -(float(sibry)-float(simag))*1e4*1e4
psi = np.linspace(0, psi_a, nw)


data = np.loadtxt('flux_coordinates_densities_temperatures.dat', skiprows=1)
s_tor = data[1:, 0]


rho_tor = data[1:, 1]
rho_pol = data[1:, 2]
psi_pol = psi_a*rho_pol**2  # First column is rho_pol = sqrt(psi_pol)
n = data[1:, 3]
nspl = UnivariateSpline(psi_pol, n)
nprspl = nspl.derivative()
Te = data[1:, 5]
Ti = data[1:, 6]
Te[Te<0] = 0.0
Ti[Ti<0] = 0.0

plt.figure()
plt.plot(rho_pol, n)
plt.xlabel(r'$\rho_\mathrm{pol}$')
plt.title('Density $n [1/cm^3]$')


plt.figure()
plt.plot(rho_pol, Te)
plt.plot(rho_pol, Ti)
plt.xlabel(r'$\rho_\mathrm{pol}$')
plt.title('Temperature $T [eV]$')
plt.legend(['Electrons', 'Ions'])

#%% read velocities and rotation frequencies and compare

# This velocity comes in km/s . We convert to cm/s
data = np.loadtxt('vrot_cez_35568_t2.6881_rhopol.dat', skiprows=2)
rho_pol_vrot = data[:, 0]
vrot = data[:, 1]*1e3*1e2
vrot_min = data[:, 2]*1e3*1e2
vrot_max = data[:, 3]*1e3*1e2

# Contravariant rotation velocity (frequency)
data = np.loadtxt('vrot_cez_35568_t2.6881_R.dat', skiprows=2)
R = data[:, 0]*1e2
Vphi = vrot/R
Vphi_min = vrot_min/R
Vphi_max = vrot_max/R
R0 = 165.0

R_spl = UnivariateSpline(rho_pol_vrot, R)
Vphi_spl = UnivariateSpline(rho_pol_vrot, Vphi)

plt.figure()
plt.plot(rho_pol, Vphi_spl(rho_pol))
plt.xlabel(r'$\rho_\mathrm{pol}$')
plt.title(r'Toroidal rotation frequency $V^\mathrm{\varphi}$ (rad/s)')

data = np.loadtxt('axisymmetric_quantities.dat', skiprows=1)
MtOvR_NEO2 = data[:, -1]
vT = np.sqrt(2*Ti*eV/mi)
Om_tE_NEO2 = MtOvR_NEO2*vT

# %%


Vshift = c*Ti*eV/qe*nprspl(psi_pol)/n
Vphi_plt = Vphi_spl(rho_pol)

plt.figure()
plt.plot(rho_pol, Vphi_plt)
plt.plot(rho_pol, Vphi_plt+Vshift)
plt.plot(rho_pol, Om_tE_NEO2, '--')
plt.errorbar(rho_pol_vrot, Vphi, (Vphi_max-Vphi_min), linestyle='', color='k')
plt.ylim([1.1*min(Vphi_plt+Vshift), 1.1*max(Vphi_plt+Vshift)])
plt.xlabel('rho_pol')
plt.title('Om_tE')
plt.legend(['Vtor', 'Vtor + Vshift = Om_tE', 'Om_tE NEO2'])

plt.savefig('test_vphi_35568_t2.6881.png')

# %%
# TODO: check k coefficient in NEO-2
# look at /temp/andfmar/NTV_050615_AUG30835rmp_Profile_Ions_lag7/RUN/aug_2_rmp_n1
#


# Poloidal connection length L_c = 2*pi*q*R, we use pi*R often
# Mean free path: l_c = v_T/nu_coll
# Banana regime: (L_c/l_c) * A^(3/2) << 1
# doesnt work for potatoes, banana regime can be down to the axis
