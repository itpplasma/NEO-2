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
filename = 'g30835.3200_ed6'
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

def dtoe(x):
  return float(x.decode('ascii').replace('D', 'E'))

data = np.loadtxt('rhotprof_aug_30835_3.2_ed6.asc', converters={0:dtoe, 1:dtoe})
rho_tor = data[:, 1]
rho_pol = data[:, 0]
psi_pol = psi_a*rho_pol**2  # First column is rho_pol = sqrt(psi_pol)
s_tor = rho_tor**2
data = np.loadtxt('vtprof_aug_30835_3.2.asc')
Vtor = data[:, 1]
data = np.loadtxt('Tiprof_aug_30835_3.2.asc')
Ti = data[:, 1]
Ti[Ti<0] = 0.0
data = np.loadtxt('neprof_aug_30835_3.2.asc')
n = data[:, 1]/1e6

# Correcting n TODO
# nsmo = 200
# #step = np.linspace(0, 1, nsmo)
# #n[:nsmo] = step*n[:nsmo] + (1.0-step)*3.595e13
# n[:nsmo] = n[nsmo]

nspl = InterpolatedUnivariateSpline(psi_pol, n)
nprspl = nspl.derivative()
plt.figure()
plt.plot(rho_pol, n)
plt.xlabel('rho_pol')
plt.title('n')
plt.figure()
plt.plot(rho_pol, nprspl(psi_pol))
plt.ylim([-1e7,1e7])
plt.xlabel('rho_pol')
plt.title('dn/dpsi')
#%%

Vshift = c*Ti*eV/qe*nprspl(psi_pol)/n


data = np.loadtxt('Mtprofile.dat')
s_NEO2 = data[:,0]
Mt_NEO2_spl = InterpolatedUnivariateSpline(s_NEO2, data[:, 1])
R0 = 165.0
vT = np.sqrt(2*Ti*eV/mi)
Om_tE_NEO2 = Mt_NEO2_spl(s_tor)*vT/R0


plt.figure()
plt.plot(rho_pol, Vtor)
plt.plot(rho_pol, Vtor+Vshift)
plt.plot(rho_pol, Om_tE_NEO2, '--')
plt.ylim([-1.1*max(Vtor), 1.1*max(Vtor)])
plt.xlabel('rho_pol')
plt.title('Om_tE')
plt.legend(['Vtor', 'Vtor + Vshift = Om_tE', 'Om_tE NEO2'])
plt.savefig('test_vphi_30835_t3.2.png')

plt.figure()
plt.plot(rho_pol, (Vtor+Vshift)/Om_tE_NEO2)
plt.ylim([0,2])
plt.title('Om_tE/Om_tE_NEO2')


# %%
# TODO: check k coefficient in NEO-2
# look at /temp/andfmar/NTV_050615_AUG30835rmp_Profile_Ions_lag7/RUN/aug_2_rmp_n1
#


# Poloidal connection length L_c = 2*pi*q*R, we use pi*R often
# Mean free path: l_c = v_T/nu_coll
# Banana regime: (L_c/l_c) * A^(3/2) << 1
# doesnt work for potatoes, banana regime can be down to the axis
