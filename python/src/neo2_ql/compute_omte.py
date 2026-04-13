"""
Compute the E x B rotation frequency Om_tE from radial force balance.

All quantities are in CGS (Gaussian) units to match NEO-2 conventions:
  - density n: 1/cm^3
  - temperature T: erg
  - charge e: statcoulomb (4.803e-10 esu)
  - magnetic field B: Gauss
  - lengths: cm
  - Om_tE: rad/s

Model levels:
  Level 0 (diamagnetic): pressure gradient only, no flows
  Level 1 (planned): adds measured toroidal rotation
  Level 2 (planned): adds neoclassical poloidal velocity estimate
"""

import numpy as np

C_CGS = 2.99792458e10      # speed of light [cm/s]
E_CGS = 4.80320427e-10     # elementary charge [statcoulomb]


def compute_omte_diamagnetic(n, T, dn_ds, dT_ds, z,
                             aiota, sqrtg_bctrvr_phi, av_nabla_stor):
    """Compute Om_tE from diamagnetic pressure gradient (Level 0).

    Uses the radial force balance with zero flow velocity:
        E_r = (1 / (Z e n)) dp/dr
    and converts to the E x B rotation frequency:
        Om_tE = c E_r / (iota * sqrt(g) B^phi)

    Parameters
    ----------
    n : float or array
        Ion density [1/cm^3].
    T : float or array
        Ion temperature [erg].
    dn_ds : float or array
        Derivative of ion density w.r.t. s_tor [1/cm^3].
    dT_ds : float or array
        Derivative of ion temperature w.r.t. s_tor [erg].
    z : float
        Ion charge number (dimensionless, e.g. 1 for D+).
    aiota : float or array
        Rotational transform iota = 1/q.
    sqrtg_bctrvr_phi : float or array
        sqrt(g) * B^phi contravariant component [G cm].
    av_nabla_stor : float or array
        Flux surface average of |nabla s_tor| [1/cm].

    Returns
    -------
    Om_tE : float or array
        E x B rotation frequency [rad/s].
    Er : float or array
        Radial electric field (w.r.t. effective radius) [statV/cm].
    """
    dp_ds = T * dn_ds + n * dT_ds
    dp_dr = dp_ds * av_nabla_stor
    Er = dp_dr / (n * z * E_CGS)
    Om_tE = C_CGS * Er / (aiota * sqrtg_bctrvr_phi)
    return Om_tE, Er
