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
  Level 1: adds measured toroidal rotation
  Level 2 (planned): adds neoclassical poloidal velocity estimate
"""

import numpy as np

C_CGS = 2.99792458e10      # speed of light [cm/s]
E_CGS = 4.80320427e-10     # elementary charge [statcoulomb]


def _prepare_common_inputs(
    n, T, dn_ds, dT_ds, z, aiota, sqrtg_bctrvr_phi, av_nabla_stor
):
    n = np.asarray(n)
    T = np.asarray(T)
    dn_ds = np.asarray(dn_ds)
    dT_ds = np.asarray(dT_ds)
    aiota = np.asarray(aiota)
    sqrtg_bctrvr_phi = np.asarray(sqrtg_bctrvr_phi)
    av_nabla_stor = np.asarray(av_nabla_stor)

    if np.any(n == 0.0):
        raise ValueError('n must be nonzero to compute E_r and Om_tE')
    if z == 0.0:
        raise ValueError('z must be nonzero to compute E_r and Om_tE')
    if np.any(aiota == 0.0):
        raise ValueError('aiota must be nonzero to compute Om_tE')
    if np.any(sqrtg_bctrvr_phi == 0.0):
        raise ValueError('sqrtg_bctrvr_phi must be nonzero to compute Om_tE')

    return n, T, dn_ds, dT_ds, aiota, sqrtg_bctrvr_phi, av_nabla_stor


def _compute_diamagnetic_er(n, T, dn_ds, dT_ds, z, av_nabla_stor):
    dp_ds = T * dn_ds + n * dT_ds
    dp_dr = dp_ds * av_nabla_stor
    return dp_dr / (n * z * E_CGS)


def compute_omte_force_balance(
    n,
    T,
    dn_ds,
    dT_ds,
    z,
    aiota,
    sqrtg_bctrvr_phi,
    av_nabla_stor,
    v_phi=None,
    b_theta=None,
):
    """Compute Om_tE from radial force balance.

    Uses the radial force balance:
        E_r = (1 / (Z e n)) dp/dr
            + v_phi * B_theta / c
    and converts to the E x B rotation frequency:
        Om_tE = c E_r / (iota * sqrt(g) B^phi)

    Parameters
    ----------
    n : float or array
        Ion density [1/cm^3] evaluated on the same radial grid as the
        derivatives.
    T : float or array
        Ion temperature [erg] evaluated on the same radial grid as the
        derivatives.
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
    v_phi : float or array, optional
        Toroidal rotation velocity [cm/s]. If omitted together with
        `b_theta`, this falls back to Level 0.
    b_theta : float or array, optional
        Poloidal magnetic field [G]. Must be provided together with `v_phi`
        for Level 1.

    Returns
    -------
    Om_tE : float or array
        E x B rotation frequency [rad/s].
    Er : float or array
        Radial electric field (w.r.t. effective radius) [statV/cm].

    Notes
    -----
    `n`, `T`, `dn_ds`, and `dT_ds` must come from the radial profile input
    used for the run. The scalar species state written to `n_spec`/`T_spec`
    in the multispecies output is not a substitute for those profiles.
    """
    (
        n,
        T,
        dn_ds,
        dT_ds,
        aiota,
        sqrtg_bctrvr_phi,
        av_nabla_stor,
    ) = _prepare_common_inputs(
        n, T, dn_ds, dT_ds, z, aiota, sqrtg_bctrvr_phi, av_nabla_stor
    )

    if (v_phi is None) != (b_theta is None):
        raise ValueError('v_phi and b_theta must be provided together')

    Er = _compute_diamagnetic_er(n, T, dn_ds, dT_ds, z, av_nabla_stor)

    if v_phi is not None:
        Er = Er + np.asarray(v_phi) * np.asarray(b_theta) / C_CGS

    Om_tE = C_CGS * Er / (aiota * sqrtg_bctrvr_phi)
    return Om_tE, Er


def compute_omte_diamagnetic(n, T, dn_ds, dT_ds, z,
                             aiota, sqrtg_bctrvr_phi, av_nabla_stor):
    """Compute Om_tE from diamagnetic pressure gradient (Level 0)."""
    return compute_omte_force_balance(
        n=n,
        T=T,
        dn_ds=dn_ds,
        dT_ds=dT_ds,
        z=z,
        aiota=aiota,
        sqrtg_bctrvr_phi=sqrtg_bctrvr_phi,
        av_nabla_stor=av_nabla_stor,
    )


def compute_omte_toroidal_rotation(
    n,
    T,
    dn_ds,
    dT_ds,
    z,
    aiota,
    sqrtg_bctrvr_phi,
    av_nabla_stor,
    v_phi,
    b_theta,
):
    """Compute Om_tE including measured toroidal rotation (Level 1)."""
    return compute_omte_force_balance(
        n=n,
        T=T,
        dn_ds=dn_ds,
        dT_ds=dT_ds,
        z=z,
        aiota=aiota,
        sqrtg_bctrvr_phi=sqrtg_bctrvr_phi,
        av_nabla_stor=av_nabla_stor,
        v_phi=v_phi,
        b_theta=b_theta,
    )
