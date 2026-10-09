"""Convert a radial electric field profile to the NEO-2 input Om_tE.

NEO-2-QL with isw_calc_Er=2 takes the ExB rotation frequency Om_tE (rad/s)
per flux surface. NEO-2 relates it to its radial electric field by

    Om_tE = c * Er / (aiota * sqrtg_bctrvr_phi) = c * Er / sqrtg_bctrvr_tht,

(ntv_mod.f90, compute_Er), i.e. Om_tE = -c dPhi/dpsi_pol. Here
sqrtg_bctrvr_tht = <|grad s|> * dpsi_pol/ds is dpsi_pol/dr for the radial
coordinate dr = ds/<|grad s|>, and Er = -dPhi/dr for the same r, i.e. the
flux-surface average of the normal component E.grad(s)/|grad s|.

With E_r in V/m, sqrtg_bctrvr_tht in G cm and Er_cgs = E_r * 1e4 / c_SI:

    Om_tE = E_r * 1e6 / sqrtg_bctrvr_tht.

The constants are exact; NEO-2 itself uses c = 2.9979e10 cm/s, which differs
from the exact value by 8e-6 relative, so NEO-2 recovers Er from this Om_tE to
that accuracy.

Requirements on the inputs:
- E_r must be the flux-surface-averaged normal field (NEO-2's convention).
  A field given at one poloidal position (e.g. the outboard midplane) differs
  from it by |grad s|(position)/<|grad s|> and must be rescaled first.
- sqrtg_bctrvr_tht is surface geometry that the multi-species HDF5 input does
  not carry. NEO-2 writes it ('sqrtg_bctrvr_tht', together with 'boozer_s')
  to neo2_multispecies_out.h5 on every surface; it does not depend on Er, so
  any run on the same equilibrium and surfaces (any isw_calc_Er) provides it.
"""

import numpy as np

_V_PER_M_TIMES_C_TO_CGS = 1.0e6  # c_cgs * (statV/cm per V/m) = 1e2 * 1e4


def er_to_om_te(er_si, sqrtg_bctrvr_tht):
    """Om_tE in rad/s from E_r in V/m and NEO-2's sqrtg_bctrvr_tht in G cm.

    Both arguments are flux-surface quantities on the same surfaces and may
    be scalars or arrays.
    """
    er_si = np.asarray(er_si, dtype=float)
    geom = np.asarray(sqrtg_bctrvr_tht, dtype=float)
    if np.any(geom == 0.0):
        raise ValueError('sqrtg_bctrvr_tht must be non-zero')
    return _V_PER_M_TIMES_C_TO_CGS * er_si / geom


def om_te_from_neo2_geometry(er_si, s_er, neo2_output_files, s_target):
    """Om_tE on the surfaces s_target from an E_r(s) profile.

    er_si, s_er: E_r in V/m (flux-surface averaged) on the toroidal-flux
        label s (boozer_s); interpolated with a cubic spline.
    neo2_output_files: neo2_multispecies_out.h5 files of a NEO-2 run on the
        same equilibrium, one per surface; they provide sqrtg_bctrvr_tht.
    s_target: surfaces (boozer_s) to evaluate, e.g. the 'boozer_s' of the
        multi-species input. Each must match one output file to 1e-10.
    """
    import h5py
    from scipy.interpolate import CubicSpline

    s_geo, geom = [], []
    for name in neo2_output_files:
        with h5py.File(name, 'r') as f:
            s_geo.append(float(np.ravel(f['boozer_s'][()])[0]))
            geom.append(float(np.ravel(f['sqrtg_bctrvr_tht'][()])[0]))
    s_geo = np.array(s_geo)
    geom = np.array(geom)

    s_target = np.atleast_1d(np.asarray(s_target, dtype=float))
    geom_target = np.empty_like(s_target)
    for k, s in enumerate(s_target):
        match = np.flatnonzero(np.abs(s_geo - s) <= 1e-10)
        if match.size == 0:
            raise ValueError(f'no NEO-2 output for boozer_s = {s}')
        geom_target[k] = geom[match[0]]

    er_target = CubicSpline(np.asarray(s_er, dtype=float),
                            np.asarray(er_si, dtype=float))(s_target)
    return er_to_om_te(er_target, geom_target)
