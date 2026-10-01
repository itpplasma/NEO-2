"""
Radial force-balance estimates of E_r and Om_tE, from simple to full NEO-2.

This module implements the benchmark ramp of issue itpplasma/NEO-2#75. Each
level adds one physical ingredient to the radial force balance of the
rotation species ``i`` (the species selected by ``species_tag_Vphi``):

Level 0  diamagnetic only (no flows)                 ``er_level0_diamagnetic``
Level 1  + rigid toroidal rotation, v_theta = 0       ``er_level1_toroidal_rotation``
Level 2  + neoclassical poloidal rotation, given k    ``er_level2_poloidal_rotation``
Level 3  full multi-species NEO-2 closure             ``er_level3_neo2_multispecies``

Level 3 replays ``compute_Er`` in ``NEO-2-QL/ntv_mod.f90`` for all
``isw_Vphi_loc`` branches (0: flux-surface average, 1, 2: Vphi given at one
point). It reproduces the stored Er up to roundoff and, in general, the
residual of NEO-2's fixed-point iteration for (Er, <E_par B>/<B^2>) (relative
tolerance ``epserr_iter``), because Er is written from the iteration's last
compute_Er call while <E_par B>/<B^2> is written after one more update.
Levels 0-2 are closed-form reductions of the same
force balance; ``poloidal_rotation_coefficient_from_neo2`` links Level 2 to
Level 3 through NEO-2's own transport coefficients.

Conventions (identical to ``compute_Er``; CGS-Gaussian units throughout)
-------------------------------------------------------------------------
- ``s`` is the normalized toroidal flux (``boozer_s``); ``r`` is the effective
  radius with ``d/dr = <|grad s|> d/ds`` and ``<|grad s|> = av_nabla_stor``
  [1/cm]. Radial derivatives are passed per unit ``s`` (``dn_spec_ov_ds``
  [1/cm^3], ``dT_spec_ov_ds`` [erg]), as NEO-2 stores them.
- ``E_r = -dPhi/dr`` [statV/cm]; multiply by ``STATV_PER_CM_TO_V_PER_M`` for
  V/m.
- ``psi_pr = sqrtg_bctrvr_tht`` [G cm] is sqrt(g) B^theta in symmetry flux
  coordinates with respect to ``r``, i.e. d psi_pol / dr with psi_pol the
  poloidal flux over 2 pi, in NEO-2's coordinate orientation.
- ``Om_tE = c E_r / psi_pr`` [rad/s] (``ntv_mod.f90``: ``c*Er/(aiota*
  sqrtg_bctrvr_phi)`` with ``sqrtg_bctrvr_phi = sqrtg_bctrvr_tht/aiota``).
- ``B_tht = bcovar_tht``, ``B_phi = bcovar_phi`` [G cm] are the covariant
  Boozer components; ``aiota`` is the rotational transform.
- ``vphi`` [rad/s] is the toroidal angular velocity ``Vphi`` of species ``i``
  in the ``isw_Vphi_loc = 0`` sense, i.e. the flux-surface average <V^phi>
  of the contravariant Boozer component (for u_i = 0 it equals the rigid
  geometric rotation frequency).
- Temperatures are energies [erg]; ``Z`` is the charge number (electrons -1).

Derivation
----------
Write the flow of species ``a`` on a flux surface as
``V_a = omega_a(psi) R^2 grad(phi) + u_a(psi) B``. Perpendicular momentum
balance (Hinton & Hazeltine 1976, Rev. Mod. Phys. 48, 239, Sec. VI; Helander
& Sigmar 2002, Collisional Transport in Magnetized Plasmas, Ch. 8) gives

    omega_a = Om_tE - c p_a' / (Z_a e n_a psi_pr),

with primes denoting d/dr. Equivalently, in laboratory components,
``E_r = p_a'/(Z_a e n_a) - v_theta B_tor / c + v_phi B_pol / c`` with
physical (laboratory) field components B_tor, B_pol [G], not to be confused
with the covariant Boozer components B_phi, B_tht [G cm] used elsewhere.
Using ``V.B = omega B_phi + u B^2`` and ``B^phi (B_phi + iota B_tht) = B^2``
(Boozer), the toroidal rotation that NEO-2 takes as input satisfies

    Vphi (B_phi + iota B_tht) = omega iota B_tht + <V_par B>.          (1)

Eliminating omega gives the exact single-species force balance

    E_r = p_i'/(Z_i e n_i) + psi_pr/(c iota B_tht) [Vphi (B_phi + iota B_tht)
          - <V_par,i B>].                                              (2)

Level 0 drops all flows, Level 1 sets u_i = 0 (so <V_par B> = omega B_phi and
omega = Vphi), Level 2 takes the neoclassical u_i from a poloidal rotation
coefficient k_i, and Level 3 inserts NEO-2's parallel flow
``<V_par,i B> = -sum_b [D31_ib A1_b + D32_ib A2_b + D33_ib A3_b]`` with
``A1_b = dln n_b/dr - 3/2 dln T_b/dr - Z_b e E_r / T_b``,
``A2_b = dln T_b/dr`` and ``A3_b = Z_b e <E_par B>/(T_b <B^2>)``
(Kernbichler et al. 2016, Plasma Phys. Control. Fusion 58, 104001). Because
A1 contains E_r, (2) becomes linear in E_r and the D31 sum enters the
denominator, exactly as in ``compute_Er``.
"""

import numpy as np

# Constants exactly as in NEO-2-QL/ntv_mod.f90 (rounded to 5 digits there), so
# that Level 3 replays compute_Er without constant-rounding differences.
C_CGS = 2.9979e10  # speed of light [cm/s]
E_CGS = 4.8032e-10  # elementary charge [statC]
STATV_PER_CM_TO_V_PER_M = 2.99792458e4  # 1 statV/cm = 29979.2458 V/m (exact)

# Asymptotic ion poloidal rotation coefficients k_i in the large-aspect-ratio
# limit (Hinton & Hazeltine 1976; Kim, Diamond & Groebner, Phys. Fluids B 3,
# 2050 (1991)). Sign convention: u_theta,i = k_i c/(Z_i e B) dT_i/dr, positive
# k in the banana regime; operationally k_i = 5/2 - D32_ii/D31_ii in NEO-2.
POLOIDAL_ROTATION_K_LIMITS = {
    'banana': 1.17,
    'plateau': -0.5,
    'pfirsch-schlueter': -2.1,
}


def diamagnetic_er(n, T, dn_ds, dT_ds, z, av_nabla_stor):
    """Pressure-gradient term p'/(Z e n) [statV/cm] of species ``(n, T, z)``."""
    n = np.asarray(n, dtype=float)
    T = np.asarray(T, dtype=float)
    dp_dr = (T * np.asarray(dn_ds) + n * np.asarray(dT_ds)) * av_nabla_stor
    return dp_dr / (z * E_CGS * n)


def omte_from_er(er, sqrtg_bctrvr_tht):
    """ExB rotation frequency Om_tE = c E_r / (sqrt(g) B^theta) [rad/s]."""
    return C_CGS * np.asarray(er) / np.asarray(sqrtg_bctrvr_tht)


def er_level0_diamagnetic(n, T, dn_ds, dT_ds, z, av_nabla_stor):
    """Level 0: E_r = p_i' / (Z_i e n_i), all flows of species i zero.

    This is (2) with Vphi = 0 and u_i = 0 (no toroidal and no poloidal
    rotation). Returns E_r [statV/cm].
    """
    return diamagnetic_er(n, T, dn_ds, dT_ds, z, av_nabla_stor)


def er_level1_toroidal_rotation(n, T, dn_ds, dT_ds, z, av_nabla_stor,
                                vphi, sqrtg_bctrvr_tht):
    """Level 1: add rigid toroidal rotation, poloidal flow neglected.

    With u_i = 0 in (1), omega_i = Vphi and (2) reduces to

        E_r = p_i'/(Z_i e n_i) + psi_pr Vphi / c,

    the familiar ``E_r = p'/(Z e n) + v_phi B_pol / c`` with
    ``v_phi = R Vphi`` and ``R B_pol = d psi_pol/dr`` (laboratory B_pol [G]). Returns E_r
    [statV/cm].
    """
    return (diamagnetic_er(n, T, dn_ds, dT_ds, z, av_nabla_stor)
            + np.asarray(sqrtg_bctrvr_tht) * np.asarray(vphi) / C_CGS)


def poloidal_rotation_er(k, T, dT_ds, z, av_nabla_stor, aiota,
                         bcovar_tht, bcovar_phi):
    """Poloidal-rotation term ``-v_theta B_tor / c`` of (2) [statV/cm].

    Neoclassical theory gives ``u_i <B^2> = k_i c B_phi T_i' / (Z_i e
    psi_pr)`` (Kim, Diamond & Groebner 1991; Hirshman & Sigmar 1981,
    Nucl. Fusion 21, 1079). Inserting <V_par B> = omega B_phi + u <B^2> into
    (1) and (2) gives

        -k_i B_phi T_i' / (Z_i e (B_phi + iota B_tht)),

    which tends to ``-k_i T_i'/(Z_i e)`` for B_phi >> iota B_tht.
    """
    dT_dr = np.asarray(dT_ds) * av_nabla_stor
    geom = bcovar_phi / (bcovar_phi + aiota * bcovar_tht)
    return -np.asarray(k) * geom * dT_dr / (z * E_CGS)


def er_level2_poloidal_rotation(n, T, dn_ds, dT_ds, z, av_nabla_stor,
                                vphi, sqrtg_bctrvr_tht, aiota,
                                bcovar_tht, bcovar_phi, k):
    """Level 2: Level 1 plus neoclassical poloidal rotation with coefficient k.

    E_r = p_i'/(Z_i e n_i) + psi_pr Vphi / c
          - k_i B_phi T_i' / (Z_i e (B_phi + iota B_tht)).

    ``k`` may be a regime value from ``POLOIDAL_ROTATION_K_LIMITS``,
    ``poloidal_rotation_coefficient_sauter`` or
    ``poloidal_rotation_coefficient_from_neo2``. With the latter, Level 2
    equals Level 3 exactly if (a) the ion row of NEO-2's coefficient matrix
    has no cross-species entries (D31_ib = D32_ib = 0 for b != i), (b) the
    diagonal entry satisfies the momentum identity of
    ``rigid_rotation_defect``, and (c) there is no inductive drive
    (sum_b D33_ib Z_b e <E_par B> / (T_b <B^2>) = 0, e.g. <E_par B> = 0).
    Cross-species D31 drives do not cancel even when
    the full row conserves momentum. Returns E_r [statV/cm].
    """
    return (er_level1_toroidal_rotation(n, T, dn_ds, dT_ds, z, av_nabla_stor,
                                        vphi, sqrtg_bctrvr_tht)
            + poloidal_rotation_er(k, T, dT_ds, z, av_nabla_stor, aiota,
                                   bcovar_tht, bcovar_phi))


def poloidal_rotation_coefficient_sauter(ftrap, nu_star_i):
    """Collisionality-dependent ion k from Sauter, Angioni & Lin-Liu (1999).

    Phys. Plasmas 6, 2834, ion poloidal-flow coefficient alpha, with the erratum Phys. Plasmas 9,
    5140 (2002); Sauter's alpha is the negative of k used here:

        alpha_0 = -1.17 (1 - f_t) / (1 - 0.22 f_t - 0.19 f_t^2)
        alpha = [(alpha_0 + 0.25 (1 - f_t^2) sqrt(nu)) / (1 + 0.5 sqrt(nu))
                 + 0.315 nu^2 f_t^6] / (1 + 0.15 nu^2 f_t^6)
        k = -alpha

    ``nu_star_i`` must be Sauter's ion collisionality, not NEO-2's
    ``nu_star_spec``.
    """
    ft = np.asarray(ftrap, dtype=float)
    nu = np.asarray(nu_star_i, dtype=float)
    alpha0 = -1.17 * (1.0 - ft) / (1.0 - 0.22 * ft - 0.19 * ft**2)
    alpha = ((alpha0 + 0.25 * (1.0 - ft**2) * np.sqrt(nu))
             / (1.0 + 0.5 * np.sqrt(nu))
             + 0.315 * nu**2 * ft**6) / (1.0 + 0.15 * nu**2 * ft**6)
    return -alpha


def _ion_row(spec_i, row, col):
    """Return (entry indices, column species) of the coefficient row spec_i."""
    row = np.asarray(row, dtype=int)
    col = np.asarray(col, dtype=int)
    entries = np.flatnonzero(row == spec_i)
    if entries.size == 0:
        raise ValueError('no transport coefficients for species %d' % spec_i)
    return entries, col[entries]


def er_level3_neo2_multispecies(spec_i, n, T, dn_ds, dT_ds, z, row, col,
                                D31, D32, vphi, aiota, av_nabla_stor,
                                sqrtg_bctrvr_tht, bcovar_tht, bcovar_phi,
                                D33=None, avEparB_ov_avb2=0.0,
                                vphi_loc_factor=None):
    """Level 3: full multi-species neoclassical E_r, replaying ``compute_Er``.

    Implements ``NEO-2-QL/ntv_mod.f90::compute_Er`` for ``isw_Vphi_loc = 0``.
    Inserting NEO-2's <V_par,i B> (module docstring) into (2) and collecting
    E_r gives

        E_r = N / D,
        D = c iota B_tht / psi_pr + sum_b D31_ib Z_b e / T_b,
        N = Vphi (iota B_tht + B_phi)
            + c iota B_tht T_i / (Z_i e psi_pr) dln p_i/dr
            + sum_b [D31_ib (dln n_b/dr + dln T_b/dr)
                     + (D32_ib - 5/2 D31_ib) dln T_b/dr
                     + D33_ib Z_b e <E_par B>/(<B^2> T_b)].

    Parameters
    ----------
    spec_i : int
        Zero-based index of the rotation species (``species_tag_Vphi``).
    n, T, dn_ds, dT_ds, z : array_like, shape (num_spec,)
        Species density [1/cm^3], temperature [erg], their s-derivatives and
        charge numbers.
    row, col : array_like of int
        Zero-based species indices of each coefficient entry.
    D31, D32, D33 : array_like
        Dimensional axisymmetric coefficients ``D31_AX``, ``D32_AX``,
        ``D33_AX`` (NEO-2 HDF5 output), same ordering as ``row``/``col``.
    vphi_loc_factor : float, optional
        ``None`` (default) for ``isw_Vphi_loc = 0``. For ``isw_Vphi_loc = 1,
        2`` (Vphi given at one point of the surface) the factor
        ``f = B^phi (1 + 2 pi iota psi_tor' G_symm_tb) / <B^2>`` at that point,
        see ``vphi_local_factor``. The local form of (1) is
        ``Vphi = omega (1 - f B_phi) + f <V_par B>``, and ``compute_Er`` then
        uses ``D = c (1 - f B_phi)/psi_pr + f sum_b D31_ib Z_b e/T_b`` and
        ``N = Vphi + c (1 - f B_phi) T_i/(Z_i e psi_pr) dln p_i/dr
        + f sum_b [...]``. ``f = 1/(B_phi + iota B_tht)`` (the flux-surface
        average <B^phi>/<B^2>) reproduces the ``isw_Vphi_loc = 0`` result.
    avEparB_ov_avb2 : float
        <E_par B>/<B^2> as written by NEO-2 (``avEparB_ov_avb2``), units
        statV/(cm G); only used with D33.

    Returns
    -------
    er : float
        E_r [statV/cm].
    terms : dict
        Additive decomposition: ``N = sum of 'num_*'``, ``D = sum of 'den_*'``.
    """
    n, T, dn_ds, dT_ds, z = (np.asarray(a, dtype=float)
                             for a in (n, T, dn_ds, dT_ds, z))
    entries, cols = _ion_row(spec_i, row, col)
    D31 = np.asarray(D31, dtype=float)[entries]
    D32 = np.asarray(D32, dtype=float)[entries]
    D33 = (np.zeros_like(D31) if D33 is None
           else np.asarray(D33, dtype=float)[entries])
    nb, Tb, zb = n[cols], T[cols], z[cols]
    dlnn_b = av_nabla_stor * dn_ds[cols] / nb
    dlnT_b = av_nabla_stor * dT_ds[cols] / Tb
    if vphi_loc_factor is None:
        # isw_Vphi_loc = 0 (ntv_mod.f90 compute_Er), multiplied through by
        # B_phi + iota B_tht.
        base = C_CGS * aiota * bcovar_tht / sqrtg_bctrvr_tht
        scale = 1.0
        num_vphi = vphi * (aiota * bcovar_tht + bcovar_phi)
    else:
        # isw_Vphi_loc = 1, 2 (ntv_mod.f90 compute_Er, fac1 branch).
        base = (C_CGS / sqrtg_bctrvr_tht) * (1.0 - vphi_loc_factor * bcovar_phi)
        scale = vphi_loc_factor
        num_vphi = vphi

    terms = {
        'den_base': base,
        'den_d31': scale * np.sum(D31 * zb * E_CGS / Tb),
        'num_vphi': num_vphi,
        'num_dia': base * diamagnetic_er(n[spec_i], T[spec_i], dn_ds[spec_i],
                                         dT_ds[spec_i], z[spec_i],
                                         av_nabla_stor),
        'num_d31': scale * np.sum(D31 * (dlnn_b + dlnT_b)),
        'num_d32': scale * np.sum((D32 - 2.5 * D31) * dlnT_b),
        'num_d33': scale * np.sum(D33 * avEparB_ov_avb2 * zb * E_CGS / Tb),
    }
    num = sum(v for key, v in terms.items() if key.startswith('num_'))
    den = terms['den_base'] + terms['den_d31']
    return num / den, terms


def vphi_local_factor(bctrvr_phi_Vphi, G_symm_tb_Vphi, aiota, psi_pr_hat,
                      Bref, avbhat2):
    """Factor f = B^phi (1 + 2 pi iota psi_tor' G_symm_tb) / <B^2>.

    Evaluated at the point where Vphi is given (``isw_Vphi_loc = 1, 2``),
    exactly as ``fac1`` in ``compute_Er``: ``bctrvr_phi_Vphi`` is the
    contravariant Boozer B^phi there, ``G_symm_tb_Vphi`` the theta derivative
    of G_symm (``compute_Gsymm``), which converts the Boozer toroidal angle
    to that of symmetry flux coordinates; psi_tor' = ``psi_pr_hat * Bref``
    and <B^2> = ``avbhat2 * Bref**2``.
    """
    return (bctrvr_phi_Vphi
            * (1.0 + 2.0 * np.pi * aiota * psi_pr_hat * Bref * G_symm_tb_Vphi)
            / (avbhat2 * Bref**2))


def poloidal_rotation_coefficient_from_neo2(spec_i, row, col, D31, D32):
    """Effective k_i = 5/2 - D32_ii / D31_ii from NEO-2's diagonal ion block.

    For a single ion species whose coefficient row has no cross-species
    entries, satisfies the momentum identity of ``rigid_rotation_defect`` and
    has no inductive drive, (2) with NEO-2's parallel flow reproduces Level 2
    with exactly this k (Kim, Diamond & Groebner 1991 write the same relation
    in terms of viscosity coefficients). Otherwise it is an effective value.
    """
    entries, cols = _ion_row(spec_i, row, col)
    diag = entries[cols == spec_i]
    return 2.5 - np.asarray(D32)[diag][0] / np.asarray(D31)[diag][0]


def rigid_rotation_defect(spec_i, T, z, row, col, D31, sqrtg_bctrvr_tht,
                          bcovar_phi):
    """Relative defect of ``sum_b D31_ib Z_b e psi_pr / (c T_b) = B_phi``.

    A rigid toroidal rotation of all species with common omega is an exact
    solution of the drift-kinetic equation with momentum-conserving
    collisions, carrying <V_par B> = omega B_phi. In NEO-2's flux-force
    form this requires the identity above; with it the Level 3 denominator
    equals c (B_phi + iota B_tht) / psi_pr, the "D31 denominator
    correction" ``1 + B_phi/(iota B_tht)`` of the issue #75 audit. Returns
    ``(lhs - B_phi) / B_phi``.
    """
    T = np.asarray(T, dtype=float)
    z = np.asarray(z, dtype=float)
    entries, cols = _ion_row(spec_i, row, col)
    lhs = (np.sum(np.asarray(D31)[entries] * z[cols] * E_CGS / T[cols])
           * sqrtg_bctrvr_tht / C_CGS)
    return (lhs - bcovar_phi) / bcovar_phi


REQUIRED_DATASETS = (
    'species_tag', 'species_tag_Vphi', 'isw_Vphi_loc', 'n_spec', 'T_spec',
    'dn_spec_ov_ds', 'dT_spec_ov_ds', 'z_spec', 'row_ind_spec',
    'col_ind_spec', 'D31_AX', 'D32_AX', 'D33_AX', 'Vphi', 'aiota',
    'av_nabla_stor', 'sqrtg_bctrvr_tht', 'bcovar_tht', 'bcovar_phi', 'Er')

LOCAL_VPHI_DATASETS = ('bctrvr_phi_Vphi', 'G_symm_tb_Vphi', 'psi_pr_hat',
                       'Bref', 'avbhat2')


def load_neo2_force_balance_inputs(path):
    """Read the inputs of ``er_level3_neo2_multispecies`` from NEO-2 output.

    ``path`` is a single-surface ``neo2_multispecies_out.h5`` written with
    ``isw_calc_Er = 1``. Returns a dict of keyword arguments for
    ``er_level3_neo2_multispecies`` plus the stored ``Er``. Species tags in
    ``row_ind_spec``/``col_ind_spec`` are mapped to zero-based indices.

    Besides geometry, coefficients and ``Er``, the replay needs
    ``dn_spec_ov_ds``, ``dT_spec_ov_ds``, ``Vphi``, ``species_tag_Vphi`` and
    ``isw_Vphi_loc``, which ``write_multispec_output_a`` writes only in
    revisions containing the issue #75 output change; older files raise
    ``KeyError``. For more than one species ``avEparB_ov_avb2`` is required
    as well. For ``isw_Vphi_loc = 1, 2`` (Vphi given at one point) the
    local B^phi and G_symm_tb at that point (``bctrvr_phi_Vphi``,
    ``G_symm_tb_Vphi``) and ``psi_pr_hat``, ``Bref``, ``avbhat2`` are
    needed to form ``vphi_loc_factor``. Mode 2 is checked against NEO-2
    runs; mode 1 shares the same algebra and differs only in how theta_B is
    found from (R_Vphi, Z_Vphi). In the tested revision, NEO-2 itself STOPs
    in mode 1 for the tested points (``calc_thetaB_RZloc`` runs separate Newton
    iterations on R and Z, which converge to different roots).
    """
    import h5py

    with h5py.File(path, 'r') as f:
        g = {key: f[key][()] for key in f.keys()
             if isinstance(f[key], h5py.Dataset)}
    required = list(REQUIRED_DATASETS)
    if len(np.atleast_1d(g.get('species_tag', []))) > 1:
        # Written only for num_spec > 1, where it enters through D33.
        required.append('avEparB_ov_avb2')
    missing = [key for key in required if key not in g]
    mode = int(g['isw_Vphi_loc']) if 'isw_Vphi_loc' in g else 0
    if mode not in (0, 1, 2):
        raise ValueError('isw_Vphi_loc must be 0, 1 or 2')
    if mode >= 1:
        missing += [key for key in LOCAL_VPHI_DATASETS if key not in g]
    if missing:
        raise KeyError('NEO-2 output lacks datasets needed for the E_r '
                       'replay: ' + ', '.join(missing))
    vphi_loc_factor = None
    if mode >= 1:
        vphi_loc_factor = float(vphi_local_factor(
            g['bctrvr_phi_Vphi'], g['G_symm_tb_Vphi'], g['aiota'],
            g['psi_pr_hat'], g['Bref'], g['avbhat2']))
    tags = np.asarray(g['species_tag'], dtype=int)
    index = {int(t): i for i, t in enumerate(tags)}
    matches = np.flatnonzero(tags == int(g['species_tag_Vphi']))
    if matches.size != 1:
        raise ValueError('species_tag_Vphi must match exactly one species')
    return {
        'spec_i': int(matches[0]),
        'n': g['n_spec'], 'T': g['T_spec'],
        'dn_ds': g['dn_spec_ov_ds'], 'dT_ds': g['dT_spec_ov_ds'],
        'z': g['z_spec'],
        'row': np.array([index[int(t)] for t in g['row_ind_spec']]),
        'col': np.array([index[int(t)] for t in g['col_ind_spec']]),
        'D31': g['D31_AX'], 'D32': g['D32_AX'], 'D33': g['D33_AX'],
        'avEparB_ov_avb2': float(g.get('avEparB_ov_avb2', 0.0)),
        'vphi': float(g['Vphi']), 'aiota': float(g['aiota']),
        'av_nabla_stor': float(g['av_nabla_stor']),
        'sqrtg_bctrvr_tht': float(g['sqrtg_bctrvr_tht']),
        'bcovar_tht': float(g['bcovar_tht']),
        'bcovar_phi': float(g['bcovar_phi']),
        'vphi_loc_factor': vphi_loc_factor,
        'Er_stored': float(g['Er']),
    }
