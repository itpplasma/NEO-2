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
  Level 1: adds toroidal rotation through a coordinate-consistent pair
  Level 2: adds a neoclassical poloidal rotation estimate
  Reduced single-ion NEO-2 limit: analytic D31 reduction of the full solver
  Reduced NEO-2 Vphi convention: isw_Vphi_loc=0 algebra without transport
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


def compute_boozer_metric(m_modes, n_modes, rmnc, rmns, zmnc, zmns,
                          nper, theta_B, phi=0.0):
    """Compute R and |e_theta| from Boozer Fourier harmonics.

    Evaluates position and poloidal basis vector magnitude at a given
    (theta_B, phi) point by analytically differentiating the Fourier
    series from a Boozer coordinate file.

    Parameters
    ----------
    m_modes, n_modes : array_like
        Poloidal and toroidal mode numbers.
    rmnc, rmns, zmnc, zmns : array_like
        Fourier harmonics for R and Z [m].
    nper : int
        Number of field periods.
    theta_B : float
        Boozer poloidal angle [rad].
    phi : float, optional
        Boozer toroidal angle [rad]. Default 0 (outboard midplane).

    Returns
    -------
    R_cm : float
        Major radius [cm].
    abs_e_theta_cm : float
        Poloidal basis vector magnitude |e_theta| [cm].
    """
    m_modes = np.asarray(m_modes, dtype=float)
    n_modes = np.asarray(n_modes, dtype=float)
    rmnc = np.asarray(rmnc, dtype=float)
    rmns = np.asarray(rmns, dtype=float)
    zmnc = np.asarray(zmnc, dtype=float)
    zmns = np.asarray(zmns, dtype=float)

    angle = m_modes * theta_B - nper * n_modes * phi
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)

    R_m = np.dot(rmnc, cos_a) + np.dot(rmns, sin_a)
    Z_m = np.dot(zmnc, cos_a) + np.dot(zmns, sin_a)

    dR_dth = np.dot(-m_modes * rmnc, sin_a) + np.dot(m_modes * rmns, cos_a)
    dZ_dth = np.dot(-m_modes * zmnc, sin_a) + np.dot(m_modes * zmns, cos_a)

    R_cm = R_m * 100.0
    abs_e_theta_cm = np.sqrt(dR_dth**2 + dZ_dth**2) * 100.0
    return R_cm, abs_e_theta_cm


def compute_boozer_metric_from_bc(bc_file_path, surface_s, theta_B=0.0,
                                  phi=0.0):
    """Load a .bc file and return R, |e_theta| at (s, theta_B, phi).

    Parameters
    ----------
    bc_file_path : str or path-like
        Path to the Boozer .bc file.
    surface_s : float or array_like
        Normalized toroidal flux label(s).
    theta_B : float, optional
        Boozer poloidal angle [rad]. Default 0.
    phi : float, optional
        Boozer toroidal angle [rad]. Default 0.

    Returns
    -------
    R_cm : float or ndarray
        Major radius [cm].
    abs_e_theta_cm : float or ndarray
        Poloidal basis vector magnitude [cm].
    """
    from libneo import BoozerFile

    bc = BoozerFile(str(bc_file_path))
    s_arr = np.asarray(bc.s, dtype=float)
    surface_s = np.atleast_1d(np.asarray(surface_s, dtype=float))
    R_list, e_list = [], []
    for s_val in surface_s:
        idx = int(np.argmin(np.abs(s_arr - s_val)))
        R_cm, e_cm = compute_boozer_metric(
            bc.m[idx], bc.n[idx], bc.rmnc[idx], bc.rmns[idx],
            bc.zmnc[idx], bc.zmns[idx], bc.nper, theta_B, phi,
        )
        R_list.append(R_cm)
        e_list.append(e_cm)
    if len(R_list) == 1:
        return R_list[0], e_list[0]
    return np.array(R_list), np.array(e_list)


def compute_boozer_metric_from_rz_profile(R_prof_cm, Z_prof_cm):
    """Compute R and |e_theta| at theta_B=0 from NEO-2 profile arrays.

    Parameters
    ----------
    R_prof_cm : array_like
        R values [cm] at uniformly-spaced Boozer theta from 0 to 2*pi
        (101 points, as stored in R_Vphi_prof).
    Z_prof_cm : array_like
        Z values [cm] at the same theta grid.

    Returns
    -------
    R_cm : float
        Major radius at theta_B=0 [cm].
    abs_e_theta_cm : float
        |e_theta| at theta_B=0 [cm], via central finite differences.
    """
    R_prof = np.asarray(R_prof_cm, dtype=float)
    Z_prof = np.asarray(Z_prof_cm, dtype=float)
    npts = R_prof.size
    dtheta = 2.0 * np.pi / (npts - 1)
    dR = (R_prof[1] - R_prof[-2]) / (2.0 * dtheta)
    dZ = (Z_prof[1] - Z_prof[-2]) / (2.0 * dtheta)
    return float(R_prof[0]), float(np.sqrt(dR**2 + dZ**2))


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
    v_theta=None,
    b_phi=None,
):
    """Compute Om_tE from radial force balance.

    Uses the radial force balance:
        E_r = (1 / (Z e n)) dp/dr
            + v_phi * B_theta / c
            - v_theta * B_phi / c
    and converts to the E x B rotation frequency:
        Om_tE = c E_r / (iota * sqrt(g) B^phi)

    Only the products ``v_phi * B_theta`` and ``v_theta * B_phi`` enter.
    These must be supplied as physical cylindrical pairs: velocities [cm/s]
    with fields [G].

    .. warning::

       Passing raw NEO-2 Boozer pairs ``(V^phi, bcovar_tht)`` here gives
       the wrong toroidal product because ``B_theta_cov / R != B_pol``;
       the correct relation is ``B_pol = B_theta_cov / |e_theta|`` where
       ``|e_theta|`` is the poloidal basis vector magnitude (approximately
       the minor radius, NOT the major radius R). Use
       ``compute_omte_toroidal_rotation_physical`` or
       ``compute_omte_neoclassical_poloidal_physical`` when the inputs are
       Boozer-coordinate quantities from NEO-2.

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
        Toroidal rotation quantity paired with `b_theta`. Use either a
        physical toroidal velocity [cm/s] or the NEO-2 contravariant
        toroidal angular frequency `Vphi` [rad/s].
    b_theta : float or array, optional
        Poloidal field quantity paired with `v_phi`. Use either a physical
        field [G] or the NEO-2 Boozer covariant component `bcovar_tht`
        [G cm]. Must be provided together with `v_phi` for Level 1.
    v_theta : float or array, optional
        Poloidal rotation quantity paired with `b_phi`. Use either a
        physical poloidal velocity [cm/s] or a Boozer contravariant
        poloidal angular frequency [rad/s].
    b_phi : float or array, optional
        Toroidal field quantity paired with `v_theta`. Use either a
        physical field [G] or the NEO-2 Boozer covariant component
        `bcovar_phi` [G cm]. Must be provided together with `v_theta`
        for Level 2.

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

    ``compute_omte_from_neo2_output`` provides exact replay paths from stored
    NEO-2 output data. Those are not part of the reduced Level 0-2 hierarchy.
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
    if (v_theta is None) != (b_phi is None):
        raise ValueError('v_theta and b_phi must be provided together')

    Er = _compute_diamagnetic_er(n, T, dn_ds, dT_ds, z, av_nabla_stor)

    if v_phi is not None:
        Er = Er + np.asarray(v_phi) * np.asarray(b_theta) / C_CGS
    if v_theta is not None:
        Er = Er - np.asarray(v_theta) * np.asarray(b_phi) / C_CGS

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
    """Compute Om_tE including toroidal rotation (Level 1).

    ``v_phi`` and ``b_theta`` must be a physical toroidal-velocity [cm/s]
    and poloidal-field [G] pair.  For Boozer-coordinate inputs from NEO-2,
    use ``compute_omte_toroidal_rotation_physical`` which applies the
    metric correction factor ``R / |e_theta|``.
    """
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


def compute_omte_toroidal_rotation_neo2_convention(
    n,
    T,
    dn_ds,
    dT_ds,
    z,
    aiota,
    sqrtg_bctrvr_phi,
    av_nabla_stor,
    vphi,
    bcovar_tht,
    bcovar_phi,
):
    """Compute the reduced isw_Vphi_loc=0 numerator/denominator without transport.

    This mirrors the algebraic structure of the Fortran ``compute_Er()``
    subroutine for ``isw_Vphi_loc=0``, but **without** the D31/D32/D33
    transport coefficient sums.  It is NOT a standalone physical model.

    In the full Fortran solver, ``Vphi`` [rad/s] (the species toroidal angular
    frequency) enters through:

        nom = Vphi * (iota*B_theta_cov + B_phi_cov)
              + c*T*B_theta_cov/(Z*e*sqrtg) * dp/dr/p
              + sum over transport terms (D31, D32, D33)
        denom = c*B_theta_cov/sqrtg + sum over D31 terms
        Er = nom / denom

    Without the transport sums, the Vphi term dominates and gives results
    that are orders of magnitude off from the self-consistent Er.  Use this
    function only as a building block inside the full-output replay path
    (``neo2_output_omte.py``), not as a standalone Om_tE estimate.

    Parameters
    ----------
    vphi : float or array
        Species toroidal angular frequency [rad/s] (from NEO-2 HDF5 ``Vphi``).
    bcovar_tht : float or array
        Boozer covariant poloidal B component [G cm].
    bcovar_phi : float or array
        Boozer covariant toroidal B component [G cm].
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
    bcovar_tht = np.asarray(bcovar_tht)
    bcovar_phi = np.asarray(bcovar_phi)
    vphi = np.asarray(vphi)

    if np.any(bcovar_tht == 0.0):
        raise ValueError('bcovar_tht must be nonzero to compute NEO-2 Vphi model')

    dp_ds = T * dn_ds + n * dT_ds
    dp_dr = dp_ds * av_nabla_stor
    pressure = n * T
    if np.any(pressure == 0.0):
        raise ValueError('n*T must be nonzero to compute NEO-2 Vphi model')
    denom_er = C_CGS * bcovar_tht / sqrtg_bctrvr_phi
    nom_er = (
        vphi * (aiota * bcovar_tht + bcovar_phi)
        + (C_CGS * T * bcovar_tht / (z * E_CGS * sqrtg_bctrvr_phi))
        * (dp_dr / pressure)
    )
    er = nom_er / denom_er
    om_tE = C_CGS * er / (aiota * sqrtg_bctrvr_phi)
    return om_tE, er


def select_poloidal_rotation_coefficient(nu_star):
    """Select K_i from a simple collisionality regime map."""
    nu_star = np.asarray(nu_star)
    k_i = np.empty_like(nu_star, dtype=float)
    k_i[nu_star < 0.1] = -1.17
    k_i[(nu_star >= 0.1) & (nu_star < 10.0)] = -0.5
    k_i[nu_star >= 10.0] = 0.5
    return k_i


def compute_poloidal_rotation_neoclassical(dT_ds, z, b_phi, av_nabla_stor, k_i):
    """Estimate the Level 2 poloidal rotation quantity from dT/dr.

    Uses the Kim-Diamond-Groebner estimate (KDG 1991, eq. 6) in CGS:

        v_theta = K_i * c / (Z e B_phi) * dT/dr

    with dT/dr = dT/ds * <|nabla s|>.

    The factor of c (speed of light) is required in CGS and cancels
    against the 1/c in the force balance term -v_theta * B_phi / c,
    giving the net contribution Er_pol = -K_i * dT/dr / (Z e).

    If `b_phi` is a physical toroidal field [G], the return value is a
    physical poloidal velocity [cm/s]. If `b_phi` is the NEO-2 Boozer
    covariant component `bcovar_phi` [G cm], the return value is the
    matching contravariant poloidal angular frequency [rad/s]. In both cases
    the product `v_theta * b_phi` entering force balance is identical
    (pair-product identity).
    """
    dT_ds = np.asarray(dT_ds)
    b_phi = np.asarray(b_phi)
    av_nabla_stor = np.asarray(av_nabla_stor)

    if z == 0.0:
        raise ValueError('z must be nonzero to compute v_theta')
    if np.any(b_phi == 0.0):
        raise ValueError('b_phi must be nonzero to compute v_theta')

    dT_dr = dT_ds * av_nabla_stor
    return k_i * C_CGS * dT_dr / (z * E_CGS * b_phi)


def compute_omte_neoclassical_poloidal(
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
    b_phi,
    k_i,
):
    """Compute Om_tE with toroidal and estimated poloidal rotation (Level 2).

    ``v_phi`` and ``b_theta`` must be physical quantities [cm/s] and [G].
    The KDG poloidal term is coordinate-independent (B_phi cancels), but
    the toroidal product is not.  For Boozer-coordinate inputs from NEO-2,
    use ``compute_omte_neoclassical_poloidal_physical``.
    """
    v_theta = compute_poloidal_rotation_neoclassical(
        dT_ds=dT_ds,
        z=z,
        b_phi=b_phi,
        av_nabla_stor=av_nabla_stor,
        k_i=k_i,
    )
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
        v_theta=v_theta,
        b_phi=b_phi,
    )


def compute_omte_neo2_single_ion_limit(
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
    b_phi,
    k_i,
):
    """Compute the reduced single-ion NEO-2 force-balance limit.

    This implements the analytic reduction documented in
    ``ripple_solver_normalizations_and_output.tex``:

        E_r =
            sqrt(g) B^theta / c * Vphi
          + (1 / Z e n) dp/dr
          - k_i * B_phi / (Z e (iota B_theta + B_phi)) * dT/dr

    It is the exact single-ion, negligible-A3 limit of the NEO-2
    ``compute_Er()`` algebra after substituting the analytic D31 bootstrap
    coefficient. Unlike the local cylindrical shortcut used in Level 2, this
    expression keeps the Boozer geometry factor that absorbs the large D31
    denominator into the reduced force balance.
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

    b_theta = np.asarray(b_theta)
    b_phi = np.asarray(b_phi)
    v_phi = np.asarray(v_phi)
    denom_geom = aiota * b_theta + b_phi
    if np.any(denom_geom == 0.0):
        raise ValueError('aiota * b_theta + b_phi must be nonzero')

    dp_ds = T * dn_ds + n * dT_ds
    dp_dr = dp_ds * av_nabla_stor
    dT_dr = dT_ds * av_nabla_stor
    sqrtg_bctrvr_tht = aiota * sqrtg_bctrvr_phi

    Er = (
        sqrtg_bctrvr_tht * v_phi / C_CGS
        + dp_dr / (n * z * E_CGS)
        - (k_i * b_phi / (z * E_CGS * denom_geom)) * dT_dr
    )
    Om_tE = C_CGS * Er / (aiota * sqrtg_bctrvr_phi)
    return Om_tE, Er


def compute_omte_neoclassical_poloidal_auto_k(
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
    b_phi,
    nu_star,
):
    """Compute Level 2 using a simple collisionality-based K_i selection."""
    return compute_omte_neoclassical_poloidal(
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
        b_phi=b_phi,
        k_i=select_poloidal_rotation_coefficient(nu_star),
    )


def compute_omte_toroidal_rotation_boozer(
    n,
    T,
    dn_ds,
    dT_ds,
    z,
    aiota,
    sqrtg_bctrvr_phi,
    av_nabla_stor,
    vphi,
):
    """Compute Om_tE including toroidal rotation in Boozer coordinates (Level 1).

    Uses the correct Boozer-algebra rotation product
    ``sqrt(g) B^theta * V^phi / c`` instead of the incorrect covariant
    pairing ``V^phi * B_theta_cov / c``.  The covariant ``B_theta_cov``
    and the contravariant ``sqrt(g) B^theta`` differ in sign and magnitude
    because the inverse metric mixes the dominant ``B_phi_cov`` component
    into ``B^theta``.
    """
    (
        n, T, dn_ds, dT_ds, aiota, sqrtg_bctrvr_phi, av_nabla_stor,
    ) = _prepare_common_inputs(
        n, T, dn_ds, dT_ds, z, aiota, sqrtg_bctrvr_phi, av_nabla_stor
    )

    sqrtg_bctrvr_tht = aiota * sqrtg_bctrvr_phi

    Er = _compute_diamagnetic_er(n, T, dn_ds, dT_ds, z, av_nabla_stor)
    Er = Er + sqrtg_bctrvr_tht * np.asarray(vphi) / C_CGS

    Om_tE = C_CGS * Er / (aiota * sqrtg_bctrvr_phi)
    return Om_tE, Er


def compute_omte_neoclassical_poloidal_boozer(
    n,
    T,
    dn_ds,
    dT_ds,
    z,
    aiota,
    sqrtg_bctrvr_phi,
    av_nabla_stor,
    vphi,
    k_i,
):
    """Compute Om_tE with Boozer-correct rotation and KDG poloidal estimate (Level 2).

    Uses ``sqrt(g) B^theta * V^phi / c`` for the toroidal rotation term
    and the coordinate-independent KDG poloidal contribution
    ``-k_i * dT/dr / (Z e)``.  This is the minimal correct Boozer-algebra
    force balance that agrees with the reduced single-ion NEO-2 limit to
    within the small geometry factor ``B_phi / (iota B_theta + B_phi)``
    on the KDG term (typically ~2% correction).
    """
    (
        n, T, dn_ds, dT_ds, aiota, sqrtg_bctrvr_phi, av_nabla_stor,
    ) = _prepare_common_inputs(
        n, T, dn_ds, dT_ds, z, aiota, sqrtg_bctrvr_phi, av_nabla_stor
    )

    sqrtg_bctrvr_tht = aiota * sqrtg_bctrvr_phi
    dT_dr = dT_ds * av_nabla_stor

    Er = _compute_diamagnetic_er(n, T, dn_ds, dT_ds, z, av_nabla_stor)
    Er = Er + sqrtg_bctrvr_tht * np.asarray(vphi) / C_CGS
    Er = Er - np.asarray(k_i) * dT_dr / (z * E_CGS)

    Om_tE = C_CGS * Er / (aiota * sqrtg_bctrvr_phi)
    return Om_tE, Er
