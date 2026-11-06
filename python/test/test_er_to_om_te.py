import numpy as np
import pytest

from neo2_ql.er_to_om_te import er_to_om_te


def test_matches_gaussian_unit_conversion():
    # Om_tE = c Er / sqrtg_bctrvr_tht in Gaussian units (ntv_mod compute_Er),
    # with Er in statV/cm = (V/m) * 1e-2 / 299.792458.
    c_cgs = 2.99792458e10
    er_si = np.array([-3.0e3, 0.0, 1.2e4])
    geom = np.array([2.5e5, -4.0e5, 1.0e6])
    expected = c_cgs * (er_si * 1.0e-2 / 299.792458) / geom
    assert np.allclose(er_to_om_te(er_si, geom), expected, rtol=1e-14)


def test_rejects_zero_geometry():
    with pytest.raises(ValueError):
        er_to_om_te(1.0, 0.0)
