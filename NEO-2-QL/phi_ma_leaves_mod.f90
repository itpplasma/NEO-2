module phi_ma_leaves_mod
    ! Generated pointwise leaves; kept free of NEO-2 dependencies so that
    ! ntv_mod can use them as well as the Phi^MA response.
    implicit none
    private
    public :: phi_ma_scales_leaf, thermodynamic_forces_leaf, field_line_theta_leaf, &
        phi_ma_phase_leaf, phi_ma_drift_leaf, phi_ma_band_leaf, phi_ma_jhat_leaf

contains

    include 'generated/phi_ma_leaves.inc'

end module phi_ma_leaves_mod
