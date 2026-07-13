module phi_crossing_alignment_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

    pure integer function nearest_phi_crossing_offset(bhat, eta_cross) &
            result(offset)
        real(dp), intent(in) :: bhat(-1:1)
        real(dp), intent(in) :: eta_cross

        real(dp) :: residual(-1:1)

        residual = abs(bhat*eta_cross - 1.0_dp)
        offset = 0
        if (residual(-1) < residual(offset)) offset = -1
        if (residual(1) < residual(offset)) offset = 1
    end function nearest_phi_crossing_offset

end module phi_crossing_alignment_mod
