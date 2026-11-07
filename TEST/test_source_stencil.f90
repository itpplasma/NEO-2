!> Stencils of the QL source and flux quadrature (source_flux, add_f01_source)
!> must stay on the field line [ibeg, iend] for either parity of iend - ibeg
!> and integrate a quadratic exactly; for even iend - ibeg they must keep the
!> original parity rule.
program test_source_stencil
    use source_stencil_mod, only: forward_stencil_centred, backward_stencil_centred
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: w_near = 5.0d0/12.0d0, w_mid = 8.0d0/12.0d0, &
        w_far = -1.0d0/12.0d0
    integer :: status, ibeg, iend

    status = 0
    do ibeg = -2, 0
        do iend = ibeg + 2, ibeg + 9
            call check_line(ibeg, iend, status)
        end do
    end do

    if (status == 0) then
        print *, "All tests passed!"
    else
        print *, "FAIL:", status, "checks failed"
        error stop 1
    end if

contains

    subroutine check_line(ibeg, iend, status)
        integer, intent(in) :: ibeg, iend
        integer, intent(inout) :: status
        integer :: istep, lo, hi
        real(dp) :: int_fwd, int_bwd, exact
        logical :: centred

        int_fwd = 0.0d0
        int_bwd = 0.0d0
        do istep = ibeg, iend
            if (istep > ibeg) then
                centred = forward_stencil_centred(istep, ibeg, iend)
                if (centred) then
                    lo = istep - 1; hi = istep + 1
                    int_fwd = int_fwd + w_near*f(istep - 1) + w_mid*f(istep) + w_far*f(istep + 1)
                else
                    lo = istep - 2; hi = istep
                    int_fwd = int_fwd + w_far*f(istep - 2) + w_mid*f(istep - 1) + w_near*f(istep)
                end if
                call check_stencil('forward', ibeg, iend, istep, lo, hi, centred, status)
            end if
            if (istep < iend) then
                centred = backward_stencil_centred(istep, ibeg, iend)
                if (centred) then
                    lo = istep - 1; hi = istep + 1
                    int_bwd = int_bwd + w_far*f(istep - 1) + w_mid*f(istep) + w_near*f(istep + 1)
                else
                    lo = istep; hi = istep + 2
                    int_bwd = int_bwd + w_near*f(istep) + w_mid*f(istep + 1) + w_far*f(istep + 2)
                end if
                call check_stencil('backward', ibeg, iend, istep, lo, hi, centred, status)
            end if
        end do

        exact = (real(iend, dp)**3 - real(ibeg, dp)**3)/3.0d0
        if (abs(int_fwd - exact) > 1.0d-12*max(1.0d0, abs(exact)) .or. &
            abs(int_bwd - exact) > 1.0d-12*max(1.0d0, abs(exact))) then
            print *, "FAIL quadrature", ibeg, iend, int_fwd, int_bwd, exact
            status = status + 1
        end if
    end subroutine check_line

    subroutine check_stencil(sweep, ibeg, iend, istep, lo, hi, centred, status)
        character(*), intent(in) :: sweep
        integer, intent(in) :: ibeg, iend, istep, lo, hi
        logical, intent(in) :: centred
        integer, intent(inout) :: status

        if (lo < ibeg .or. hi > iend) then
            print *, "FAIL ", sweep, " stencil leaves line: ibeg, iend, istep, lo, hi =", &
                ibeg, iend, istep, lo, hi
            status = status + 1
        end if
        if (mod(iend - ibeg, 2) == 0 .and. (centred .neqv. (mod(istep - ibeg, 2) == 1))) then
            print *, "FAIL ", sweep, " parity rule changed: ibeg, iend, istep =", &
                ibeg, iend, istep
            status = status + 1
        end if
    end subroutine check_stencil

    pure real(dp) function f(x)
        integer, intent(in) :: x

        f = real(x, dp)**2
    end function f

end program test_source_stencil
