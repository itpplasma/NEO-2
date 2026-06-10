! Adaptive 1d quadrature built on the libneo math kit (neo_gauss_kronrod),
! replacing the former FGSL/GSL wrappers (issue #88). fint1d_cquad keeps the
! interface of the old CQUAD wrapper but integrates with adaptive
! Gauss-Kronrod 21. fint1d_qagiu maps [x_low, infinity) onto (0, 1] via
! x = x_low + (1 - t)/t, the same transformation QUADPACK QAGI uses, and
! integrates the transformed integrand with the 15-point rule.
module integration_routines_mod

    use nrtype, only : dp
    use neo_gauss_kronrod, only : integrate_gk

    implicit none

    public fint1d_cquad
    private fint1d_param0_cquad, fint1d_param1_cquad
    interface fint1d_cquad
        module procedure fint1d_param0_cquad, fint1d_param1_cquad
    end interface fint1d_cquad

    public fint1d_qagiu
    private fint1d_param0_qagiu, fint1d_param1_qagiu
    interface fint1d_qagiu
        module procedure fint1d_param0_qagiu, fint1d_param1_qagiu
    end interface fint1d_qagiu

    ! Accumulated QUADPACK error codes (0 = success, 1 = subdivision limit,
    ! 2 = roundoff, 3 = bad integrand behaviour, 6 = invalid input),
    ! reported and reset by disp_integration_error, as the old
    ! disp_gsl_integration_error did for GSL codes.
    public disp_integration_error
    private record_error
    integer, parameter, private :: max_err_code = 6
    logical, dimension(0:max_err_code), private :: err_detected = .false.

    integer, parameter, private :: gk_limit = 1000

contains

    recursive function fint1d_param0_cquad(func1d_param0_user, x_low, x_up, &
            epsabs, epsrel) result(res)
        interface
            function func1d_param0_user(x)
                use nrtype, only : dp
                real(dp) :: func1d_param0_user
                real(dp) :: x
            end function func1d_param0_user
        end interface
        real(dp) :: x_low, x_up
        real(dp) :: epsabs, epsrel
        real(dp), dimension(2) :: res

        integer :: ierr

        call integrate_gk(eval, x_low, x_up, epsabs, epsrel, res(1), res(2), &
                          ierr, key=21, limit=gk_limit)
        call record_error(ierr)

    contains

        recursive function eval(x) result(fx)
            real(dp), intent(in) :: x
            real(dp) :: fx
            real(dp) :: x_loc

            x_loc = x
            fx = func1d_param0_user(x_loc)
        end function eval

    end function fint1d_param0_cquad


    recursive function fint1d_param1_cquad(func1d_param1_user, param1, x_low, &
            x_up, epsabs, epsrel) result(res)
        interface
            function func1d_param1_user(x, param1)
                use nrtype, only : dp
                real(dp) :: func1d_param1_user
                real(dp) :: x, param1
            end function func1d_param1_user
        end interface
        real(dp) :: param1
        real(dp) :: x_low, x_up
        real(dp) :: epsabs, epsrel
        real(dp), dimension(2) :: res

        integer :: ierr

        call integrate_gk(eval, x_low, x_up, epsabs, epsrel, res(1), res(2), &
                          ierr, key=21, limit=gk_limit)
        call record_error(ierr)

    contains

        recursive function eval(x) result(fx)
            real(dp), intent(in) :: x
            real(dp) :: fx
            real(dp) :: x_loc

            x_loc = x
            fx = func1d_param1_user(x_loc, param1)
        end function eval

    end function fint1d_param1_cquad


    recursive function fint1d_param0_qagiu(func1d_param0_user, x_low, epsabs, &
            epsrel) result(res)
        interface
            function func1d_param0_user(x)
                use nrtype, only : dp
                real(dp) :: func1d_param0_user
                real(dp) :: x
            end function func1d_param0_user
        end interface
        real(dp) :: x_low
        real(dp) :: epsabs, epsrel
        real(dp), dimension(2) :: res

        integer :: ierr

        call integrate_gk(eval_mapped, 0.0d0, 1.0d0, epsabs, epsrel, res(1), &
                          res(2), ierr, key=15, limit=gk_limit)
        call record_error(ierr)

    contains

        recursive function eval_mapped(t) result(ft)
            real(dp), intent(in) :: t
            real(dp) :: ft
            real(dp) :: x_loc

            x_loc = x_low + (1.0d0 - t)/t
            ft = func1d_param0_user(x_loc)/(t*t)
        end function eval_mapped

    end function fint1d_param0_qagiu


    recursive function fint1d_param1_qagiu(func1d_param1_user, param1, x_low, &
            epsabs, epsrel) result(res)
        interface
            function func1d_param1_user(x, param1)
                use nrtype, only : dp
                real(dp) :: func1d_param1_user
                real(dp) :: x, param1
            end function func1d_param1_user
        end interface
        real(dp) :: param1
        real(dp) :: x_low
        real(dp) :: epsabs, epsrel
        real(dp), dimension(2) :: res

        integer :: ierr

        call integrate_gk(eval_mapped, 0.0d0, 1.0d0, epsabs, epsrel, res(1), &
                          res(2), ierr, key=15, limit=gk_limit)
        call record_error(ierr)

    contains

        recursive function eval_mapped(t) result(ft)
            real(dp), intent(in) :: t
            real(dp) :: ft
            real(dp) :: x_loc

            x_loc = x_low + (1.0d0 - t)/t
            ft = func1d_param1_user(x_loc, param1)/(t*t)
        end function eval_mapped

    end function fint1d_param1_qagiu


    subroutine record_error(ierr)
        integer, intent(in) :: ierr

        if (ierr .eq. 0) return

        if (ierr .lt. 0 .or. ierr .gt. max_err_code) then
            print *, 'integration_routines_mod.f90: unknown quadrature error code!'
            stop
        end if

        err_detected(ierr) = .true.
    end subroutine record_error


    subroutine disp_integration_error()
        integer :: k

        if (any(err_detected)) then
            print *, '-------------------------------------------------'
            do k = lbound(err_detected, 1), ubound(err_detected, 1)
                if (err_detected(k)) then
                    print *, "integration_routines_mod.f90: &
                        &Possible Warning from Integration Routines - Code = ", k
                end if
            end do
            print *, '-------------------------------------------------'
        end if

        err_detected = .false.
    end subroutine disp_integration_error

end module integration_routines_mod
