program test_periodic_spline_underflow
    use nrtype, only: dp
    use spline_mod, only: splper
    use, intrinsic :: ieee_exceptions, only: ieee_set_halting_mode, ieee_underflow
    implicit none

    call ieee_set_halting_mode(ieee_underflow, .true.)
    call check_grid(17)
    call check_grid(1025)
contains
    subroutine check_grid(n)
        integer, intent(in) :: n
        real(dp) :: h, theta, eigenvalue1, eigenvalue2, expected, tolerance
        real(dp) :: values(n), bi(n), ci(n), di(n)
        integer :: j

        h = 8.0_dp*atan(1.0_dp)/real(n-1, dp)
        do j = 1, n
            theta = real(j-1, dp)*h
            values(j) = cos(theta) + 0.25_dp*sin(2.0_dp*theta)
        end do
        values(n) = values(1)
        call splper(n, h, values, bi, ci, di)
        ! Independent Fourier diagonalization of the periodic cubic system:
        ! (c[j-1] + 4*c[j] + c[j+1])*h**2 = 3*(y[j+1]-2*y[j]+y[j-1]).
        eigenvalue1 = -6.0_dp*sin(h/2.0_dp)**2/(h**2*(2.0_dp+cos(h)))
        eigenvalue2 = -6.0_dp*sin(h)**2/(h**2*(2.0_dp+cos(2.0_dp*h)))
        tolerance = 64.0_dp*epsilon(1.0_dp)/h**2
        do j = 1, n
            theta = real(j-1, dp)*h
            expected = eigenvalue1*cos(theta) + 0.25_dp*eigenvalue2*sin(2.0_dp*theta)
            if (abs(ci(j)-expected) > tolerance) error stop "periodic cubic Fourier oracle"
        end do
        if (maxval(abs([bi(n)-bi(1), ci(n)-ci(1), di(n)-di(1)])) > epsilon(1.0_dp)) &
            error stop "periodic spline seam"
    end subroutine check_grid
end program test_periodic_spline_underflow
