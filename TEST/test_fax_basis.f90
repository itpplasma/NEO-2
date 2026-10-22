program test_fax_basis
    use nrtype, only: dp
    use fax_export_mod, only: fax_basis
    use collop, only: collop_base_exp, phi_exp, d_phi_exp, dd_phi_exp
    use collop_laguerre, only: init_phi_laguerre, phi_laguerre
    use collop_bspline, only: init_phi_bspline, phi_bspline, &
        d_phi_bspline, dd_phi_bspline
    use collisionality_mod, only: isw_lorentz, isw_relativistic, &
        collop_bspline_order, collop_bspline_dist, collop_bspline_taylor, phi_x_max
    use rkstep_mod, only: lag
    implicit none
    real(dp), allocatable :: coeff(:, :, :), edges(:), centers(:)
    real(dp) :: x, value, ref, maxerr
    integer :: basis, i, m, seg, power
    logical :: extend

    lag = 6
    isw_lorentz = 0
    isw_relativistic = 0
    collop_bspline_order = 3
    collop_bspline_dist = 1.0_dp
    collop_bspline_taylor = .true.
    phi_x_max = 4.0_dp
    maxerr = 0.0_dp
    do basis = 0, 11, 11
        collop_base_exp = basis
        if (basis == 0) then
            call init_phi_laguerre(lag, 1)
            phi_exp => phi_laguerre
        else
            call init_phi_bspline(lag, 1)
            phi_exp => phi_bspline
            d_phi_exp => d_phi_bspline
            dd_phi_exp => dd_phi_bspline
        end if
        call fax_basis(coeff, edges, centers, extend)
        do i = 0, 100
            x = 0.05_dp*real(i, dp)
            seg = 1
            do while (seg < size(centers))
                if (x < edges(seg+1)) exit
                seg = seg+1
            end do
            do m = 0, lag
                value = coeff(ubound(coeff, 1), m, seg)
                do power = ubound(coeff, 1)-1, 0, -1
                    value = value*(x-centers(seg))+coeff(power, m, seg)
                end do
                ref = phi_exp(m, x)
                maxerr = max(maxerr, abs(value-ref)/max(1.0_dp, abs(ref)))
            end do
        end do
        deallocate(coeff, edges, centers)
    end do
    print '(a, es24.16)', 'f_AX basis max normalized error: ', maxerr
    if (maxerr > 2.0e-11_dp) error stop 'f_AX basis reconstruction'
    print '(a)', 'All tests passed!'
end program test_fax_basis
