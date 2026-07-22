program neo2_coordinate_probe
    use nrtype, only : dp
    use neo_input, only : psi_pr
    use neo_magfie, only : neo_magfie_calc, magfie_sarray, magfie_spline, &
        & magfie_result, iota_array, curr_pol_array, curr_tor_array
    implicit none

    real(dp) :: x(3), bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(dp) :: bcovar_s_hat_der(3), r_cgs, z_cgs, phi_cyl, b_cart(3)
    integer :: i

    x = [0.5_dp, 0.0_dp, 0.0_dp]
    do i = 1, 3
        call read_real_argument(i, x(i))
    end do

    allocate(magfie_sarray(1))
    magfie_sarray(1) = x(1)
    magfie_spline = 1
    magfie_result = 0

    call neo_magfie_calc(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl, &
        & bcovar_s_hat_der, r_cgs, z_cgs, phi_cyl, b_cart)

    write(*, '(a,es24.16)') 'NEO2_PROBE s=', x(1)
    write(*, '(a,es24.16)') 'NEO2_PROBE phi_b=', x(2)
    write(*, '(a,es24.16)') 'NEO2_PROBE theta_b=', x(3)
    write(*, '(a,es24.16)') 'NEO2_PROBE R_m=', r_cgs * 1.0e-2_dp
    write(*, '(a,es24.16)') 'NEO2_PROBE Z_m=', z_cgs * 1.0e-2_dp
    write(*, '(a,es24.16)') 'NEO2_PROBE phi_cyl_rad=', phi_cyl
    write(*, '(a,es24.16)') 'NEO2_PROBE x_m=', r_cgs * 1.0e-2_dp * cos(phi_cyl)
    write(*, '(a,es24.16)') 'NEO2_PROBE y_m=', r_cgs * 1.0e-2_dp * sin(phi_cyl)
    write(*, '(a,es24.16)') 'NEO2_PROBE z_m=', z_cgs * 1.0e-2_dp
    write(*, '(a,es24.16)') 'NEO2_PROBE B_T=', bmod
    write(*, '(a,es24.16)') 'NEO2_PROBE Bx_T=', b_cart(1)
    write(*, '(a,es24.16)') 'NEO2_PROBE By_T=', b_cart(2)
    write(*, '(a,es24.16)') 'NEO2_PROBE Bz_T=', b_cart(3)
    write(*, '(a,es24.16)') 'NEO2_PROBE sqrtg_native=', sqrtg
    write(*, '(a,es24.16)') 'NEO2_PROBE iota=', iota_array(1)
    write(*, '(a,es24.16)') 'NEO2_PROBE Bcov_phi_pre_handedness_Tm=', curr_pol_array(1)
    write(*, '(a,es24.16)') 'NEO2_PROBE Bcov_theta_pre_handedness_Tm=', curr_tor_array(1)
    write(*, '(a,es24.16)') 'NEO2_PROBE psi_pr_native=', psi_pr
    do i = 1, 3
        write(*, '(a,i0,a,es24.16)') 'NEO2_PROBE hctrvr(', i, ')=', hctrvr(i)
        write(*, '(a,i0,a,es24.16)') 'NEO2_PROBE hcovar(', i, ')=', hcovar(i)
    end do

contains

    subroutine read_real_argument(index, value)
        integer, intent(in) :: index
        real(dp), intent(inout) :: value
        character(len=128) :: text
        integer :: status

        call get_command_argument(index, text, status=status)
        if (status == 0 .and. len_trim(text) > 0) then
            read(text, *, iostat=status) value
            if (status /= 0) error stop 'invalid real-valued probe argument'
        end if
    end subroutine read_real_argument

end program neo2_coordinate_probe
