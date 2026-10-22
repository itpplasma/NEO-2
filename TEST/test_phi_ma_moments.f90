program test_phi_ma_moments
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use rkstep_mod, only: lag, leg
    use collop_compute, only: init_collop, compute_Minv, compute_exb_source, &
        M_transform, M_transform_inv, make_ortho, isw_relativistic
    use phi_ma_response_mod, only: normalize_phi_ma_j_hat
    implicit none
    real(dp) :: raw(0:2), projected(0:2), gram(0:2, 0:2), pi
    complex(dp) :: current(1), jhat(1)
    logical :: valid
    integer :: m, k

    lag = 2
    leg = 2
    isw_relativistic = 0
    call init_collop(1, 1, 0.0_dp, 0.0_dp)
    allocate(M_transform(0:2, 0:2), M_transform_inv(0:2, 0:2))
    make_ortho = .false.
    call compute_Minv(M_transform_inv)
    call compute_exb_source(raw)
    pi = acos(-1.0_dp)
    do m = 0, 2
        ! Independent Gaussian integrals for phi_m=x^m, not an endpoint probe.
        if (abs(raw(m) - gamma((m+4.0_dp)/2)/ &
            (2*pi**1.5_dp)) > 1.0e-11_dp) error stop 'x^-1 raw moment'
        do k = 0, 2
            gram(m,k) = gamma((m+k+5.0_dp)/2)/(2*pi**1.5_dp)
        end do
    end do
    make_ortho = .true.
    call compute_exb_source(projected)
    if (maxval(abs(matmul(gram, projected) - raw)) > 1.0e-10_dp) &
        error stop 'x^-1 orthogonalisation'

    ! Independently probed conversion from reference/L1/jhat_conversion.json.
    current = cmplx(1.0e7_dp, -2.0e6_dp, dp)
    call normalize_phi_ma_j_hat(current, -1.0_dp, 1.0e13_dp, &
        1.602176634e-9_dp, 9.10938356e-28_dp, 0.01_dp, &
        cmplx(1.0e-6_dp, 0.0_dp, dp), 0.22039696849234186_dp, jhat, valid)
    if (.not. valid) error stop 'J_hat unexpectedly undefined'
    if (abs(jhat(1) - cmplx(-34.59905551723764_dp, 6.919811103447528_dp, dp)) &
        > 1.0e-10_dp) error stop 'J_hat signed charge/thermal normalization'
    print *, 'All tests passed!'
end program test_phi_ma_moments
