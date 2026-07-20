program test_helical_response
    use helical_response_mod, only: D_CURRENT_UNIT, &
        dimensional_current_coefficients, reconstruct_complex_response
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    real(dp) :: gamma(3, 3), gamma_current(2), d_current(2)
    real(dp) :: real_run(2), imaginary_run(2), baseline(2)
    complex(dp) :: response(2)

    gamma = 0.0_dp
    gamma(3, 1:2) = [2.0_dp, -3.0_dp]
    call dimensional_current_coefficients(gamma, 4.0_dp, -2.0_dp, &
        gamma_current, d_current)

    if (any(gamma_current /= [2.0_dp, -3.0_dp])) &
        error stop 'FAIL: current channels do not map to gamma rows 3 and columns 1:2'
    if (maxval(abs(d_current - [-4.9931712191872093e20_dp, &
        7.489756828780813e20_dp])) > 1.0e6_dp) &
        error stop 'FAIL: dimensional D31/D32 normalization changed'
    if (D_CURRENT_UNIT /= 'cm^2 G/s') &
        error stop 'FAIL: D31/D32 unit metadata is wrong'
    real_run = [5.0_dp, 8.0_dp]
    imaginary_run = [1.0_dp, 7.0_dp]
    baseline = [2.0_dp, 3.0_dp]
    call reconstruct_complex_response(real_run, imaginary_run, baseline, response)
    if (maxval(abs(response - [cmplx(3.0_dp, 1.0_dp, dp), &
        cmplx(5.0_dp, -4.0_dp, dp)])) > 1.0e-14_dp) &
        error stop 'FAIL: real-imaginary quadrature reconstruction has the wrong sign'

    print *, 'All tests passed!'
end program test_helical_response
