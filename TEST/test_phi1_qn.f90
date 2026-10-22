program test_phi1_qn
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use phi1_mod, only: phi1_basis, phi1_fieldline_weights, phi1_gram, &
        phi1_assemble_qn, phi1_solve_qn
    implicit none
    real(dp), parameter :: tol = 2.0e-14_dp
    real(dp) :: charge(2), density(2), temperature(2), chi_expected, chi
    real(dp) :: residual, rcond, weights(3)
    complex(dp) :: gram(2,2), response(2,2,2), dB(2,2), matrix_expected(2,2)
    complex(dp) :: expected(2), matrix(2,2), rhs(2), phi(2), phase(2)
    complex(dp) :: psi(3,2), small_psi(2,2), measured_gram(2,2)
    integer :: unit, info, l, j
    character(len=1024) :: fixture
    namelist /qn_fixture/ charge, density, temperature, gram, response, dB, &
        chi_expected, matrix_expected, expected

    call get_command_argument(1, fixture)
    open(newunit=unit, file=trim(fixture), status='old', action='read')
    read(unit,nml=qn_fixture)
    close(unit)

    ! The reference case retains adiabatic electrons in chi, with zero response.
    call phi1_assemble_qn(charge, density, temperature, gram, response, dB, &
                          matrix, rhs, chi, info)
    call require(info == 0, 'reference assembly')
    call require(abs(chi - chi_expected) < tol, 'adiabatic susceptibility')
    call require(maxval(abs(matrix - matrix_expected)) < tol, 'reference matrix')
    call require(maxval(abs(rhs - [-1.0_dp,2.0_dp])) < tol, 'RHS sign')
    call check_solution('adiabatic reference')

    ! Both species now contribute, preserving exactly the same charge sums.
    response(:,:,2) = -0.5_dp*response(:,:,1)
    response(:,:,1) = 0.5_dp*response(:,:,1)
    dB(:,2) = -0.5_dp*dB(:,1)
    dB(:,1) = 0.5_dp*dB(:,1)
    call phi1_assemble_qn(charge, density, temperature, gram, response, dB, &
                          matrix, rhs, chi, info)
    call require(info == 0, 'two kinetic species assembly')
    call require(maxval(abs(matrix - matrix_expected)) < tol, 'charge sum')
    call check_solution('two kinetic species')

    ! Change basis by complex phases: a non-Hermitian complex closure with a
    ! known solution, checking that the solver keeps both real and imaginary parts.
    phase = [cmplx(0.6_dp,0.8_dp,dp), cmplx(0.0_dp,-1.0_dp,dp)]
    do j = 1, 2
        do l = 1, 2
            matrix(l,j) = conjg(phase(l))*matrix_expected(l,j)*phase(j)
        end do
    end do
    rhs = conjg(phase)*rhs
    expected = conjg(phase)*expected
    call check_solution('complex phase covariance')

    rhs = 0.0_dp
    expected = 0.0_dp
    call check_solution('zero forcing')

    matrix = 0.0_dp
    rhs = 1.0_dp
    call phi1_solve_qn(matrix, rhs, 1.0e-12_dp, 1.0e-12_dp, &
                       phi, residual, rcond, info)
    call require(info > 0, 'singular system rejection')
    call require(all(abs(phi) == 0.0_dp), 'no coefficients after rejection')

    matrix(1,1) = 1.0_dp
    matrix(2,2) = 1.0e-14_dp
    call phi1_solve_qn(matrix, rhs, 1.0e-12_dp, 1.0e-12_dp, &
                       phi, residual, rcond, info)
    call require(info == -2, 'ill-conditioned system rejection')
    call require(rcond < 1.0e-12_dp, 'small reciprocal condition estimate')

    call phi1_basis([0.3_dp,0.8_dp,1.3_dp], 0.2_dp, 0.3_dp, 0.7_dp, &
                     [-1,2], -2, psi)
    call require(maxval(abs(psi(:,1) - &
        exp(cmplx(0.0_dp,[-0.8_dp,-2.15_dp,-3.5_dp],dp)))) < tol, &
        'field-line basis first mode and starting angles')
    call require(maxval(abs(psi(:,2) - &
        exp(cmplx(0.0_dp,[-0.2_dp,-0.5_dp,-0.8_dp],dp)))) < tol, &
        'field-line basis second mode')
    call phi1_fieldline_weights([0.0_dp,1.0_dp,3.0_dp], &
        [1.0_dp,2.0_dp,1.0_dp], [-1.0_dp,-1.0_dp,-1.0_dp], weights, info)
    call require(info == 0, 'positive field-line quadrature')
    call require(maxval(abs(weights - [2.0_dp,3.0_dp,4.0_dp]/9.0_dp)) < tol, &
        'nonuniform geometry weights')
    small_psi(:,1) = 1.0_dp
    small_psi(:,2) = [cmplx(1.0_dp,0.0_dp,dp), cmplx(0.0_dp,1.0_dp,dp)]
    call phi1_gram(small_psi, [1.0_dp,3.0_dp], measured_gram, info)
    call require(info == 0, 'Gram assembly')
    call require(abs(measured_gram(1,2) - cmplx(0.25_dp,0.75_dp,dp)) < tol, &
        'nonorthogonal weighted Gram entry')
    call require(abs(measured_gram(2,1) - conjg(measured_gram(1,2))) < tol, &
        'Gram conjugation')
    call require(abs(measured_gram(1,1) - 1.0_dp) < tol, 'Gram normalization')
    call require(abs(measured_gram(2,2) - 1.0_dp) < tol, 'Gram normalization')

    print *, 'All tests passed!'

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (condition) return
        print *, 'FAIL: ', message
        error stop 1
    end subroutine require

    subroutine check_solution(message)
        character(len=*), intent(in) :: message

        call phi1_solve_qn(matrix, rhs, 1.0e-12_dp, 1.0e-12_dp, &
                           phi, residual, rcond, info)
        call require(info == 0, message)
        call require(maxval(abs(phi - expected)) < tol, 'reference coefficients')
        call require(maxval(abs(matmul(matrix,phi) - rhs)) < tol, &
            'independent unscaled residual')
        call require(residual < tol, 'reported residual')
        call require(rcond > 1.0e-2_dp, 'well-conditioned reference system')
    end subroutine check_solution

end program test_phi1_qn
