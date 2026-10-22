program test_phi1_qn
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use phi1_mod, only: phi1_assemble_qn, phi1_solve_qn
    implicit none
    real(dp) :: charge(2), density(2), temperature(2), chi_expected, chi
    real(dp) :: residual, rcond
    complex(dp) :: gram(2,2), response(2,2,2), dB(2,2), matrix_expected(2,2)
    complex(dp) :: expected(2), matrix(2,2), rhs(2), phi(2)
    integer :: unit, info
    character(1024) :: fixture
    namelist /qn_fixture/ charge, density, temperature, gram, response, dB, &
        chi_expected, matrix_expected, expected
    call get_command_argument(1,fixture)
    open(newunit=unit,file=trim(fixture),status='old')
    read(unit,nml=qn_fixture)
    close(unit)
    call phi1_assemble_qn(charge,density,temperature,gram,response,dB, &
        matrix,rhs,chi,info)
    if (info /= 0 .or. abs(chi-chi_expected) > 2e-14_dp) error stop 'assembly'
    if (maxval(abs(matrix-matrix_expected)) > 2e-14_dp) error stop 'matrix'
    call solve()
    if (maxval(abs(phi-expected)) > 2e-14_dp) error stop 'complex solution'
    if (maxval(abs(matmul(matrix,phi)-rhs)) > 2e-14_dp) error stop 'residual'
    if (residual > 2e-14_dp) error stop 'reported residual'
    rhs = 0.0_dp
    call solve()
    if (any(phi /= (0.0_dp,0.0_dp))) error stop 'zero forcing'
    matrix = 0.0_dp
    call phi1_solve_qn(matrix,rhs,1e-12_dp,1e-12_dp,phi,residual,rcond,info)
    if (info <= 0) error stop 'singular rejection'
    matrix(1,1) = 1.0_dp
    matrix(2,2) = 1e-14_dp
    call phi1_solve_qn(matrix,rhs,1e-12_dp,1e-12_dp,phi,residual,rcond,info)
    if (info /= -2) error stop 'conditioning rejection'
contains
    subroutine solve()
        call phi1_solve_qn(matrix,rhs,1e-12_dp,1e-12_dp,phi,residual,rcond,info)
        if (info /= 0) error stop 'solve'
    end subroutine
end program
