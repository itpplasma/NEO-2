program test_sparse_method_controls
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use sparse_mod, only: sparse_solve, sparse_solve_method, &
        sparse_method_from_environment
    implicit none
    real(dp) :: matrix(3,3), rhs(3), exact(3), product(3), backward_error
    complex(dp) :: complex_matrix(3,3), complex_rhs(3), complex_exact(3)
    complex(dp) :: complex_product(3)
    integer :: expected
    character(len=32) :: argument

    call get_command_argument(1, argument)
    sparse_solve_method = 3
    call sparse_method_from_environment()
    ! An accepted invalid selector must exit successfully so WILL_FAIL rejects it.
    if (trim(argument) == 'invalid') stop 0
    read(argument, *) expected
    if (sparse_solve_method /= expected) error stop 'wrong sparse method'

    ! Independent dense matvecs bind fidelity of the real and complex solve paths.
    matrix = reshape([4.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 5.0_dp, &
                      1.0_dp, 0.0_dp, 3.0_dp, 6.0_dp], [3,3])
    exact = [1.0_dp, -2.0_dp, 3.0_dp]
    rhs = matmul(matrix, exact)
    product = rhs
    call sparse_solve(matrix, rhs)
    backward_error = sum(abs(product-matmul(matrix,rhs))) / &
        (sum(abs(product))+sum(abs(matmul(matrix,rhs))))
    if (.not. ieee_is_finite(backward_error)) error stop 'nonfinite real residual'
    if (backward_error > 1.0e-12_dp) error stop 'real inverse scaled residual'
    if (maxval(abs(rhs-exact)) > 1.0e-12_dp) error stop 'real inverse solution'

    complex_matrix = cmplx(matrix, 0.1_dp*transpose(matrix), dp)
    complex_exact = cmplx(exact, [0.5_dp, 1.0_dp, -0.5_dp], dp)
    complex_rhs = matmul(complex_matrix, complex_exact)
    complex_product = complex_rhs
    call sparse_solve(complex_matrix, complex_rhs)
    backward_error = sum(abs(complex_product-matmul(complex_matrix,complex_rhs))) / &
        (sum(abs(complex_product))+sum(abs(matmul(complex_matrix,complex_rhs))))
    if (.not. ieee_is_finite(backward_error)) error stop 'nonfinite complex residual'
    if (backward_error > 1.0e-12_dp) error stop 'complex inverse scaled residual'
    if (maxval(abs(complex_rhs-complex_exact)) > 1.0e-12_dp) &
        error stop 'complex inverse solution'
    print *, 'All tests passed!'
end program test_sparse_method_controls
