program test_spline_coverage
    !> Comprehensive coverage test for spline implementation
    !> This test exercises all code paths in the sparse spline implementation
    !> to ensure adequate test coverage for codecov
    use nrtype, only: I4B, DP
    use splinecof3_direct_sparse_mod, only: splinecof3_direct_sparse, splinecof3_assemble_matrix
    use ieee_arithmetic, only: ieee_is_finite
    implicit none
    
    logical :: all_tests_passed = .true.
    
    write(*,'(A)') '=== Comprehensive Spline Coverage Tests ==='
    write(*,'(A)') 'Testing all code paths in sparse spline implementation'
    write(*,'(A)') ''
    
    ! Test 1: Matrix assembly function
    call test_matrix_assembly()
    
    ! Test 2: All boundary condition combinations
    call test_all_boundary_conditions()
    
    ! Test 3: Edge cases and error conditions
    call test_edge_cases()
    
    ! Test 4: Different lambda weight scenarios
    call test_lambda_scenarios()
    
    ! Test 5: Non-zero m parameter scenarios
    call test_m_parameter_scenarios()
    
    ! Test 6: Various data patterns
    call test_data_patterns()
    
    if (all_tests_passed) then
        write(*,'(A)') ''
        write(*,'(A)') 'All coverage tests PASSED!'
        write(*,'(A)') 'Comprehensive test coverage achieved.'
        stop 0
    else
        write(*,'(A)') ''
        write(*,'(A)') 'Some coverage tests FAILED!'
        stop 1
    end if
    
contains

    function test_function(x, m) result(f_val)
        real(DP), intent(in) :: x, m
        real(DP) :: f_val
        f_val = exp(-x*x) * (1.0_DP + m*x)
    end function test_function

    subroutine test_matrix_assembly()
        integer(I4B), parameter :: n = 5
        real(DP) :: x(n), y(n)
        real(DP) :: c1, cn, m
        real(DP) :: lambda1(n)
        integer(I4B) :: indx(n), sw1, sw2
        integer(I4B) :: nrow, ncol, nnz
        integer(I4B), allocatable :: irow_coo(:), icol_coo(:)
        real(DP), allocatable :: val_coo(:), rhs(:)
        integer(I4B) :: i
        
        write(*,'(A)') 'Test 1: Matrix assembly function'
        
        ! Setup test data
        x = [0.0_DP, 0.5_DP, 1.0_DP, 1.5_DP, 2.0_DP]
        y = x*x + 0.1_DP*sin(x)
        indx = [(i, i=1,n)]
        lambda1 = 1.0_DP
        c1 = 0.0_DP
        cn = 0.0_DP
        sw1 = 2  ! Natural start
        sw2 = 4  ! Natural end
        m = 0.5_DP
        
        ! Test matrix assembly
        call splinecof3_assemble_matrix(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                                       m, test_function, nrow, ncol, nnz, &
                                       irow_coo, icol_coo, val_coo, rhs)
        
        ! Verify matrix dimensions
        if (nrow /= ncol .or. nrow /= 33) then
            write(*,'(A,3I0)') '  FAILED: Unexpected matrix dimensions: ', nrow, ncol, nnz
            all_tests_passed = .false.
        else
            write(*,'(A,I0,A,I0)') '  PASSED: Matrix assembly (', nrow, 'x', ncol, ')'
        end if
        
        ! Verify non-zero count is reasonable
        if (nnz < 50 .or. nnz > 500) then
            write(*,'(A,I0)') '  FAILED: Unexpected number of non-zeros: ', nnz
            all_tests_passed = .false.
        else
            write(*,'(A,I0)') '  PASSED: Non-zero count reasonable: ', nnz
        end if
        
        deallocate(irow_coo, icol_coo, val_coo, rhs)
    end subroutine test_matrix_assembly

    subroutine test_all_boundary_conditions()
        integer(I4B), parameter :: n = 4
        real(DP) :: x(n), y(n)
        real(DP) :: c1, cn, m
        real(DP) :: lambda1(n)
        integer(I4B) :: indx(n)
        real(DP) :: a(n), b(n), c(n), d(n)
        integer(I4B) :: sw1, sw2, i, test_num
        
        write(*,'(A)') ''
        write(*,'(A)') 'Test 2: All boundary condition combinations'
        
        ! Setup test data
        x = [0.0_DP, 1.0_DP, 2.0_DP, 3.0_DP]
        y = [1.0_DP, 2.0_DP, 1.5_DP, 3.0_DP]
        indx = [(i, i=1,n)]
        lambda1 = 1.0_DP
        m = 0.0_DP
        
        test_num = 0
        
        ! Test all valid boundary condition combinations
        do sw1 = 1, 4
            do sw2 = 1, 4
                if (sw1 == sw2) cycle  ! Skip invalid combinations
                
                test_num = test_num + 1
                
                ! Set appropriate boundary values
                select case(sw1)
                case(1); c1 = 0.5_DP   ! First derivative at start
                case(2); c1 = 0.0_DP   ! Second derivative at start  
                case(3); c1 = -0.3_DP  ! First derivative at end (swap)
                case(4); c1 = 0.2_DP   ! Second derivative at end (swap)
                end select
                
                select case(sw2)
                case(1); cn = 0.7_DP   ! First derivative at start (swap)
                case(2); cn = 0.0_DP   ! Second derivative at start (swap)
                case(3); cn = -0.5_DP  ! First derivative at end
                case(4); cn = 0.0_DP   ! Second derivative at end
                end select
                
                call splinecof3_direct_sparse(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                                             a, b, c, d, m, test_function)
                
                ! Basic sanity checks
                if (any(.not. ieee_is_finite(a(1:n-1))) .or. &
                    any(.not. ieee_is_finite(b(1:n-1))) .or. &
                    any(.not. ieee_is_finite(c(1:n-1))) .or. &
                    any(.not. ieee_is_finite(d(1:n-1)))) then
                    write(*,'(A,I0,A,2I0)') '  FAILED: Test ', test_num, ' (sw1,sw2)=(', sw1, sw2, ') - non-finite coefficients'
                    all_tests_passed = .false.
                else
                    write(*,'(A,I0,A,2I0)') '  PASSED: Test ', test_num, ' (sw1,sw2)=(', sw1, sw2, ')'
                end if
            end do
        end do
    end subroutine test_all_boundary_conditions

    subroutine test_edge_cases()
        integer(I4B), parameter :: n = 3  ! Minimum size
        real(DP) :: x(n), y(n)
        real(DP) :: c1, cn, m
        real(DP) :: lambda1(n)
        integer(I4B) :: indx(n)
        real(DP) :: a(n), b(n), c(n), d(n)
        integer(I4B) :: i
        
        write(*,'(A)') ''
        write(*,'(A)') 'Test 3: Edge cases and boundary scenarios'
        
        ! Test minimum size problem
        x = [0.0_DP, 0.5_DP, 1.0_DP]
        y = [1.0_DP, 1.5_DP, 2.0_DP]
        indx = [(i, i=1,n)]
        lambda1 = 1.0_DP
        c1 = 0.0_DP
        cn = 0.0_DP
        m = 0.0_DP
        
        call splinecof3_direct_sparse(x, y, c1, cn, lambda1, indx, 2, 4, &
                                     a, b, c, d, m, test_function)
        
        if (any(.not. ieee_is_finite(a(1:n-1)))) then
            write(*,'(A)') '  FAILED: Minimum size problem'
            all_tests_passed = .false.
        else
            write(*,'(A)') '  PASSED: Minimum size problem'
        end if
        
        ! Test with very large boundary values that should be reset
        c1 = 1.0e35_DP
        cn = 1.0e35_DP
        call splinecof3_direct_sparse(x, y, c1, cn, lambda1, indx, 1, 3, &
                                     a, b, c, d, m, test_function)
        
        ! Should be reset to 0 due to the >1e30 check in the sparse implementation
        if (abs(c1) < 1.0e30_DP .and. abs(cn) < 1.0e30_DP) then
            write(*,'(A)') '  PASSED: Large boundary value reset'
        else
            write(*,'(A)') '  FAILED: Large boundary value reset'
            all_tests_passed = .false.
        end if
    end subroutine test_edge_cases

    subroutine test_lambda_scenarios()
        integer(I4B), parameter :: n = 5
        real(DP) :: x(n), y(n)
        real(DP) :: c1, cn, m
        real(DP) :: lambda1(n)
        integer(I4B) :: indx(n)
        real(DP) :: a(n), b(n), c(n), d(n)
        integer(I4B) :: i
        
        write(*,'(A)') ''
        write(*,'(A)') 'Test 4: Different lambda weight scenarios'
        
        x = [0.0_DP, 0.25_DP, 0.75_DP, 1.25_DP, 2.0_DP]
        y = sin(x) + 0.1_DP*x*x
        indx = [(i, i=1,n)]
        c1 = 0.0_DP
        cn = 0.0_DP
        m = 0.0_DP
        
        ! Test 1: Negative lambda (triggers optimal lambda calculation)
        lambda1 = -1.0_DP
        call splinecof3_direct_sparse(x, y, c1, cn, lambda1, indx, 2, 4, &
                                     a, b, c, d, m, test_function)
        
        if (any(.not. ieee_is_finite(a(1:n-1)))) then
            write(*,'(A)') '  FAILED: Negative lambda (optimal weights)'
            all_tests_passed = .false.
        else
            write(*,'(A)') '  PASSED: Negative lambda (optimal weights)'
        end if
        
        ! Test 2: Non-uniform lambda weights
        lambda1 = [0.1_DP, 0.3_DP, 0.7_DP, 0.9_DP, 0.5_DP]
        call splinecof3_direct_sparse(x, y, c1, cn, lambda1, indx, 2, 4, &
                                     a, b, c, d, m, test_function)
        
        if (any(.not. ieee_is_finite(a(1:n-1)))) then
            write(*,'(A)') '  FAILED: Non-uniform lambda weights'
            all_tests_passed = .false.
        else
            write(*,'(A)') '  PASSED: Non-uniform lambda weights'
        end if
        
        ! Test 3: Very small lambda weights (near pure interpolation)
        lambda1 = 1.0e-6_DP
        call splinecof3_direct_sparse(x, y, c1, cn, lambda1, indx, 2, 4, &
                                     a, b, c, d, m, test_function)
        
        if (any(.not. ieee_is_finite(a(1:n-1)))) then
            write(*,'(A)') '  FAILED: Very small lambda weights'
            all_tests_passed = .false.
        else
            write(*,'(A)') '  PASSED: Very small lambda weights'
        end if
    end subroutine test_lambda_scenarios

    subroutine test_m_parameter_scenarios()
        integer(I4B), parameter :: n = 6
        real(DP) :: x(n), y(n)
        real(DP) :: c1, cn, m
        real(DP) :: lambda1(n)
        integer(I4B) :: indx(n)
        real(DP) :: a(n), b(n), c(n), d(n)
        integer(I4B) :: i
        
        write(*,'(A)') ''
        write(*,'(A)') 'Test 5: Non-zero m parameter scenarios'
        
        x = [0.0_DP, 0.2_DP, 0.6_DP, 1.0_DP, 1.4_DP, 2.0_DP]
        y = x*x + 0.5_DP*cos(x)
        indx = [(i, i=1,n)]
        lambda1 = 0.8_DP
        c1 = 0.5_DP
        cn = -0.3_DP
        
        ! Test different m values
        do i = 1, 5
            m = real(i-3, DP) * 0.5_DP  ! -1.0, -0.5, 0.0, 0.5, 1.0
            
            call splinecof3_direct_sparse(x, y, c1, cn, lambda1, indx, 1, 3, &
                                         a, b, c, d, m, test_function)
            
            if (any(.not. ieee_is_finite(a(1:n-1)))) then
                write(*,'(A,F4.1)') '  FAILED: m parameter = ', m
                all_tests_passed = .false.
            else
                write(*,'(A,F4.1)') '  PASSED: m parameter = ', m
            end if
        end do
    end subroutine test_m_parameter_scenarios

    subroutine test_data_patterns()
        integer(I4B), parameter :: n = 7
        real(DP) :: x(n), y(n)
        real(DP) :: c1, cn, m
        real(DP) :: lambda1(n)
        integer(I4B) :: indx(n)
        real(DP) :: a(n), b(n), c(n), d(n)
        integer(I4B) :: i
        ! Variables for smaller subset test
        integer(I4B), parameter :: n_small = 5
        real(DP) :: x_small(n_small), y_small(n_small), lambda1_small(n_small)
        integer(I4B) :: indx_small(n_small)
        real(DP) :: a_small(n_small), b_small(n_small), c_small(n_small), d_small(n_small)
        
        write(*,'(A)') ''
        write(*,'(A)') 'Test 6: Various data patterns'
        
        x = [0.0_DP, 0.1_DP, 0.3_DP, 0.7_DP, 1.2_DP, 1.8_DP, 2.5_DP]
        lambda1 = 1.0_DP
        c1 = 0.0_DP
        cn = 0.0_DP
        m = 0.0_DP
        
        ! Test 1: Linear data
        y = 2.0_DP*x + 1.0_DP
        indx = [(i, i=1,n)]
        call splinecof3_direct_sparse(x, y, c1, cn, lambda1, indx, 2, 4, &
                                     a, b, c, d, m, test_function)
        
        if (any(.not. ieee_is_finite(a(1:n-1)))) then
            write(*,'(A)') '  FAILED: Linear data pattern'
            all_tests_passed = .false.
        else
            write(*,'(A)') '  PASSED: Linear data pattern'
        end if
        
        ! Test 2: Oscillatory data
        y = sin(2.0_DP*x) * exp(-0.5_DP*x)
        call splinecof3_direct_sparse(x, y, c1, cn, lambda1, indx, 2, 4, &
                                     a, b, c, d, m, test_function)
        
        if (any(.not. ieee_is_finite(a(1:n-1)))) then
            write(*,'(A)') '  FAILED: Oscillatory data pattern'
            all_tests_passed = .false.
        else
            write(*,'(A)') '  PASSED: Oscillatory data pattern'
        end if
        
        ! Test 3: Subset of data points
        x_small = x(1:n_small)
        y_small = x_small*x_small + 0.1_DP*x_small
        lambda1_small = 1.0_DP
        indx_small = [(i, i=1,n_small)]
        call splinecof3_direct_sparse(x_small, y_small, c1, cn, lambda1_small, indx_small, 2, 4, &
                                     a_small, b_small, c_small, d_small, m, test_function)
        
        if (any(.not. ieee_is_finite(a_small(1:n_small-1)))) then
            write(*,'(A)') '  FAILED: Subset of data points'
            all_tests_passed = .false.
        else
            write(*,'(A)') '  PASSED: Subset of data points'
        end if
        
        ! Test 4: Constant data
        y = 3.14_DP
        indx = [(i, i=1,n)]
        call splinecof3_direct_sparse(x, y, c1, cn, lambda1, indx, 2, 4, &
                                     a, b, c, d, m, test_function)
        
        if (any(.not. ieee_is_finite(a(1:n-1)))) then
            write(*,'(A)') '  FAILED: Constant data pattern'
            all_tests_passed = .false.
        else
            write(*,'(A)') '  PASSED: Constant data pattern'
        end if
    end subroutine test_data_patterns

end program test_spline_coverage