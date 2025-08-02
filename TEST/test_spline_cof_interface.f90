program test_spline_cof_interface
    !> Test spline_cof.f90 interface functions to improve patch coverage
    !> Targets the modified spline_cof.f90 (15.95% coverage) with real usage scenarios
    use nrtype, only: I4B, DP
    use spline_mod, only: splinecof3_a
    implicit none
    
    logical :: all_tests_passed = .true.
    
    write(*,'(A)') '=== Spline Interface Tests ==='
    write(*,'(A)') 'Testing spline_cof.f90 interface for patch coverage improvement'
    write(*,'(A)') ''
    
    ! Test 1: Natural spline boundary conditions
    call test_natural_splines()
    
    ! Test 2: Clamped spline boundary conditions
    call test_clamped_splines()
    
    ! Test 3: Mixed boundary conditions
    call test_mixed_boundaries()
    
    ! Test 4: Edge case with minimum data points
    call test_minimum_points()
    
    if (all_tests_passed) then
        write(*,'(A)') ''
        write(*,'(A)') 'All spline interface tests PASSED!'
        stop 0
    else
        write(*,'(A)') ''
        write(*,'(A)') 'Some spline interface tests FAILED!'
        stop 1
    end if
    
contains

    function weight_function(x, m) result(w)
        real(DP), intent(in) :: x, m
        real(DP) :: w
        w = 1.0_DP  ! Uniform weights
    end function weight_function

    subroutine test_natural_splines()
        integer(I4B), parameter :: n = 5
        real(DP) :: x(n), y(n), lambda(n)
        integer(I4B) :: indx(n)
        real(DP) :: a(n), b(n), c(n), d(n)
        real(DP) :: c1, cn, m
        integer(I4B) :: i
        
        write(*,'(A)') 'Test 1: Natural spline boundary conditions (sw1=2, sw2=4)'
        
        ! Test data: smooth curve
        x = [0.0_DP, 1.0_DP, 2.0_DP, 3.0_DP, 4.0_DP]
        y = [1.0_DP, 1.5_DP, 2.2_DP, 2.8_DP, 3.1_DP]
        indx = [(i, i=1,n)]
        lambda = 0.7_DP
        c1 = 0.0_DP  ! Natural boundary (second derivative = 0)
        cn = 0.0_DP  ! Natural boundary (second derivative = 0)
        m = 0.0_DP
        
        call splinecof3_a(x, y, c1, cn, lambda, indx, 2, 4, a, b, c, d, m, weight_function)
        
        ! Verify spline coefficients are finite
        if (all_finite(a) .and. all_finite(b) .and. &
            all_finite(c) .and. all_finite(d)) then
            write(*,'(A)') '  PASSED: Natural spline computed successfully'
        else
            write(*,'(A)') '  FAILED: Invalid spline coefficients'
            all_tests_passed = .false.
        end if
    end subroutine test_natural_splines

    subroutine test_clamped_splines()
        integer(I4B), parameter :: n = 4
        real(DP) :: x(n), y(n), lambda(n)
        integer(I4B) :: indx(n)
        real(DP) :: a(n), b(n), c(n), d(n)
        real(DP) :: c1, cn, m
        integer(I4B) :: i
        
        write(*,'(A)') 'Test 2: Clamped spline boundary conditions (sw1=1, sw2=3)'
        
        ! Test data: quadratic-like curve
        x = [0.0_DP, 1.0_DP, 2.0_DP, 3.0_DP]
        y = [0.0_DP, 1.0_DP, 4.0_DP, 9.0_DP]
        indx = [(i, i=1,n)]
        lambda = 0.9_DP
        c1 = 1.0_DP   ! First derivative at start
        cn = 5.0_DP   ! First derivative at end
        m = 0.5_DP
        
        call splinecof3_a(x, y, c1, cn, lambda, indx, 1, 3, a, b, c, d, m, weight_function)
        
        if (all_finite(a) .and. all_finite(b) .and. &
            all_finite(c) .and. all_finite(d)) then
            write(*,'(A)') '  PASSED: Clamped spline computed successfully'
        else
            write(*,'(A)') '  FAILED: Invalid clamped spline coefficients'
            all_tests_passed = .false.
        end if
    end subroutine test_clamped_splines

    subroutine test_mixed_boundaries()
        integer(I4B), parameter :: n = 6
        real(DP) :: x(n), y(n), lambda(n)
        integer(I4B) :: indx(n)
        real(DP) :: a(n), b(n), c(n), d(n)
        real(DP) :: c1, cn, m
        integer(I4B) :: i
        
        write(*,'(A)') 'Test 3: Mixed boundary conditions (sw1=1, sw2=4)'
        
        ! Test data: exponential-like decay
        x = [0.0_DP, 0.5_DP, 1.0_DP, 1.5_DP, 2.0_DP, 2.5_DP]
        y = [3.0_DP, 2.2_DP, 1.6_DP, 1.2_DP, 0.9_DP, 0.7_DP]
        indx = [(i, i=1,n)]
        lambda = 0.6_DP
        c1 = -1.2_DP  ! First derivative at start
        cn = 0.0_DP   ! Second derivative at end
        m = 0.2_DP
        
        call splinecof3_a(x, y, c1, cn, lambda, indx, 1, 4, a, b, c, d, m, weight_function)
        
        if (all_finite(a) .and. all_finite(b) .and. &
            all_finite(c) .and. all_finite(d)) then
            write(*,'(A)') '  PASSED: Mixed boundary spline computed successfully'
        else
            write(*,'(A)') '  FAILED: Invalid mixed boundary spline coefficients'
            all_tests_passed = .false.
        end if
    end subroutine test_mixed_boundaries

    subroutine test_minimum_points()
        integer(I4B), parameter :: n = 3
        real(DP) :: x(n), y(n), lambda(n)
        integer(I4B) :: indx(n)
        real(DP) :: a(n), b(n), c(n), d(n)
        real(DP) :: c1, cn, m
        integer(I4B) :: i
        
        write(*,'(A)') 'Test 4: Minimum points case (3 points, sw1=2, sw2=2)'
        
        ! Test data: minimum case with 3 points
        x = [0.0_DP, 1.0_DP, 2.0_DP]
        y = [1.0_DP, 2.0_DP, 1.5_DP]
        indx = [(i, i=1,n)]
        lambda = 1.0_DP  ! Pure fitting
        c1 = 0.0_DP      ! Natural boundaries
        cn = 0.0_DP
        m = 0.0_DP
        
        call splinecof3_a(x, y, c1, cn, lambda, indx, 2, 2, a, b, c, d, m, weight_function)
        
        if (all_finite(a) .and. all_finite(b) .and. &
            all_finite(c) .and. all_finite(d)) then
            write(*,'(A)') '  PASSED: Minimum points case handled successfully'
        else
            write(*,'(A)') '  FAILED: Minimum points case failed'
            all_tests_passed = .false.
        end if
    end subroutine test_minimum_points

    !> Helper function to check if all elements in array are finite
    logical function all_finite(x)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        real(DP), intent(in) :: x(:)
        integer :: i
        all_finite = .true.
        do i = 1, size(x)
            if (.not. ieee_is_finite(x(i))) then
                all_finite = .false.
                return
            end if
        end do
    end function all_finite

end program test_spline_cof_interface