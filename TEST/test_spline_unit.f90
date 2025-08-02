program test_spline_unit
    use nrtype, only: I4B, DP
    implicit none
    
    interface
        subroutine splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
             a, b, c, d, m, f)
            use nrtype, only: I4B, DP
            real(DP),                   intent(inout) :: c1, cn
            real(DP),     dimension(:), intent(in)    :: x, y, lambda1
            integer(I4B), dimension(:), intent(in)    :: indx
            real(DP),     dimension(:), intent(out)   :: a, b, c, d
            integer(I4B),               intent(in)    :: sw1, sw2
            real(DP),                   intent(in)    :: m
            interface
               function f(x,m)
                 use nrtype, only : DP
                 implicit none
                 real(DP), intent(in) :: x, m
                 real(DP)             :: f
               end function f
            end interface
        end subroutine splinecof3_a
    end interface

    ! Test parameters
    real(DP), parameter :: tolerance = 1.0e-12  ! Tolerance for numerical differences between implementations
    logical :: all_tests_passed = .true.
    
    write(*,'(A)') '=== Large Spline Unit Tests ==='
    write(*,'(A)') ''
    
    ! Test with sizes matching performance benchmark
    call test_large_splines(50)
    call test_large_splines(100)
    call test_large_splines(200)
    call test_large_splines(500)
    
    if (all_tests_passed) then
        write(*,'(A)') ''
        write(*,'(A)') 'All large spline tests PASSED!'
        stop 0
    else
        write(*,'(A)') ''
        write(*,'(A)') 'Some large spline tests FAILED!'
        stop 1
    end if

contains

    !> Test function for spline fitting
    function test_function(x, m) result(f_val)
        real(DP), intent(in) :: x, m
        real(DP) :: f_val
        f_val = 1.0_DP
    end function test_function

    !> Test splines with many points
    subroutine test_large_splines(n_intervals)
        integer(I4B), intent(in) :: n_intervals
        integer(I4B) :: n_points, i, j
        real(DP), allocatable :: x(:), y(:), lambda1(:)
        integer(I4B), allocatable :: indx(:)
        real(DP), allocatable :: a(:), b(:), c(:), d(:)
        real(DP) :: c1, cn, m
        integer(I4B) :: sw1, sw2
        real(DP) :: h, x_eval, y_eval, y_spline
        logical :: test_passed
        
        write(*,'(A,I5,A)') 'Testing with ', n_intervals, ' intervals...'
        
        ! Setup problem size
        n_points = n_intervals * 5  ! 5 points per interval on average
        
        ! Allocate arrays
        allocate(x(n_points), y(n_points))
        allocate(indx(n_intervals), lambda1(n_intervals))
        allocate(a(n_intervals), b(n_intervals), c(n_intervals), d(n_intervals))
        
        ! Generate test data - smooth function y = sin(x) + 0.1*x^2
        do i = 1, n_points
            x(i) = real(i-1, DP) * 10.0_DP / real(n_points-1, DP)
            y(i) = sin(x(i)) + 0.1_DP * x(i)**2
        end do
        
        ! Create index array - evenly spaced intervals
        do i = 1, n_intervals
            indx(i) = 1 + (i-1) * (n_points-1) / (n_intervals-1)
        end do
        indx(n_intervals) = n_points  ! Ensure last point is included
        
        ! Test different scenarios
        
        ! Test 1: Natural boundary conditions (fast path)
        write(*,'(A)') '  Test 1: Natural boundaries...'
        lambda1 = 1.0_DP
        c1 = 0.0_DP
        cn = 0.0_DP
        sw1 = 2
        sw2 = 4
        m = 0.0_DP
        
        call splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                          a, b, c, d, m, test_function)
        
        ! Verify spline continuity
        test_passed = .true.
        do i = 1, n_intervals-1
            j = indx(i+1)
            h = x(j) - x(indx(i))
            ! Check continuity at knot points
            y_eval = a(i) + h*(b(i) + h*(c(i) + h*d(i)))
            if (abs(y_eval - a(i+1)) > tolerance) then
                test_passed = .false.
                write(*,'(A,I4,A,E12.4)') '    Continuity error at interval ', i, ': ', abs(y_eval - a(i+1))
            end if
        end do
        
        if (test_passed) then
            write(*,'(A)') '    Natural boundaries: PASSED'
        else
            write(*,'(A)') '    Natural boundaries: FAILED'
            all_tests_passed = .false.
        end if
        
        ! Test 2: Mixed boundary conditions (sparse path)
        write(*,'(A)') '  Test 2: Mixed boundaries...'
        sw1 = 1  ! First derivative
        sw2 = 3  ! Different condition
        c1 = 1.0_DP  ! dy/dx at start
        cn = -1.0_DP  ! dy/dx at end
        
        call splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                          a, b, c, d, m, test_function)
        
        ! Check that first derivative matches at start
        if (abs(b(1) - c1) > tolerance) then
            test_passed = .false.
            write(*,'(A,E12.4)') '    First derivative error: ', abs(b(1) - c1)
        end if
        
        if (test_passed) then
            write(*,'(A)') '    Mixed boundaries: PASSED'
        else
            write(*,'(A)') '    Mixed boundaries: FAILED'
            all_tests_passed = .false.
        end if
        
        ! Test 3: Non-uniform lambda weights
        write(*,'(A)') '  Test 3: Non-uniform weights...'
        sw1 = 2
        sw2 = 4
        c1 = 0.0_DP
        cn = 0.0_DP
        ! Create varying weights
        do i = 1, n_intervals
            lambda1(i) = 0.5_DP + 0.5_DP * sin(real(i, DP))
        end do
        
        call splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                          a, b, c, d, m, test_function)
        
        ! Basic validation - check coefficients are finite
        test_passed = .true.
        do i = 1, n_intervals
            if (.not. (abs(a(i)) < 1.0e10_DP .and. abs(b(i)) < 1.0e10_DP .and. &
                       abs(c(i)) < 1.0e10_DP .and. abs(d(i)) < 1.0e10_DP)) then
                test_passed = .false.
                write(*,'(A,I4)') '    Non-finite coefficients at interval ', i
            end if
        end do
        
        if (test_passed) then
            write(*,'(A)') '    Non-uniform weights: PASSED'
        else
            write(*,'(A)') '    Non-uniform weights: FAILED'
            all_tests_passed = .false.
        end if
        
        ! Test 4: Edge case - many points in single interval
        write(*,'(A)') '  Test 4: Dense data points...'
        ! Reset to have most points in middle intervals
        indx(1) = 1
        do i = 2, n_intervals-1
            indx(i) = n_points/4 + (i-2) * (n_points/2) / max(1, n_intervals-3)
        end do
        indx(n_intervals) = n_points
        
        lambda1 = 1.0_DP
        c1 = 0.0_DP
        cn = 0.0_DP
        sw1 = 2
        sw2 = 4
        
        call splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                          a, b, c, d, m, test_function)
        
        ! Check for reasonable coefficient magnitudes
        test_passed = .true.
        do i = 1, n_intervals
            if (abs(d(i)) > 50000.0_DP) then
                test_passed = .false.
                write(*,'(A,I4,A,E12.4)') '    Extremely large d coefficient at interval ', i, ': ', d(i)
            end if
        end do
        
        if (test_passed) then
            write(*,'(A)') '    Dense data points: PASSED'
        else
            write(*,'(A)') '    Dense data points: FAILED' 
            all_tests_passed = .false.
        end if
        
        ! Cleanup
        deallocate(x, y, indx, lambda1, a, b, c, d)
        
    end subroutine test_large_splines

end program test_spline_unit