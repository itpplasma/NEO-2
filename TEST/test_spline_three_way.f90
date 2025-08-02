program test_spline_three_way
    use nrtype, only: I4B, DP
    use splinecof3_fast_mod, only: splinecof3_general_fast
    use spline_test_control, only: disable_fast_path
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
        
        subroutine splinecof3_original_dense(x, y, c1, cn, lambda1, indx, sw1, sw2, &
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
        end subroutine splinecof3_original_dense
    end interface
    
    real(DP), parameter :: tolerance = 1.0e-11
    logical :: all_tests_passed = .true.
    
    write(*,'(A)') '=== Three-Way Spline Implementation Comparison ==='
    write(*,'(A)') ''
    
    ! Test 1: True fast path case (consecutive indices, natural BC)
    call test_fast_path_natural()
    
    ! Test 2: Non-consecutive indices (sparse path only)
    call test_sparse_path_natural()
    
    ! Test 3: Different boundary conditions
    call test_sparse_path_mixed()
    
    ! Test 4: Force sparse path for fast-path-eligible case
    call test_forced_sparse_path()
    
    if (all_tests_passed) then
        write(*,'(A)') ''
        write(*,'(A)') 'All tests PASSED!'
        stop 0
    else
        write(*,'(A)') ''
        write(*,'(A)') 'Some tests FAILED!'
        stop 1
    end if
    
contains

    function test_function(x, m) result(f_val)
        real(DP), intent(in) :: x, m
        real(DP) :: f_val
        f_val = 1.0_DP
    end function test_function
    
    subroutine test_fast_path_natural()
        integer(I4B), parameter :: n = 5
        real(DP) :: x(n), y(n)
        integer(I4B) :: indx(n)
        real(DP) :: lambda1(n)
        real(DP) :: a_sparse(n), b_sparse(n), c_sparse(n), d_sparse(n)
        real(DP) :: a_fast(n), b_fast(n), c_fast(n), d_fast(n)
        real(DP) :: a_orig(n), b_orig(n), c_orig(n), d_orig(n)
        real(DP) :: c1, cn, m
        integer(I4B) :: sw1, sw2, i
        logical :: test_passed
        
        write(*,'(A)') 'Test 1: Fast path eligible (consecutive indices, natural BC)'
        
        ! Setup data
        x = [0.0_DP, 1.0_DP, 2.0_DP, 3.0_DP, 4.0_DP]
        y = [0.0_DP, 1.0_DP, 4.0_DP, 9.0_DP, 16.0_DP]  ! x^2
        indx = [(i, i=1,n)]  ! Consecutive indices for fast path
        lambda1 = 1.0_DP
        c1 = 0.0_DP
        cn = 0.0_DP
        sw1 = 2  ! Natural
        sw2 = 4  ! Natural
        m = 0.0_DP
        
        ! Original dense
        call splinecof3_original_dense(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                                      a_orig, b_orig, c_orig, d_orig, m, test_function)
        
        ! Fast implementation
        c1 = 0.0_DP; cn = 0.0_DP
        call splinecof3_general_fast(x, y, c1, cn, sw1, sw2, a_fast, b_fast, c_fast, d_fast)
        
        ! Wrapper with fast path enabled (should use fast path)
        c1 = 0.0_DP; cn = 0.0_DP
        disable_fast_path = .false.
        call splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                          a_sparse, b_sparse, c_sparse, d_sparse, m, test_function)
        
        ! Compare all three
        test_passed = .true.
        
        ! Compare fast vs original
        if (any(abs(a_fast(1:n-1) - a_orig(1:n-1)) > tolerance)) then
            write(*,'(A)') '  FAILED: Fast vs Original - a coefficients differ'
            test_passed = .false.
        end if
        
        ! Compare sparse vs original
        if (any(abs(a_sparse(1:n-1) - a_orig(1:n-1)) > tolerance)) then
            write(*,'(A)') '  FAILED: Sparse vs Original - a coefficients differ'
            write(*,'(A,5E12.5)') '    Original a:', a_orig(1:n-1)
            write(*,'(A,5E12.5)') '    Sparse a:  ', a_sparse(1:n-1)
            write(*,'(A,5E12.5)') '    Diff:      ', abs(a_sparse(1:n-1) - a_orig(1:n-1))
            test_passed = .false.
        end if
        
        if (test_passed) then
            write(*,'(A)') '  PASSED: All three implementations agree'
        end if
        
        if (.not. test_passed) all_tests_passed = .false.
        
    end subroutine test_fast_path_natural
    
    subroutine test_sparse_path_natural()
        integer(I4B), parameter :: n = 5
        real(DP) :: x(n), y(n)
        integer(I4B) :: indx(3)
        real(DP) :: lambda1(3)
        real(DP) :: a_sparse(3), b_sparse(3), c_sparse(3), d_sparse(3)
        real(DP) :: a_orig(3), b_orig(3), c_orig(3), d_orig(3)
        real(DP) :: c1, cn, m
        integer(I4B) :: sw1, sw2
        logical :: test_passed
        
        write(*,'(A)') ''
        write(*,'(A)') 'Test 2: Sparse path only (non-consecutive indices)'
        
        ! Setup data
        x = [0.0_DP, 1.0_DP, 2.0_DP, 3.0_DP, 4.0_DP]
        y = [0.0_DP, 1.0_DP, 4.0_DP, 9.0_DP, 16.0_DP]  ! x^2
        indx = [1, 3, 5]  ! Non-consecutive - forces sparse path
        lambda1 = 1.0_DP
        c1 = 0.0_DP
        cn = 0.0_DP
        sw1 = 2  ! Natural
        sw2 = 4  ! Natural
        m = 0.0_DP
        
        ! Original dense
        call splinecof3_original_dense(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                                      a_orig, b_orig, c_orig, d_orig, m, test_function)
        
        ! Wrapper (should use sparse path)
        c1 = 0.0_DP; cn = 0.0_DP
        call splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                          a_sparse, b_sparse, c_sparse, d_sparse, m, test_function)
        
        ! Compare
        test_passed = .true.
        
        if (any(abs(a_sparse - a_orig) > tolerance) .or. &
            any(abs(b_sparse - b_orig) > tolerance) .or. &
            any(abs(c_sparse - c_orig) > tolerance) .or. &
            any(abs(d_sparse - d_orig) > tolerance)) then
            write(*,'(A)') '  FAILED: Sparse vs Original differ'
            write(*,'(A,3E12.5)') '    a diff:', abs(a_sparse - a_orig)
            write(*,'(A,3E12.5)') '    b diff:', abs(b_sparse - b_orig)
            write(*,'(A,3E12.5)') '    c diff:', abs(c_sparse - c_orig)
            write(*,'(A,3E12.5)') '    d diff:', abs(d_sparse - d_orig)
            test_passed = .false.
        else
            write(*,'(A)') '  PASSED: Sparse and Original agree'
        end if
        
        if (.not. test_passed) all_tests_passed = .false.
        
    end subroutine test_sparse_path_natural
    
    subroutine test_sparse_path_mixed()
        integer(I4B), parameter :: n = 5
        real(DP) :: x(n), y(n)
        integer(I4B) :: indx(3)
        real(DP) :: lambda1(3)
        real(DP) :: a_sparse(3), b_sparse(3), c_sparse(3), d_sparse(3)
        real(DP) :: a_orig(3), b_orig(3), c_orig(3), d_orig(3)
        real(DP) :: c1, cn, m
        integer(I4B) :: sw1, sw2
        logical :: test_passed
        
        write(*,'(A)') ''
        write(*,'(A)') 'Test 3: Mixed boundary conditions'
        
        ! Setup data
        x = [0.0_DP, 1.0_DP, 2.0_DP, 3.0_DP, 4.0_DP]
        y = [0.0_DP, 1.0_DP, 4.0_DP, 9.0_DP, 16.0_DP]  ! x^2
        indx = [1, 3, 5]
        lambda1 = 1.0_DP
        c1 = 2.0_DP  ! First derivative
        cn = 8.0_DP  ! Last derivative
        sw1 = 1  ! First derivative
        sw2 = 3  ! Last derivative
        m = 0.0_DP
        
        ! Original dense
        call splinecof3_original_dense(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                                      a_orig, b_orig, c_orig, d_orig, m, test_function)
        
        ! Wrapper (should use sparse path)
        c1 = 2.0_DP; cn = 8.0_DP
        call splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                          a_sparse, b_sparse, c_sparse, d_sparse, m, test_function)
        
        ! Compare
        test_passed = .true.
        
        if (any(abs(a_sparse - a_orig) > tolerance) .or. &
            any(abs(b_sparse - b_orig) > tolerance) .or. &
            any(abs(c_sparse - c_orig) > tolerance) .or. &
            any(abs(d_sparse - d_orig) > tolerance)) then
            write(*,'(A)') '  FAILED: Sparse vs Original differ'
            write(*,'(A,3E12.5)') '    a diff:', abs(a_sparse - a_orig)
            write(*,'(A,3E12.5)') '    b diff:', abs(b_sparse - b_orig)
            write(*,'(A,3E12.5)') '    c diff:', abs(c_sparse - c_orig)
            write(*,'(A,3E12.5)') '    d diff:', abs(d_sparse - d_orig)
            write(*,'(A)') '  Debug: b values'
            write(*,'(A,3F10.6)') '    Original b:', b_orig
            write(*,'(A,3F10.6)') '    Sparse b:  ', b_sparse
            write(*,'(A,F10.6)') '    cn value: ', cn
            test_passed = .false.
        else
            write(*,'(A)') '  PASSED: Sparse and Original agree'
        end if
        
        if (.not. test_passed) all_tests_passed = .false.
        
    end subroutine test_sparse_path_mixed
    
    subroutine test_forced_sparse_path()
        integer(I4B), parameter :: n = 5
        real(DP) :: x(n), y(n)
        integer(I4B) :: indx(n)
        real(DP) :: lambda1(n)
        real(DP) :: a_sparse(n), b_sparse(n), c_sparse(n), d_sparse(n)
        real(DP) :: a_forced(n), b_forced(n), c_forced(n), d_forced(n)
        real(DP) :: a_orig(n), b_orig(n), c_orig(n), d_orig(n)
        real(DP) :: c1, cn, m
        integer(I4B) :: sw1, sw2, i
        logical :: test_passed
        
        write(*,'(A)') ''
        write(*,'(A)') 'Test 4: Force sparse path for fast-path-eligible case'
        
        ! Setup data (fast path eligible)
        x = [0.0_DP, 1.0_DP, 2.0_DP, 3.0_DP, 4.0_DP]
        y = [0.0_DP, 1.0_DP, 4.0_DP, 9.0_DP, 16.0_DP]  ! x^2
        indx = [(i, i=1,n)]  ! Consecutive indices
        lambda1 = 1.0_DP
        c1 = 0.0_DP
        cn = 0.0_DP
        sw1 = 2  ! Natural
        sw2 = 4  ! Natural
        m = 0.0_DP
        
        ! Original dense
        call splinecof3_original_dense(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                                      a_orig, b_orig, c_orig, d_orig, m, test_function)
        
        ! Normal call (should use fast path)
        c1 = 0.0_DP; cn = 0.0_DP
        disable_fast_path = .false.
        call splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                          a_sparse, b_sparse, c_sparse, d_sparse, m, test_function)
        
        ! Forced sparse path
        c1 = 0.0_DP; cn = 0.0_DP
        disable_fast_path = .true.
        call splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                          a_forced, b_forced, c_forced, d_forced, m, test_function)
        
        ! Reset flag
        disable_fast_path = .false.
        
        ! Compare all three
        test_passed = .true.
        
        ! Compare normal vs original
        if (any(abs(a_sparse(1:n-1) - a_orig(1:n-1)) > tolerance)) then
            write(*,'(A)') '  FAILED: Normal (fast path) vs Original differ'
            test_passed = .false.
        end if
        
        ! Compare forced sparse vs original
        if (any(abs(a_forced(1:n-1) - a_orig(1:n-1)) > tolerance) .or. &
            any(abs(b_forced(1:n-1) - b_orig(1:n-1)) > tolerance) .or. &
            any(abs(c_forced(1:n-1) - c_orig(1:n-1)) > tolerance) .or. &
            any(abs(d_forced(1:n-1) - d_orig(1:n-1)) > tolerance)) then
            write(*,'(A)') '  FAILED: Forced sparse vs Original differ'
            write(*,'(A,5E12.5)') '    a diff:', abs(a_forced(1:n-1) - a_orig(1:n-1))
            write(*,'(A,5E12.5)') '    b diff:', abs(b_forced(1:n-1) - b_orig(1:n-1))
            write(*,'(A,5E12.5)') '    c diff:', abs(c_forced(1:n-1) - c_orig(1:n-1))
            write(*,'(A,5E12.5)') '    d diff:', abs(d_forced(1:n-1) - d_orig(1:n-1))
            test_passed = .false.
        else
            write(*,'(A)') '  PASSED: Forced sparse path matches original'
        end if
        
        if (test_passed) then
            write(*,'(A)') '  All paths produce identical results'
        end if
        
        if (.not. test_passed) all_tests_passed = .false.
        
    end subroutine test_forced_sparse_path

end program test_spline_three_way