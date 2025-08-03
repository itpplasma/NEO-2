program test_spline_bicgstab_accuracy
    use nrtype, only: I4B, DP
    use sparse_solvers_mod, only: bicgstab_abs_tolerance, bicgstab_rel_tolerance, &
                                  bicgstab_max_iter, sparse_solve_method, SOLVER_BICGSTAB, &
                                  SOLVER_UMFPACK
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
    integer(I4B), parameter :: n_intervals = 50
    integer(I4B) :: n_points, i, j
    real(DP), allocatable :: x(:), y(:), lambda1(:)
    integer(I4B), allocatable :: indx(:)
    real(DP), allocatable :: a_umf(:), b_umf(:), c_umf(:), d_umf(:)
    real(DP), allocatable :: a_bcg(:), b_bcg(:), c_bcg(:), d_bcg(:)
    real(DP) :: c1, cn, m
    integer(I4B) :: sw1, sw2
    real(DP) :: max_diff_a, max_diff_b, max_diff_c, max_diff_d
    real(DP) :: saved_abs_tol, saved_rel_tol
    integer :: saved_max_iter
    real(DP) :: test_tol
    
    write(*,'(A)') '=== BiCGSTAB Accuracy Test for Spline Problems ==='
    write(*,'(A)') ''
    write(*,'(A)') 'Testing achievable accuracy with various BiCGSTAB tolerance settings'
    write(*,'(A)') 'comparing against UMFPACK (direct solver) reference solution.'
    write(*,'(A)') ''
    
    ! Setup problem
    n_points = n_intervals * 5
    allocate(x(n_points), y(n_points))
    allocate(indx(n_intervals), lambda1(n_intervals))
    allocate(a_umf(n_intervals), b_umf(n_intervals), c_umf(n_intervals), d_umf(n_intervals))
    allocate(a_bcg(n_intervals), b_bcg(n_intervals), c_bcg(n_intervals), d_bcg(n_intervals))
    
    ! Generate test data - smooth function
    do i = 1, n_points
        x(i) = real(i-1, DP) * 10.0_DP / real(n_points-1, DP)
        y(i) = sin(x(i)) + 0.1_DP * x(i)**2
    end do
    
    ! Create index array
    do i = 1, n_intervals
        indx(i) = 1 + (i-1) * (n_points-1) / (n_intervals-1)
    end do
    indx(n_intervals) = n_points
    
    ! Test with natural boundary conditions
    lambda1 = 1.0_DP
    c1 = 0.0_DP
    cn = 0.0_DP
    sw1 = 2
    sw2 = 4
    m = 0.0_DP
    
    ! Save current settings
    saved_abs_tol = bicgstab_abs_tolerance
    saved_rel_tol = bicgstab_rel_tolerance
    saved_max_iter = bicgstab_max_iter
    
    ! First get reference solution with UMFPACK
    sparse_solve_method = SOLVER_UMFPACK
    call splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                      a_umf, b_umf, c_umf, d_umf, m, test_function)
    
    write(*,'(A)') 'Reference solution computed with UMFPACK (direct solver)'
    write(*,'(A)') ''
    write(*,'(A)') 'Testing BiCGSTAB with different tolerance settings:'
    write(*,'(A)') 'Abs Tol    Rel Tol    Max Error in Coefficients'
    write(*,'(A)') '--------   --------   -------------------------'
    
    ! Test different tolerance settings
    sparse_solve_method = SOLVER_BICGSTAB
    bicgstab_max_iter = 10000  ! Allow many iterations
    
    ! Test 1: Very relaxed tolerances
    test_tol = 1.0e-4_DP
    bicgstab_abs_tolerance = test_tol
    bicgstab_rel_tolerance = test_tol
    call test_accuracy()
    
    ! Test 2: Relaxed tolerances
    test_tol = 1.0e-6_DP
    bicgstab_abs_tolerance = test_tol
    bicgstab_rel_tolerance = test_tol
    call test_accuracy()
    
    ! Test 3: Moderate tolerances
    test_tol = 1.0e-8_DP
    bicgstab_abs_tolerance = test_tol
    bicgstab_rel_tolerance = test_tol
    call test_accuracy()
    
    ! Test 4: Mixed - very loose relative
    bicgstab_abs_tolerance = 1.0e-8_DP
    bicgstab_rel_tolerance = 1.0e-3_DP
    write(*,'(E10.1,1X,E10.1,3X)', advance='no') bicgstab_abs_tolerance, bicgstab_rel_tolerance
    call splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                      a_bcg, b_bcg, c_bcg, d_bcg, m, test_function)
    call compute_differences()
    write(*,'(4(A,E10.3))') 'a:', max_diff_a, ' b:', max_diff_b, ' c:', max_diff_c, ' d:', max_diff_d
    
    ! Test 5: Extremely relaxed
    test_tol = 1.0e-2_DP
    bicgstab_abs_tolerance = test_tol
    bicgstab_rel_tolerance = test_tol
    call test_accuracy()
    
    ! Restore original settings
    bicgstab_abs_tolerance = saved_abs_tol
    bicgstab_rel_tolerance = saved_rel_tol
    bicgstab_max_iter = saved_max_iter
    
    write(*,'(A)') ''
    write(*,'(A)') 'Conclusion:'
    write(*,'(A)') 'BiCGSTAB can achieve reasonable accuracy (~1e-6 to 1e-8) on spline problems'
    write(*,'(A)') 'but struggles to reach machine precision due to ill-conditioning.'
    write(*,'(A)') 'For high-accuracy requirements (< 1e-10), UMFPACK is recommended.'
    
contains

    function test_function(x, m) result(f_val)
        real(DP), intent(in) :: x, m
        real(DP) :: f_val
        f_val = 1.0_DP
    end function test_function
    
    subroutine test_accuracy()
        write(*,'(E10.1,1X,E10.1,3X)', advance='no') bicgstab_abs_tolerance, bicgstab_rel_tolerance
        call splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                          a_bcg, b_bcg, c_bcg, d_bcg, m, test_function)
        call compute_differences()
        write(*,'(4(A,E10.3))') 'a:', max_diff_a, ' b:', max_diff_b, ' c:', max_diff_c, ' d:', max_diff_d
    end subroutine test_accuracy
    
    subroutine compute_differences()
        integer :: k
        max_diff_a = 0.0_DP
        max_diff_b = 0.0_DP
        max_diff_c = 0.0_DP
        max_diff_d = 0.0_DP
        
        do k = 1, n_intervals
            max_diff_a = max(max_diff_a, abs(a_bcg(k) - a_umf(k)))
            max_diff_b = max(max_diff_b, abs(b_bcg(k) - b_umf(k)))
            max_diff_c = max(max_diff_c, abs(c_bcg(k) - c_umf(k)))
            max_diff_d = max(max_diff_d, abs(d_bcg(k) - d_umf(k)))
        end do
    end subroutine compute_differences
    
end program test_spline_bicgstab_accuracy