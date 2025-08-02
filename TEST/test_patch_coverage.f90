program test_patch_coverage
    !> Simple targeted test to improve patch coverage
    !> Focuses on spline_matrix_assembly.f90 which has 0% coverage
    use nrtype, only: I4B, DP
    use spline_matrix_assembly_mod, only: assemble_spline_matrix_sparse_coo
    ! Focus on matrix assembly only for now
    implicit none
    
    logical :: all_tests_passed = .true.
    
    write(*,'(A)') '=== Patch Coverage Tests ==='
    write(*,'(A)') 'Testing key functions to improve PR diff coverage'
    write(*,'(A)') ''
    
    ! Test 1: Exercise spline_matrix_assembly.f90 (0% coverage)
    call test_matrix_assembly()
    
    if (all_tests_passed) then
        write(*,'(A)') ''
        write(*,'(A)') 'All patch coverage tests PASSED!'
        stop 0
    else
        write(*,'(A)') ''
        write(*,'(A)') 'Some patch coverage tests FAILED!'
        stop 1
    end if
    
contains

    function test_function(x, m) result(f_val)
        real(DP), intent(in) :: x, m
        real(DP) :: f_val
        f_val = 1.0_DP + 0.1_DP*x  ! Simple linear + constant
    end function test_function

    subroutine test_matrix_assembly()
        integer(I4B), parameter :: n = 4
        real(DP) :: x(n), y(n)
        real(DP) :: c1, cn, m
        real(DP) :: lambda1(n)
        integer(I4B) :: indx(n)
        integer(I4B) :: nrow, ncol, nnz
        integer(I4B), allocatable :: irow_coo(:), icol_coo(:)
        real(DP), allocatable :: val_coo(:), rhs(:)
        integer(I4B) :: i
        
        write(*,'(A)') 'Test 1: Matrix assembly functions (spline_matrix_assembly.f90)'
        
        ! Setup test data
        x = [0.0_DP, 1.0_DP, 2.0_DP, 3.0_DP]
        y = [1.0_DP, 1.2_DP, 1.5_DP, 1.8_DP]
        indx = [(i, i=1,n)]
        lambda1 = 0.8_DP
        c1 = 0.0_DP
        cn = 0.0_DP
        m = 0.0_DP
        
        call assemble_spline_matrix_sparse_coo(x, y, c1, cn, lambda1, indx, 2, 4, &
                                              m, test_function, nrow, ncol, nnz, &
                                              irow_coo, icol_coo, val_coo, rhs)
        
        if (allocated(irow_coo) .and. nrow > 0 .and. nnz > 0) then
            write(*,'(A,I0,A,I0,A,I0)') '  PASSED: Matrix assembled (', nrow, 'x', ncol, ', nnz=', nnz, ')'
            deallocate(irow_coo, icol_coo, val_coo, rhs)
        else
            write(*,'(A)') '  FAILED: Matrix assembly failed'
            all_tests_passed = .false.
        end if
    end subroutine test_matrix_assembly

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

end program test_patch_coverage