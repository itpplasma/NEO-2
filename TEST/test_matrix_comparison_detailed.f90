program test_matrix_comparison_detailed
    ! Detailed comparison of matrix elements from dense vs sparse implementations
    use nrtype, only: I4B, DP
    use sparse_mod, only: full2sparse, sparse2full
    use spline_matrix_assembly_mod, only: assemble_spline_matrix_sparse_coo
    implicit none
    
    ! Test parameters
    integer(I4B), parameter :: n = 5  ! Small size for detailed analysis
    real(DP) :: x(n), y(n)
    integer(I4B) :: indx(n)
    real(DP) :: lambda1(n)
    real(DP) :: c1, cn, m
    integer(I4B) :: sw1, sw2, i, j
    
    ! Dense matrix from original
    real(DP), allocatable :: MA_dense(:,:), MA_reconstructed(:,:)
    real(DP), allocatable :: rhs_dense(:)
    
    ! Sparse matrix from direct implementation  
    integer(I4B) :: nrow, ncol, nnz_direct
    integer(I4B), allocatable :: irow_direct(:), icol_direct(:)
    real(DP), allocatable :: val_direct(:), rhs_direct(:)
    
    ! Sparse matrix from dense conversion
    integer(I4B) :: nnz_converted
    integer(I4B), allocatable :: irow_converted(:), pcol_converted(:)
    real(DP), allocatable :: val_converted(:)
    
    ! Interface for test function
    interface
        function test_function(x, m) result(f_val)
            use nrtype, only : DP
            implicit none
            real(DP), intent(in) :: x, m
            real(DP) :: f_val
        end function test_function
    end interface
    
    interface
        subroutine splinecof3_original_dense(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                                             a, b, c, d, m, f)
            use nrtype, only : I4B, DP
            implicit none
            real(DP), dimension(:), intent(in) :: x, y, lambda1
            real(DP), intent(inout) :: c1, cn
            integer(I4B), dimension(:), intent(in) :: indx
            real(DP), dimension(:), intent(out) :: a, b, c, d
            integer(I4B), intent(in) :: sw1, sw2
            real(DP), intent(in) :: m
            interface
                function f(x,m)
                    use nrtype, only : DP
                    implicit none
                    real(DP), intent(in) :: x, m
                    real(DP) :: f
                end function f
            end interface
        end subroutine splinecof3_original_dense
    end interface
    
    write(*,'(A)') '=== Detailed Matrix Comparison ==='
    write(*,'(A)') ''
    
    ! Setup test data
    do i = 1, n
        x(i) = real(i-1, DP) * 0.5_DP
    end do
    y = x**2  ! Simple quadratic
    indx = [(i, i=1,n)]
    lambda1 = 1.0_DP
    m = 0.0_DP
    c1 = 0.0_DP
    cn = 0.0_DP
    sw1 = 2  ! Natural
    sw2 = 4  ! Natural
    
    write(*,'(A)') 'Test configuration:'
    write(*,'(A,I0)') '  Number of points: ', n
    write(*,'(A)') '  Boundary conditions: Natural (sw1=2, sw2=4)'
    write(*,'(A)') ''
    
    ! Get the dense matrix from original implementation
    call get_dense_matrix(x, y, c1, cn, lambda1, indx, sw1, sw2, m, &
                         MA_dense, rhs_dense)
    
    nrow = size(MA_dense, 1)
    ncol = size(MA_dense, 2)
    
    write(*,'(A,I0,A,I0)') 'Dense matrix size: ', nrow, ' x ', ncol
    
    ! Convert dense to sparse
    call full2sparse(MA_dense, irow_converted, pcol_converted, val_converted, nrow, ncol, nnz_converted)
    write(*,'(A,I0)') 'Non-zeros after conversion: ', nnz_converted
    
    ! Get sparse matrix from direct implementation
    call assemble_spline_matrix_sparse_coo(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                                          m, test_function, nrow, ncol, nnz_direct, &
                                          irow_direct, icol_direct, val_direct, rhs_direct)
    write(*,'(A,I0)') 'Non-zeros from direct sparse: ', nnz_direct
    write(*,'(A)') ''
    
    ! Compare number of non-zeros
    if (nnz_converted /= nnz_direct) then
        write(*,'(A)') '✗ Different number of non-zeros!'
        write(*,'(A,I0,A,I0)') '  Dense->sparse: ', nnz_converted, ', Direct sparse: ', nnz_direct
        
        ! Show which elements are different
        call show_matrix_differences(nrow, ncol, nnz_converted, irow_converted, pcol_converted, val_converted, &
                                    nnz_direct, irow_direct, icol_direct, val_direct)
    else
        write(*,'(A)') '✓ Same number of non-zeros'
        
        ! Check if elements match
        call compare_matrix_elements(nnz_converted, irow_converted, pcol_converted, val_converted, &
                                    nnz_direct, irow_direct, icol_direct, val_direct)
    end if
    
    ! Show first few rows of dense matrix for inspection
    write(*,'(A)') ''
    write(*,'(A)') 'First 5x5 block of dense matrix:'
    do i = 1, min(5, nrow)
        write(*,'(I3,A)', advance='no') i, ': '
        do j = 1, min(5, ncol)
            if (abs(MA_dense(i,j)) > 1e-15) then
                write(*,'(F10.6)', advance='no') MA_dense(i,j)
            else
                write(*,'(A10)', advance='no') '    0     '
            end if
        end do
        write(*,*)
    end do
    
    ! Clean up
    deallocate(MA_dense, rhs_dense)
    deallocate(irow_converted, pcol_converted, val_converted)
    deallocate(irow_direct, icol_direct, val_direct, rhs_direct)
    
contains

    subroutine get_dense_matrix(x, y, c1, cn, lambda1, indx, sw1, sw2, m, MA, rhs)
        real(DP), dimension(:), intent(in) :: x, y, lambda1
        real(DP), intent(in) :: c1, cn, m
        integer(I4B), dimension(:), intent(in) :: indx
        integer(I4B), intent(in) :: sw1, sw2
        real(DP), allocatable, intent(out) :: MA(:,:), rhs(:)
        
        ! This is a hack to extract the matrix from the original implementation
        ! We'll call it but interrupt before solving
        real(DP), dimension(size(x)) :: a_dummy, b_dummy, c_dummy, d_dummy
        real(DP) :: c1_local, cn_local
        
        ! For now, just create a dummy matrix to show the concept
        integer(I4B) :: VAR = 7
        integer(I4B) :: size_dimension
        
        size_dimension = VAR * size(indx) - 2
        allocate(MA(size_dimension, size_dimension))
        allocate(rhs(size_dimension))
        
        MA = 0.0_DP
        rhs = 0.0_DP
        
        ! In practice, we'd need to modify splinecof3_original_dense to export MA
        write(*,'(A)') 'Note: Matrix extraction from original would require code modification'
        
    end subroutine get_dense_matrix
    
    subroutine show_matrix_differences(nrow, ncol, nnz1, irow1, pcol1, val1, &
                                      nnz2, irow2, icol2, val2)
        integer(I4B), intent(in) :: nrow, ncol, nnz1, nnz2
        integer(I4B), dimension(:), intent(in) :: irow1, pcol1, irow2, icol2
        real(DP), dimension(:), intent(in) :: val1, val2
        
        integer(I4B) :: i, j, found
        logical :: in_first, in_second
        
        write(*,'(A)') ''
        write(*,'(A)') 'Elements only in dense->sparse conversion:'
        do i = 1, nnz1
            found = 0
            do j = 1, nnz2
                ! Note: pcol1 is column pointer format, icol2 is direct column indices
                ! This comparison would need proper conversion
            end do
        end do
        
        write(*,'(A)') 'Analysis requires proper CSC to COO conversion'
        
    end subroutine show_matrix_differences
    
    subroutine compare_matrix_elements(nnz1, irow1, pcol1, val1, &
                                      nnz2, irow2, icol2, val2)
        integer(I4B), intent(in) :: nnz1, nnz2
        integer(I4B), dimension(:), intent(in) :: irow1, pcol1, irow2, icol2
        real(DP), dimension(:), intent(in) :: val1, val2
        
        write(*,'(A)') 'Element-by-element comparison requires format conversion'
        
    end subroutine compare_matrix_elements
    
end program test_matrix_comparison_detailed

! Test function implementation
function test_function(x, m) result(f_val)
    use nrtype, only : DP
    implicit none
    real(DP), intent(in) :: x, m
    real(DP) :: f_val
    f_val = 1.0_DP
end function test_function