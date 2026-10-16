! Native provider for the existing sparse_mod interface. NEO-2 owns the
! factor/solve/free policy; FortSparse owns sparse mathematics and refinement.
module neo2_native_sparse
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, csc_t, csc_z_t, csc_from_triplet, sparse_solver_t, &
        fortsparse_status_t, FORTSPARSE_OK, FORTSPARSE_BACKEND_NATIVE, &
        FORTSPARSE_SYM_GENERAL, sparse_factor, sparse_refactor, sparse_solve, &
        sparse_destroy
    implicit none
    private
    public :: native_sparse_solve
    type(sparse_solver_t), save :: solver
    integer, save :: scalar_kind = 0, factored_method = 0
    interface native_sparse_solve
        module procedure real_vector, complex_vector, real_block, complex_block
    end interface
contains
    subroutine fail(message)
        character(*), intent(in) :: message
        call release()
        write(error_unit,'(a)') 'NEO-2 native sparse: '//message
        error stop 1
    end subroutine

    subroutine checked(status)
        type(fortsparse_status_t), intent(in) :: status
        if (status%code /= FORTSPARSE_OK) then
            write(error_unit,'(a,i0)') 'FortSparse status=',status%code
            call fail(trim(status%msg))
        end if
    end subroutine

    subroutine release()
        call sparse_destroy(solver)
        scalar_kind = 0
        factored_method = 0
    end subroutine

    integer function operation(method, iopt) result(op)
        integer, intent(in) :: method
        integer, optional, intent(in) :: iopt
        if (method /= 2 .and. method /= 3) call fail('only methods 2/3 are supported')
        op = 0
        if (present(iopt)) op = iopt
        if (op < 0 .or. op > 3) call fail('invalid factor/solve/free operation')
    end function

    subroutine columns(nrow, ncol, nz, rows, pointer, result)
        integer, intent(in) :: nrow, ncol, nz, rows(:), pointer(:)
        integer, allocatable, intent(out) :: result(:)
        integer :: j
        logical :: compressed
        if (nrow < 1 .or. ncol /= nrow .or. nz < 0) call fail('invalid dimensions')
        if (ncol == huge(ncol) .or. nz == huge(nz)) call fail('index overflow')
        if (size(rows) /= nz) call fail('row count differs from nz')
        compressed = .false.
        if (size(pointer) == ncol+1) then
            if (pointer(1) == 1 .and. pointer(ncol+1) == nz+1) then
                compressed = all(pointer(2:) >= pointer(:ncol))
            end if
        end if
        allocate(result(nz))
        if (compressed) then
            do j = 1, ncol
                result(pointer(j):pointer(j+1)-1) = j
            end do
        else
            if (size(pointer) /= nz) call fail('invalid CSC/COO column array')
            result = pointer
        end if
        ! FortSparse validates row/column bounds and canonicalizes COO, including
        ! duplicates, unsorted input and trailing empty columns. A valid CSC
        ! pointer ends at nz+1 and cannot be confused with valid COO when
        ! nz=ncol+1 (COO indices are <=ncol).
    end subroutine

    subroutine prepare_real(nrow,ncol,nz,rows,pointer,values,method,verbose,op)
        integer, intent(in) :: nrow,ncol,nz,rows(:),pointer(:),method,op
        real(dp), intent(in) :: values(:)
        logical, intent(in) :: verbose
        type(csc_t) :: matrix
        type(fortsparse_status_t) :: status
        integer, allocatable :: cols(:)
        logical :: fresh
        fresh = op == 0 .or. op == 1 .or. .not. solver%factored
        if (.not. fresh) then
            if (scalar_kind /= 1) call fail('retained factor scalar kind mismatch')
            if (method == factored_method) return
        end if
        if (size(values) /= nz) call fail('value count differs from nz')
        call columns(nrow,ncol,nz,rows,pointer,cols)
        call csc_from_triplet(nrow,ncol,rows,cols,values,matrix,status)
        call checked(status)
        if (fresh) call release()
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = FORTSPARSE_SYM_GENERAL
        solver%options%verbose = verbose
        solver%refine = method == 2
        if (fresh) then
            call sparse_factor(solver,matrix,status)
        else
            ! Refinement is retained at factor time. A method change uses the
            ! library's refactor contract; no refinement formula is copied here.
            call sparse_refactor(solver,matrix,status)
        end if
        call checked(status)
        scalar_kind = 1
        factored_method = method
    end subroutine

    subroutine prepare_complex(nrow,ncol,nz,rows,pointer,values,method,verbose,op)
        integer, intent(in) :: nrow,ncol,nz,rows(:),pointer(:),method,op
        complex(dp), intent(in) :: values(:)
        logical, intent(in) :: verbose
        type(csc_z_t) :: matrix
        type(fortsparse_status_t) :: status
        integer, allocatable :: cols(:)
        logical :: fresh
        fresh = op == 0 .or. op == 1 .or. .not. solver%factored
        if (.not. fresh) then
            if (scalar_kind /= 2) call fail('retained factor scalar kind mismatch')
            if (method == factored_method) return
        end if
        if (size(values) /= nz) call fail('value count differs from nz')
        call columns(nrow,ncol,nz,rows,pointer,cols)
        call csc_from_triplet(nrow,ncol,rows,cols,values,matrix,status)
        call checked(status)
        if (fresh) call release()
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = FORTSPARSE_SYM_GENERAL
        solver%options%verbose = verbose
        solver%refine = method == 2
        if (fresh) then
            call sparse_factor(solver,matrix,status)
        else
            call sparse_refactor(solver,matrix,status)
        end if
        call checked(status)
        scalar_kind = 2
        factored_method = method
    end subroutine

    subroutine real_vector(nrow,ncol,nz,rows,pointer,values,b,method,verbose,iopt)
        integer, intent(in) :: nrow,ncol,nz,rows(:),pointer(:),method
        real(dp), intent(in) :: values(:)
        real(dp), intent(inout) :: b(:)
        logical, intent(in) :: verbose
        integer, optional, intent(in) :: iopt
        type(fortsparse_status_t) :: status
        integer :: op
        op = operation(method,iopt)
        if (op == 3) then
            call release()
            return
        end if
        if (size(b) /= nrow) call fail('incompatible RHS dimensions')
        call prepare_real(nrow,ncol,nz,rows,pointer,values,method,verbose,op)
        if (op == 0 .or. op == 2) then
            call sparse_solve(solver,b,status)
            call checked(status)
        end if
        if (op == 0) call release()
    end subroutine

    subroutine complex_vector(nrow,ncol,nz,rows,pointer,values,b,method,verbose,iopt)
        integer, intent(in) :: nrow,ncol,nz,rows(:),pointer(:),method
        complex(dp), intent(in) :: values(:)
        complex(dp), intent(inout) :: b(:)
        logical, intent(in) :: verbose
        integer, optional, intent(in) :: iopt
        type(fortsparse_status_t) :: status
        integer :: op
        op = operation(method,iopt)
        if (op == 3) then
            call release()
            return
        end if
        if (size(b) /= nrow) call fail('incompatible RHS dimensions')
        call prepare_complex(nrow,ncol,nz,rows,pointer,values,method,verbose,op)
        if (op == 0 .or. op == 2) then
            call sparse_solve(solver,b,status)
            call checked(status)
        end if
        if (op == 0) call release()
    end subroutine

    subroutine real_block(nrow,ncol,nz,rows,pointer,values,b,method,verbose,iopt)
        integer, intent(in) :: nrow,ncol,nz,rows(:),pointer(:),method
        real(dp), intent(in) :: values(:)
        real(dp), intent(inout) :: b(:,:)
        logical, intent(in) :: verbose
        integer, optional, intent(in) :: iopt
        type(fortsparse_status_t) :: status
        integer :: op,j
        op = operation(method,iopt)
        if (op == 3) then
            call release()
            return
        end if
        if (size(b,1) /= nrow) call fail('incompatible block RHS dimensions')
        call prepare_real(nrow,ncol,nz,rows,pointer,values,method,verbose,op)
        if (op == 0 .or. op == 2) then
            do j = 1,size(b,2)
                call sparse_solve(solver,b(:,j),status)
                call checked(status)
            end do
        end if
        if (op == 0) call release()
    end subroutine

    subroutine complex_block(nrow,ncol,nz,rows,pointer,values,b,method,verbose,iopt)
        integer, intent(in) :: nrow,ncol,nz,rows(:),pointer(:),method
        complex(dp), intent(in) :: values(:)
        complex(dp), intent(inout) :: b(:,:)
        logical, intent(in) :: verbose
        integer, optional, intent(in) :: iopt
        type(fortsparse_status_t) :: status
        integer :: op,j
        op = operation(method,iopt)
        if (op == 3) then
            call release()
            return
        end if
        if (size(b,1) /= nrow) call fail('incompatible block RHS dimensions')
        call prepare_complex(nrow,ncol,nz,rows,pointer,values,method,verbose,op)
        if (op == 0 .or. op == 2) then
            do j = 1,size(b,2)
                call sparse_solve(solver,b(:,j),status)
                call checked(status)
            end do
        end if
        if (op == 0) call release()
    end subroutine
end module neo2_native_sparse
