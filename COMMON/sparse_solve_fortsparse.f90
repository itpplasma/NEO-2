module sparse_solve_fortsparse
    ! Thin adapter that routes NEO-2's sparse solves through fortsparse's
    ! UMFPACK-IPC backend. Same UMFPACK numerics as the old in-process umf4
    ! wrapper, but reached out-of-process so neo_2_*.x no longer links
    ! SuiteSparse/UMFPACK. Two module-global solver handles (real, complex)
    ! retain a factorization across the factor/solve/free lifecycle that the
    ! callers in sparse_mod drive through iopt. Factorization streams the CSC
    ! arrays straight into the helper's shared memory (raw sparse_factor), so
    ! this adapter adds no matrix-sized copy of its own.
    use fortsparse, only: sparse_solver_t, sparse_factor, &
        sparse_solve, sparse_free, sparse_vector, fortsparse_status_t, &
        status_ok, FORTSPARSE_BACKEND_UMFPACK_IPC
    implicit none
    private

    integer, parameter :: dp = kind(1.0d0)

    public :: fs_solve_real
    public :: fs_solve_complex
    public :: fs_vector_real
    public :: fs_solve_real_2

    type(sparse_solver_t), save :: real_solver = &
        sparse_solver_t(backend_id=FORTSPARSE_BACKEND_UMFPACK_IPC)
    type(sparse_solver_t), save :: complex_solver = &
        sparse_solver_t(backend_id=FORTSPARSE_BACKEND_UMFPACK_IPC)

contains

    ! Solve A x = b for a real 1-D RHS. iopt drives the lifecycle: 1 factor,
    ! 2 solve, 3 free, 0 all three. Solution is returned in b.
    subroutine fs_solve_real(nrow, ncol, nz, irow, pcol, val, b, iopt, refine)
        integer,       intent(in)    :: nrow, ncol, nz
        integer,       intent(in)    :: irow(:), pcol(:)
        real(dp),      intent(in)    :: val(:)
        real(dp),      intent(inout) :: b(:)
        integer,       intent(in)    :: iopt
        logical,       intent(in)    :: refine

        type(fortsparse_status_t) :: status

        real_solver%refine = refine
        if (iopt == 1 .or. iopt == 0) then
            call sparse_factor(real_solver, nrow, ncol, nz, pcol, irow, val, &
                status)
            call check(status, 'fs_solve_real: factor')
        end if
        if (iopt == 2 .or. iopt == 0) then
            call sparse_solve(real_solver, b, status)
            call check(status, 'fs_solve_real: solve')
        end if
        if (iopt == 3 .or. iopt == 0) call sparse_free(real_solver)
    end subroutine fs_solve_real

    ! Solve A x = b for a complex 1-D RHS. Lifecycle as in fs_solve_real.
    subroutine fs_solve_complex(nrow, ncol, nz, irow, pcol, val, b, iopt, refine)
        integer,       intent(in)    :: nrow, ncol, nz
        integer,       intent(in)    :: irow(:), pcol(:)
        complex(dp),   intent(in)    :: val(:)
        complex(dp),   intent(inout) :: b(:)
        integer,       intent(in)    :: iopt
        logical,       intent(in)    :: refine

        type(fortsparse_status_t) :: status

        complex_solver%refine = refine
        if (iopt == 1 .or. iopt == 0) then
            call sparse_factor(complex_solver, nrow, ncol, nz, pcol, irow, &
                val, status)
            call check(status, 'fs_solve_complex: factor')
        end if
        if (iopt == 2 .or. iopt == 0) then
            call sparse_solve(complex_solver, b, status)
            call check(status, 'fs_solve_complex: solve')
        end if
        if (iopt == 3 .or. iopt == 0) call sparse_free(complex_solver)
    end subroutine fs_solve_complex

    ! A length-n real solve vector owned by the (already factored) real solver,
    ! living in the helper's shared mapping. Filling it as the RHS and reading
    ! the solution from such a vector makes the solve copy nothing across the
    ! process boundary. Valid until the next factor or the iopt=3 free (which
    ! releases the mapping); the caller never frees it.
    function fs_vector_real(n) result(p)
        integer, intent(in) :: n
        real(dp), pointer   :: p(:)

        p => sparse_vector(real_solver, n)
    end function fs_vector_real

    ! Two-vector real solve: A sol = rhs, reusing the resident factorization.
    ! When rhs and sol are fs_vector_real vectors the helper reads and writes
    ! them in shared memory directly, with no marshalling. The matrix must
    ! already be factored (fs_solve_real with iopt=1).
    subroutine fs_solve_real_2(rhs, sol, refine)
        real(dp), target, contiguous, intent(in)  :: rhs(:)
        real(dp), target, contiguous, intent(out) :: sol(:)
        logical,                      intent(in)  :: refine

        type(fortsparse_status_t) :: status

        real_solver%refine = refine
        call sparse_solve(real_solver, rhs, sol, status)
        call check(status, 'fs_solve_real_2: solve')
    end subroutine fs_solve_real_2

    ! Abort with a diagnostic on any non-ok fortsparse status, mirroring the
    ! old umf4 PRINT-on-error-then-stop behavior.
    subroutine check(status, where)
        type(fortsparse_status_t), intent(in) :: status
        character(*),              intent(in) :: where

        if (.not. status_ok(status)) then
            print *, 'fortsparse error in ', where, ': ', trim(status%msg)
            stop 1
        end if
    end subroutine check

end module sparse_solve_fortsparse
