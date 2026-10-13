module k_ilu_backend_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
#ifdef NEO2_HAVE_FORTSPARSE
    use fortsparse, only: sparse_solver_t, fortsparse_status_t, status_ok, &
        sparse_factor, sparse_solve, sparse_free, csc_t, csc_matvec, &
        FORTSPARSE_BACKEND_SUPERLU_ILU
#endif
    implicit none
    private
    public :: k_ilu_backend_t
    type :: k_ilu_backend_t
        logical :: active = .false.
        integer :: solve_count = 0
        real(dp) :: solve_cpu = 0.0_dp
#ifdef NEO2_HAVE_FORTSPARSE
        type(sparse_solver_t) :: solver
        type(csc_t) :: matrix
#endif
    contains
        procedure :: factor => ilu_factor
        procedure :: solve => ilu_solve
        procedure :: matvec => ilu_matvec
        procedure :: free => ilu_free
    end type
contains
    subroutine ilu_factor(self, k_only, iter_solver, method, n, nz, cp, ri, val)
        class(k_ilu_backend_t), intent(inout) :: self
        logical, intent(in) :: k_only
        integer, intent(in) :: iter_solver, method, n, nz, cp(:), ri(:)
        real(dp), intent(in) :: val(:)
        character(16) :: backend
        integer :: env_status, env_length
#ifdef NEO2_HAVE_FORTSPARSE
        type(fortsparse_status_t) :: status
#endif
        if (.not. k_only) return
        call get_environment_variable('NEO2_K_SPARSE_BACKEND', backend, &
            length=env_length, status=env_status)
        if (env_status == 1 .or. env_length == 0) return
        if (env_status /= 0) error stop 'Invalid NEO2_K_SPARSE_BACKEND.'
        if (trim(backend) /= 'superlu-ilu') error stop 'Expected superlu-ilu.'
        if (iter_solver /= 1 .or. method /= 3) &
            error stop 'K ILU requires GMRES and sparse method 3.'
#ifdef NEO2_HAVE_FORTSPARSE
        self%solver%backend_id = FORTSPARSE_BACKEND_SUPERLU_ILU
        self%solver%refine = .false.
        call sparse_factor(self%solver, n, n, nz, cp, ri, val, status)
        call check(status)
        self%matrix%nrow = n
        self%matrix%ncol = n
        self%matrix%nnz = nz
        self%matrix%col_ptr = cp
        self%matrix%row_idx = ri
        self%matrix%val = val
        self%active = .true.
        write (*,*) 'K ILU: drop=1e-4 fill=10, equilibrated COLAMD, GMRES(30)'
#else
        error stop 'K ILU requires NEO2_USE_FORTSPARSE.'
#endif
    end subroutine
    subroutine ilu_solve(self, rhs)
        class(k_ilu_backend_t), intent(inout) :: self
        real(dp), intent(inout) :: rhs(:)
#ifdef NEO2_HAVE_FORTSPARSE
        type(fortsparse_status_t) :: status
        real(dp) :: before, after
        call cpu_time(before)
        call sparse_solve(self%solver, rhs, status)
        call cpu_time(after)
        self%solve_count = self%solve_count + 1
        self%solve_cpu = self%solve_cpu + after - before
        call check(status)
#else
        error stop 'K ILU requires NEO2_USE_FORTSPARSE.'
#endif
    end subroutine
    subroutine ilu_matvec(self, x, y)
        class(k_ilu_backend_t), intent(in) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
#ifdef NEO2_HAVE_FORTSPARSE
        y = csc_matvec(self%matrix, x)
#else
        error stop 'K ILU requires NEO2_USE_FORTSPARSE.'
#endif
    end subroutine
    subroutine ilu_free(self)
        class(k_ilu_backend_t), intent(inout) :: self
#ifdef NEO2_HAVE_FORTSPARSE
        call sparse_free(self%solver)
        self%matrix = csc_t()
#endif
        self%active = .false.
    end subroutine
#ifdef NEO2_HAVE_FORTSPARSE
    subroutine check(status)
        type(fortsparse_status_t), intent(in) :: status
        if (status_ok(status)) then
            if (len_trim(status%msg) > 0) write (*,*) trim(status%msg)
            return
        end if
        write (error_unit,*) 'K ILU error:', status%code, trim(status%msg)
        error stop 'K ILU failed.'
    end subroutine
#endif
end module
