module k_sparse_backend_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
#ifdef NEO2_HAVE_FORTSPARSE
    use fortsparse, only: sparse_solver_t, fortsparse_status_t, status_ok, &
        fs_factor => sparse_factor, fs_solve => sparse_solve, &
        fs_free => sparse_free, sparse_factor_stats, sparse_factor_stats_t, &
        FORTSPARSE_BACKEND_SUPERLU, FORTSPARSE_BACKEND_UMFPACK_IPC, &
        FORTSPARSE_BACKEND_MUMPS, FORTSPARSE_SYM_GENERAL
#endif
    implicit none
    private
    public :: k_sparse_backend_t

    type :: k_sparse_backend_t
        logical :: active = .false.
#ifdef NEO2_HAVE_FORTSPARSE
        type(sparse_solver_t) :: solver
#endif
    contains
        procedure :: factor => k_factor
        procedure :: solve => k_solve
        procedure :: free => k_free
    end type k_sparse_backend_t

contains

    subroutine k_factor(self, k_only, method, nrow, ncol, nz, pcol, row, val)
        class(k_sparse_backend_t), intent(inout) :: self
        logical, intent(in) :: k_only
        integer, intent(in) :: method, nrow, ncol, nz, pcol(:), row(:)
        real(dp), intent(in) :: val(:)
        character(16) :: backend
        integer :: env_status, env_length
#ifdef NEO2_HAVE_FORTSPARSE
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: stats
#endif

        if (.not. k_only) return
        call get_environment_variable('NEO2_K_SPARSE_BACKEND', backend, &
            length=env_length, status=env_status)
        if (env_status == 1 .or. env_length == 0) return
        if (env_status /= 0) error stop 'Invalid NEO2_K_SPARSE_BACKEND.'
        if (method /= 3) error stop 'K sparse backend requires method 3.'
#ifdef NEO2_HAVE_FORTSPARSE
        select case (trim(backend))
        case ('superlu')
            self%solver%backend_id = FORTSPARSE_BACKEND_SUPERLU
        case ('mumps')
            self%solver%backend_id = FORTSPARSE_BACKEND_MUMPS
        case ('umfpack-ipc')
            self%solver%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC
        case default
            error stop 'Unknown NEO2_K_SPARSE_BACKEND.'
        end select
        self%solver%refine = .false.
        self%solver%options%symmetry = FORTSPARSE_SYM_GENERAL
        call fs_factor(self%solver, nrow, ncol, nz, pcol, row, val, status)
        call check_status(self, status)
        self%active = .true.
        write (*,*) 'K sparse backend: ', trim(backend)
        if (self%solver%backend_id == FORTSPARSE_BACKEND_MUMPS) then
            call sparse_factor_stats(self%solver, stats, status)
            call check_status(self, status)
            write (*,*) 'K sparse factor statistics:', stats
        end if
#else
        error stop 'K sparse backend requires NEO2_USE_FORTSPARSE.'
#endif
    end subroutine k_factor

    subroutine k_solve(self, rhs)
        class(k_sparse_backend_t), intent(inout) :: self
        real(dp), intent(inout) :: rhs(:)
#ifdef NEO2_HAVE_FORTSPARSE
        type(fortsparse_status_t) :: status

        call fs_solve(self%solver, rhs, status)
        call check_status(self, status)
#else
        error stop 'K sparse backend requires NEO2_USE_FORTSPARSE.'
#endif
    end subroutine k_solve

    subroutine k_free(self)
        class(k_sparse_backend_t), intent(inout) :: self

#ifdef NEO2_HAVE_FORTSPARSE
        call fs_free(self%solver)
#endif
        self%active = .false.
    end subroutine k_free

#ifdef NEO2_HAVE_FORTSPARSE
    subroutine check_status(self, status)
        class(k_sparse_backend_t), intent(inout) :: self
        type(fortsparse_status_t), intent(in) :: status

        if (status_ok(status)) return
        write (error_unit, *) 'K sparse backend error:', status%code, &
            trim(status%msg)
        call self%free()
        error stop 'K sparse backend failed.'
    end subroutine check_status
#endif
end module k_sparse_backend_mod
