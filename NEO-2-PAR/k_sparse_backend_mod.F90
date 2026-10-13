module k_sparse_backend_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
#ifdef NEO2_HAVE_FORTSPARSE
    use fortsparse, only: sparse_solver_t, fortsparse_status_t, status_ok, &
        fs_factor => sparse_factor, fs_solve => sparse_solve, &
        fs_free => sparse_free, sparse_factor_stats, sparse_factor_stats_t, &
        FORTSPARSE_BACKEND_SUPERLU, FORTSPARSE_BACKEND_MUMPS, &
        FORTSPARSE_SYM_GENERAL
#endif
    implicit none
    private
    public :: k_sparse_backend_t, require_k_sparse_backend

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

    subroutine require_k_sparse_backend(backend)
        integer, intent(in) :: backend

        select case (backend)
        case (0)
            return
        case (1)
#ifndef NEO2_HAVE_FORTSPARSE_SUPERLU
            error stop 'k_sparse_backend=1 requires fortsparse SuperLU support.'
#endif
        case (2)
#ifndef NEO2_HAVE_FORTSPARSE_MUMPS
            error stop 'k_sparse_backend=2 requires fortsparse MUMPS support.'
#endif
        case default
            error stop 'k_sparse_backend must be 0 (UMFPACK), 1 (SuperLU), or 2 (MUMPS).'
        end select
    end subroutine require_k_sparse_backend

    subroutine k_factor(self, backend, nrow, ncol, nz, pcol, row, val)
        class(k_sparse_backend_t), intent(inout) :: self
        integer, intent(in) :: backend, nrow, ncol, nz, pcol(:), row(:)
        real(dp), intent(in) :: val(:)
#ifdef NEO2_HAVE_FORTSPARSE
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: stats
#endif

        call require_k_sparse_backend(backend)
        if (backend == 0) return
#ifdef NEO2_HAVE_FORTSPARSE
        select case (backend)
        case (1)
            self%solver%backend_id = FORTSPARSE_BACKEND_SUPERLU
        case (2)
            self%solver%backend_id = FORTSPARSE_BACKEND_MUMPS
        end select
        self%solver%refine = .false.
        self%solver%options%symmetry = FORTSPARSE_SYM_GENERAL
        call fs_factor(self%solver, nrow, ncol, nz, pcol, row, val, status)
        call check_status(self, status)
        self%active = .true.
        write (*,*) 'K sparse backend:', backend
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
