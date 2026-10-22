module k_sparse_backend_mod
    ! Optional fortsparse MUMPS factor reused for all K-only GMRES solves.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
#ifdef NEO2_HAVE_FORTSPARSE_MUMPS
    use fortsparse, only: sparse_solver_t, fortsparse_status_t, status_ok, &
        fs_factor => sparse_factor, fs_solve => sparse_solve, &
        fs_free => sparse_free, FORTSPARSE_BACKEND_MUMPS, FORTSPARSE_SYM_GENERAL
#endif
    implicit none
    private
    public :: k_sparse_backend_t, require_k_sparse_backend

    type :: k_sparse_backend_t
#ifdef NEO2_HAVE_FORTSPARSE_MUMPS
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
        case (2)
#ifndef NEO2_HAVE_FORTSPARSE_MUMPS
            error stop 'k_sparse_backend=2 requires NEO2_USE_FORTSPARSE with MUMPS.'
#endif
        case default
            error stop 'k_sparse_backend must be 0 (UMFPACK) or 2 (MUMPS).'
        end select
    end subroutine require_k_sparse_backend

    subroutine k_factor(self, nrow, ncol, nz, pcol, row, val)
        class(k_sparse_backend_t), intent(inout) :: self
        integer, intent(in) :: nrow, ncol, nz, pcol(:), row(:)
        real(dp), intent(in) :: val(:)
#ifdef NEO2_HAVE_FORTSPARSE_MUMPS
        type(fortsparse_status_t) :: status

        self%solver%backend_id = FORTSPARSE_BACKEND_MUMPS
        self%solver%refine = .false.
        self%solver%options%symmetry = FORTSPARSE_SYM_GENERAL
        call fs_factor(self%solver, nrow, ncol, nz, pcol, row, val, status)
        call check_status(status)
#else
        error stop 'k_sparse_backend requires NEO2_USE_FORTSPARSE with MUMPS.'
#endif
    end subroutine k_factor

    subroutine k_solve(self, rhs)
        class(k_sparse_backend_t), intent(inout) :: self
        real(dp), intent(inout) :: rhs(:)
#ifdef NEO2_HAVE_FORTSPARSE_MUMPS
        type(fortsparse_status_t) :: status

        call fs_solve(self%solver, rhs, status)
        call check_status(status)
#else
        error stop 'k_sparse_backend requires NEO2_USE_FORTSPARSE with MUMPS.'
#endif
    end subroutine k_solve

    subroutine k_free(self)
        class(k_sparse_backend_t), intent(inout) :: self

#ifdef NEO2_HAVE_FORTSPARSE_MUMPS
        call fs_free(self%solver)
#endif
    end subroutine k_free

#ifdef NEO2_HAVE_FORTSPARSE_MUMPS
    subroutine check_status(status)
        type(fortsparse_status_t), intent(in) :: status

        if (status_ok(status)) return
        write (error_unit, *) 'K sparse backend error:', status%code, trim(status%msg)
        error stop 'K sparse backend failed.'
    end subroutine check_status
#endif
end module k_sparse_backend_mod
