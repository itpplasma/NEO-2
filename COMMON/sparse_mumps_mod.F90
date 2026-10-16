module sparse_mumps_mod
    ! Optional MUMPS direct solver (via fortsparse) behind sparse_mod's
    ! sparse_solve_method = 4. iopt follows sparse_mod: 0 factor, solve and
    ! free; 1 factor; 2 solve with the stored factor; 3 free.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
#ifdef NEO2_HAVE_FORTSPARSE_MUMPS
    use fortsparse, only: sparse_solver_t, fortsparse_status_t, status_ok, &
                          fs_factor => sparse_factor, fs_solve => sparse_solve, &
                          fs_solve_block => sparse_solve_block, fs_free => sparse_free, &
                          FORTSPARSE_BACKEND_MUMPS, FORTSPARSE_SYM_GENERAL
#endif
    implicit none
    private
    public :: sparse_mumps_available, sparse_solve_mumps

    interface sparse_solve_mumps
        module procedure solve_real_b1, solve_real_b2
        module procedure solve_complex_b1, solve_complex_b2
    end interface sparse_solve_mumps

#ifdef NEO2_HAVE_FORTSPARSE_MUMPS
    ! One retained factor per scalar type, as sparse_mod keeps UMFPACK's.
    type(sparse_solver_t), save :: real_solver, complex_solver
#endif

contains

    logical function sparse_mumps_available()
#ifdef NEO2_HAVE_FORTSPARSE_MUMPS
        sparse_mumps_available = .true.
#else
        sparse_mumps_available = .false.
#endif
    end function sparse_mumps_available

    subroutine solve_real_b1(nrow, ncol, nz, irow, pcol, val, b, iopt)
        integer, intent(in) :: nrow, ncol, nz, irow(:), pcol(:), iopt
        real(dp), intent(in) :: val(:)
        real(dp), intent(inout) :: b(:)
#ifdef NEO2_HAVE_FORTSPARSE_MUMPS
        type(fortsparse_status_t) :: status

        if (iopt == 0 .or. iopt == 1) then
            call prepare(real_solver)
            call fs_factor(real_solver, nrow, ncol, nz, pcol, irow, val, status)
            call check_status(status, 'factor')
        end if
        if (iopt == 0 .or. iopt == 2) then
            call fs_solve(real_solver, b, status)
            call check_status(status, 'solve')
        end if
        if (iopt == 0 .or. iopt == 3) call fs_free(real_solver)
#else
        call unavailable()
#endif
    end subroutine solve_real_b1

    subroutine solve_real_b2(nrow, ncol, nz, irow, pcol, val, b, iopt)
        integer, intent(in) :: nrow, ncol, nz, irow(:), pcol(:), iopt
        real(dp), intent(in) :: val(:)
        real(dp), intent(inout) :: b(:, :)
#ifdef NEO2_HAVE_FORTSPARSE_MUMPS
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: x(:, :)

        if (iopt == 0 .or. iopt == 1) then
            call prepare(real_solver)
            call fs_factor(real_solver, nrow, ncol, nz, pcol, irow, val, status)
            call check_status(status, 'factor')
        end if
        if (iopt == 0 .or. iopt == 2) then
            allocate (x(size(b, 1), size(b, 2)))
            call fs_solve_block(real_solver, b, x, status)
            call check_status(status, 'block solve')
            b = x
        end if
        if (iopt == 0 .or. iopt == 3) call fs_free(real_solver)
#else
        call unavailable()
#endif
    end subroutine solve_real_b2

    subroutine solve_complex_b1(nrow, ncol, nz, irow, pcol, val, b, iopt)
        integer, intent(in) :: nrow, ncol, nz, irow(:), pcol(:), iopt
        complex(dp), intent(in) :: val(:)
        complex(dp), intent(inout) :: b(:)
#ifdef NEO2_HAVE_FORTSPARSE_MUMPS
        type(fortsparse_status_t) :: status

        if (iopt == 0 .or. iopt == 1) then
            call prepare(complex_solver)
            call fs_factor(complex_solver, nrow, ncol, nz, pcol, irow, val, status)
            call check_status(status, 'factor')
        end if
        if (iopt == 0 .or. iopt == 2) then
            call fs_solve(complex_solver, b, status)
            call check_status(status, 'solve')
        end if
        if (iopt == 0 .or. iopt == 3) call fs_free(complex_solver)
#else
        call unavailable()
#endif
    end subroutine solve_complex_b1

    subroutine solve_complex_b2(nrow, ncol, nz, irow, pcol, val, b, iopt)
        integer, intent(in) :: nrow, ncol, nz, irow(:), pcol(:), iopt
        complex(dp), intent(in) :: val(:)
        complex(dp), intent(inout) :: b(:, :)
#ifdef NEO2_HAVE_FORTSPARSE_MUMPS
        type(fortsparse_status_t) :: status
        complex(dp), allocatable :: x(:, :)

        if (iopt == 0 .or. iopt == 1) then
            call prepare(complex_solver)
            call fs_factor(complex_solver, nrow, ncol, nz, pcol, irow, val, status)
            call check_status(status, 'factor')
        end if
        if (iopt == 0 .or. iopt == 2) then
            allocate (x(size(b, 1), size(b, 2)))
            call fs_solve_block(complex_solver, b, x, status)
            call check_status(status, 'block solve')
            b = x
        end if
        if (iopt == 0 .or. iopt == 3) call fs_free(complex_solver)
#else
        call unavailable()
#endif
    end subroutine solve_complex_b2

#ifdef NEO2_HAVE_FORTSPARSE_MUMPS
    subroutine prepare(solver)
        type(sparse_solver_t), intent(inout) :: solver

        ! Each MPI rank factors its own systems (MUMPS runs on MPI_COMM_SELF).
        call fs_free(solver)
        solver%backend_id = FORTSPARSE_BACKEND_MUMPS
        solver%refine = .false.
        solver%options%symmetry = FORTSPARSE_SYM_GENERAL
    end subroutine prepare

    subroutine check_status(status, step)
        type(fortsparse_status_t), intent(in) :: status
        character(*), intent(in) :: step

        if (status_ok(status)) return
        write (error_unit, *) 'MUMPS ', step, ' failed:', status%code, trim(status%msg)
        error stop 'sparse_solve_method = 4 (MUMPS) failed.'
    end subroutine check_status
#else
    subroutine unavailable()
        error stop 'sparse_solve_method = 4 (MUMPS) requires NEO2_USE_FORTSPARSE=ON.'
    end subroutine unavailable
#endif

end module sparse_mumps_mod
