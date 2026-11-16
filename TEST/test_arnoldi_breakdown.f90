program test_arnoldi_breakdown
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use arnoldi_mod, only: arnoldi, eigvecs, f_init_arnoldi, ierr, iterator, &
        & ngrow, ntol, ritznum, tol
    use collisionality_mod, only: num_spec
    use mpiprovider_module, only: mpro
    implicit none

    integer :: applications, status, trial
    logical :: converged
    complex(kind=kind(1d0)) :: solution(2)
    double precision :: neutral_eigenvalue
    double precision, parameter :: eigenvalue = 2.0d0
    double precision, parameter :: tolerance = 100.0d0 * epsilon(1.0d0)

    applications = 0
    status = 0
    num_spec = 1
    tol = 0.5d0
    ntol = 2
    allocate(f_init_arnoldi(2), ritznum(2))
    f_init_arnoldi = [(1.0d0, 0.0d0), (0.0d0, 0.0d0)]

    call mpro%init()
    call arnoldi(2, 2, 0, rank_one_operator)

    if (ierr /= 0) then
        print *, 'FAIL: exact Krylov breakdown returned an error'
        status = status + 1
    end if
    if (applications /= 1) then
        print *, 'FAIL: operator applications:', applications, ' expected 1'
        status = status + 1
    end if
    if (ngrow /= 1) then
        print *, 'FAIL: retained Ritz values:', ngrow, ' expected 1'
        status = status + 1
    else
        if (.not. ieee_is_finite(real(ritznum(1), kind=kind(1d0))) .or. &
            & .not. ieee_is_finite(aimag(ritznum(1)))) then
            print *, 'FAIL: Ritz value has non-finite components'
            status = status + 1
        end if
        if (.not. all(ieee_is_finite(real(eigvecs(:, 1), kind=kind(1d0)))) .or. &
            & .not. all(ieee_is_finite(aimag(eigvecs(:, 1))))) then
            print *, 'FAIL: Ritz vector has non-finite components'
            status = status + 1
        end if
        if (abs(ritznum(1) - eigenvalue) > tolerance) then
            print *, 'FAIL: Ritz value:', ritznum(1), ' expected:', eigenvalue
            status = status + 1
        end if
        if (abs(abs(eigvecs(1, 1)) - 1.0d0) > tolerance) then
            print *, 'FAIL: Ritz vector does not span the known eigenspace'
            status = status + 1
        end if
        if (abs(eigvecs(2, 1)) > tolerance) then
            print *, 'FAIL: Ritz vector has a component outside the known eigenspace'
            status = status + 1
        end if
    end if

    applications = 0
    f_init_arnoldi = (0.0d0, 0.0d0)
    call arnoldi(2, 2, 0, rank_one_operator)
    if (ierr /= 0) then
        print *, 'FAIL: zero initial vector returned an error'
        status = status + 1
    end if
    if (applications /= 0) then
        print *, 'FAIL: zero initial vector applied the operator'
        status = status + 1
    end if
    if (ngrow /= 0) then
        print *, 'FAIL: zero initial vector retained Ritz values'
        status = status + 1
    end if

    applications = 0
    f_init_arnoldi = (0.0d0, 0.0d0)
    f_init_arnoldi(1) = cmplx(ieee_value(0.0d0, ieee_quiet_nan), 0.0d0, &
        & kind=kind(1d0))
    call arnoldi(2, 2, 0, rank_one_operator)
    if (ierr == 0) then
        print *, 'FAIL: non-finite initial vector was accepted'
        status = status + 1
    end if
    if (applications /= 0) then
        print *, 'FAIL: non-finite initial vector applied the operator'
        status = status + 1
    end if

    ! Singular but compatible L f = q with L = I - A, A = diag(1, 1/4):
    ! the unit eigenvalue of A is a null-space mode and must not be
    ! subtracted. The second component of the solution is exactly one.
    neutral_eigenvalue = 1.0d0
    f_init_arnoldi = [(1.0d0, 0.0d0), (1.0d0, 0.0d0)]
    solution = [(0.0d0, 0.0d0), (0.75d0, 0.0d0)]
    call iterator(1, 2, 4, 1.0d-12, 100, solution, 0, neutral_operator)
    if (abs(solution(2) - 1.0d0) > 1.0d-10) then
        print *, 'FAIL: singular compatible solve:', solution(2), ' expected 1'
        status = status + 1
    end if
    call iterator(3, 2, 4, 1.0d-12, 100, solution, 0, neutral_operator)

    ! A unit mode split by less than the Ritz accuracy must stay out of
    ! the subtraction solve, on either side of one. The compatible source
    ! excites only the well-separated eigenvalue 1/4, with solution two = 1.
    do trial = 1, 2
        neutral_eigenvalue = 1.0d0 + (-1.0d0)**trial * 1.05d-6
        f_init_arnoldi = [(1.0d0, 0.0d0), (1.0d0, 0.0d0)]
        solution = [(0.0d0, 0.0d0), (0.75d0, 0.0d0)]
        call iterator(1, 2, 4, 1.0d-12, 100, solution, 0, neutral_operator)
        if (ngrow /= 0 .or. abs(solution(2) - 1.0d0) > 1.0d-10) then
            print *, 'FAIL: unresolved unit mode:', neutral_eigenvalue, ngrow
            status = status + 1
        end if
        call iterator(3, 2, 4, 1.0d-12, 100, solution, 0, neutral_operator)
    end do

    ! A resolved unstable eigenvalue still needs subtraction: (I-A)^-1 q
    ! is exactly (-10, 1) for A = diag(1.1, 1/4), q = (1, 3/4).
    neutral_eigenvalue = 1.1d0
    solution = [(1.0d0, 0.0d0), (0.75d0, 0.0d0)]
    call iterator(1, 2, 4, 1.0d-12, 100, solution, 0, neutral_operator)
    if (ngrow /= 1 .or. &
        & maxval(abs(solution - [(-10.0d0, 0.0d0), (1.0d0, 0.0d0)])) > 1.0d-10) then
        print *, 'FAIL: resolved unstable mode:', solution, ngrow
        status = status + 1
    end if
    call iterator(3, 2, 4, 1.0d-12, 100, solution, 0, neutral_operator)

    ! Direct iterations contract by 0.9 per step: too few must be reported
    ! as not converged instead of returning the unfinished iterate.
    neutral_eigenvalue = 0.9d0
    solution = [(1.0d0, 0.0d0), (0.75d0, 0.0d0)]
    call iterator(0, 2, 4, 1.0d-12, 5, solution, 0, neutral_operator, converged)
    if (converged) then
        print *, 'FAIL: exhausted iteration reported as converged'
        status = status + 1
    end if
    solution = [(1.0d0, 0.0d0), (0.75d0, 0.0d0)]
    call iterator(0, 2, 4, 1.0d-12, 400, solution, 0, neutral_operator, converged)
    if (.not. converged .or. &
        & maxval(abs(solution - [(10.0d0, 0.0d0), (1.0d0, 0.0d0)])) > 1.0d-9) then
        print *, 'FAIL: converged direct iteration:', converged, solution
        status = status + 1
    end if

    call mpro%deinit()
    deallocate(f_init_arnoldi)
    if (allocated(ritznum)) deallocate(ritznum)
    if (allocated(eigvecs)) deallocate(eigvecs)

    if (status /= 0) error stop
    print *, 'All tests passed!'

contains

    subroutine rank_one_operator(n, input, output)
        integer :: n
        complex(kind=kind(1d0)), dimension(n) :: input, output

        applications = applications + 1
        output = (0.0d0, 0.0d0)
        output(1) = eigenvalue * input(1)
    end subroutine rank_one_operator

    subroutine neutral_operator(n, input, output)
        use arnoldi_mod, only: fzero, mode
        integer :: n
        complex(kind=kind(1d0)), dimension(n) :: input, output

        output = [neutral_eigenvalue * input(1), 0.25d0 * input(2)]
        if (mode /= 2) output = output + fzero
    end subroutine neutral_operator

end program test_arnoldi_breakdown
