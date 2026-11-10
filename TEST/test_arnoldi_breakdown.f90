program test_arnoldi_breakdown
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use arnoldi_mod, only: arnoldi, eigvecs, f_init_arnoldi, ierr, iterator, &
        & ngrow, ntol, ritznum, tol
    use collisionality_mod, only: num_spec
    use mpiprovider_module, only: mpro
    implicit none

    integer :: applications, status
    complex(kind=kind(1d0)) :: solution(2)
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
    f_init_arnoldi = [(1.0d0, 0.0d0), (1.0d0, 0.0d0)]
    solution = [(0.0d0, 0.0d0), (0.75d0, 0.0d0)]
    call iterator(1, 2, 4, 1.0d-12, 100, solution, 0, neutral_operator)
    if (abs(solution(2) - 1.0d0) > 1.0d-10) then
        print *, 'FAIL: singular compatible solve:', solution(2), ' expected 1'
        status = status + 1
    end if
    call iterator(3, 2, 4, 1.0d-12, 100, solution, 0, neutral_operator)

    call nonclosed_projection()

    call mpro%deinit()
    deallocate(f_init_arnoldi)
    if (allocated(ritznum)) deallocate(ritznum)
    if (allocated(eigvecs)) deallocate(eigvecs)

    if (status /= 0) error stop
    print *, 'All tests passed!'

contains

    subroutine nonclosed_projection()
        integer :: j
        double precision :: expected(2)
        complex(kind=kind(1d0)) :: residual(3), q(3,2)

        deallocate(f_init_arnoldi)
        if (allocated(ritznum)) deallocate(ritznum)
        allocate(f_init_arnoldi(3), ritznum(2))
        f_init_arnoldi = (1.0d0, 0.0d0)
        applications = 0
        call arnoldi(3, 2, 0, diagonal_operator)
        ! For diag(2,3,4) and the constant start, Q* A Q has diagonal
        ! (3,3) and off-diagonal sqrt(2/3); the residual is nonzero.
        expected = [3.0d0 + sqrt(2.0d0/3.0d0), 3.0d0 - sqrt(2.0d0/3.0d0)]
        q(:,1) = 1.0d0/sqrt(3.0d0)
        q(:,2) = [-1.0d0, 0.0d0, 1.0d0]/sqrt(2.0d0)
        if (ierr /= 0 .or. ngrow /= 2 .or. applications /= 2) then
            print *, 'FAIL: nonclosed projection:', ierr, ngrow, applications
            status = status + 1
            return
        end if
        do j = 1,2
            if (.not.ieee_is_finite(abs(ritznum(j))) .or. &
                & abs(ritznum(j)-expected(j)) > 1.0d-12) then
                print *, 'FAIL: nonclosed Ritz value:', ritznum(j), expected(j)
                status = status + 1
            end if
            residual = [2.0d0,3.0d0,4.0d0]*eigvecs(:,j)-ritznum(j)*eigvecs(:,j)
            if (.not.all(ieee_is_finite(abs(residual))) .or. &
                & maxval(abs(matmul(conjg(transpose(q)),residual))) > 1.0d-12 .or. &
                & maxval(abs(residual)) < 1.0d-3) then
                print *, 'FAIL: nonclosed Galerkin residual:', residual
                status = status + 1
            end if
        end do
    end subroutine nonclosed_projection

    subroutine diagonal_operator(n, input, output)
        integer :: n
        complex(kind=kind(1d0)), dimension(n) :: input, output

        applications = applications + 1
        output = [2.0d0,3.0d0,4.0d0]*input
    end subroutine diagonal_operator

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

        output = [input(1), 0.25d0 * input(2)]
        if (mode /= 2) output = output + fzero
    end subroutine neutral_operator

end program test_arnoldi_breakdown
