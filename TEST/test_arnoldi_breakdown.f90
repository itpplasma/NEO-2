program test_arnoldi_breakdown
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use arnoldi_mod, only: arnoldi, eigvecs, f_init_arnoldi, ierr, ngrow, &
        & ntol, ritznum, tol
    use collisionality_mod, only: num_spec
    use mpiprovider_module, only: mpro
    implicit none

    integer :: applications, status
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

    call mpro%deinit()
    deallocate(f_init_arnoldi, ritznum)
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

end program test_arnoldi_breakdown
