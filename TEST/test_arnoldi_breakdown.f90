program test_arnoldi_breakdown
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use arnoldi_mod, only: arnoldi, eigvecs, f_init_arnoldi, ierr, ngrow, &
        & ntol, ritznum, tol
    use collisionality_mod, only: num_spec
    use mpiprovider_module, only: mpro
    implicit none

    integer :: applications, status
    complex(kind=kind(1d0)) :: test_matrix(3,3)
    character(32) :: argument
    double precision, parameter :: eigenvalue = 2.0d0
    double precision, parameter :: tolerance = 100.0d0 * epsilon(1.0d0)

    applications = 0
    status = 0
    tol = 0.5d0
    ntol = 2
    call get_command_argument(1,argument)
    if (trim(argument) == '--distributed') then
        call distributed_projection()
        stop
    end if
    num_spec = 1
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

    call complete_diagonal_projection()
    call nonnormal_projection()

    call mpro%deinit()
    deallocate(f_init_arnoldi, ritznum)
    if (allocated(eigvecs)) deallocate(eigvecs)

    if (status /= 0) error stop
    print *, 'All tests passed!'

contains

    subroutine distributed_projection()
        integer :: rank,j
        complex(kind=kind(1d0)) :: expected(2),residual

        call mpro%init()
        rank = mpro%getRank()
        num_spec = mpro%getNumProcs()
        if (num_spec /= 2) then
            call mpro%deinit()
            error stop 'distributed projection requires two MPI processes'
        end if
        allocate(f_init_arnoldi(1),ritznum(2))
        f_init_arnoldi = (1.0d0,0.0d0)
        ! One row per species: the global operator is diag(2,3), although
        ! each local vector has length one.  Krylov capacity is global.
        call arnoldi(1,2,rank,distributed_operator)
        expected = [cmplx(3.0d0,0.0d0,kind(1d0)), &
            & cmplx(2.0d0,0.0d0,kind(1d0))]
        if (ierr /= 0 .or. ngrow /= 2 .or. applications /= 2) then
            print *, 'FAIL: distributed complete projection:', ierr,ngrow,applications
            status = status + 1
        else
            do j = 1,2
                residual = ((2.0d0+rank)-ritznum(j))*eigvecs(1,j)
                if (.not.ieee_is_finite(real(residual,kind(1d0))) .or. &
                    & .not.ieee_is_finite(aimag(residual))) then
                    print *, 'FAIL: non-finite distributed true residual:',rank,j
                    status = status + 1
                else if (abs(ritznum(j)-expected(j)) > 1.0d-12 .or. &
                    & abs(residual) > 1.0d-12) then
                    print *, 'FAIL: distributed spectrum or true residual:',rank,j
                    status = status + 1
                end if
            end do
        end if
        call mpro%deinit()
        deallocate(f_init_arnoldi,ritznum)
        if (allocated(eigvecs)) deallocate(eigvecs)
        if (status /= 0) error stop
        print *, 'All tests passed!'
    end subroutine distributed_projection

    subroutine distributed_operator(n,input,output)
        integer :: n
        complex(kind=kind(1d0)), dimension(n) :: input,output

        applications = applications + 1
        output = (2.0d0+mpro%getRank())*input
    end subroutine distributed_operator

    subroutine complete_diagonal_projection()
        integer :: j
        complex(kind=kind(1d0)) :: expected(2)

        applications = 0
        f_init_arnoldi = (1.0d0, 0.0d0)
        test_matrix = (0.0d0, 0.0d0)
        test_matrix(1,1) = 2.0d0
        test_matrix(2,2) = 3.0d0
        expected = [cmplx(3.0d0,0.0d0,kind(1d0)), cmplx(2.0d0,0.0d0,kind(1d0))]
        call arnoldi(2, 2, 0, matrix_operator)
        if (ierr /= 0 .or. ngrow /= 2 .or. applications /= 2) then
            print *, 'FAIL: complete two-dimensional projection:', &
                & ierr, ngrow, applications
            print *, 'FAIL: diagonal Ritz spectrum:', ritznum
            status = status + 1
            return
        end if
        do j = 1,2
            if (abs(ritznum(j)-expected(j)) > 1.0d-12) then
                print *, 'FAIL: diagonal Ritz value:', ritznum(j), &
                    & ' expected:', expected(j)
                status = status + 1
            end if
            call check_true_residual(2,j)
        end do

        ! A capacity of one still requires the first, complete projection.
        applications = 0
        call arnoldi(2, 1, 0, matrix_operator)
        if (ierr /= 0 .or. ngrow /= 1 .or. applications /= 1) then
            print *, 'FAIL: one-dimensional projection:', ierr, ngrow, applications
            status = status + 1
        else if (abs(ritznum(1)-2.5d0) > 1.0d-12) then
            print *, 'FAIL: first Rayleigh quotient:', ritznum(1), ' expected 2.5'
            status = status + 1
        end if
    end subroutine complete_diagonal_projection

    subroutine nonnormal_projection()
        integer :: j
        complex(kind=kind(1d0)) :: expected(3)

        deallocate(f_init_arnoldi,ritznum)
        allocate(f_init_arnoldi(3),ritznum(3))
        f_init_arnoldi = (1.0d0,0.0d0)
        test_matrix = (0.0d0,0.0d0)
        expected = [cmplx(2.0d0,1.0d0,kind(1d0)), &
            & cmplx(3.0d0,-2.0d0,kind(1d0)), cmplx(-1.0d0,0.5d0,kind(1d0))]
        do j = 1,3
            test_matrix(j,j) = expected(j)
        end do
        test_matrix(1,2) = cmplx(16.0d0,8.0d0,kind(1d0))
        test_matrix(1,3) = cmplx(-4.0d0,2.0d0,kind(1d0))
        test_matrix(2,3) = cmplx(1.0d0,0.5d0,kind(1d0))
        applications = 0
        call arnoldi(3, 3, 0, matrix_operator)
        if (ierr /= 0 .or. ngrow /= 3 .or. applications /= 3) then
            print *, 'FAIL: complete nonnormal projection:', ierr, ngrow, applications
            status = status + 1
            return
        end if
        do j = 1,3
            if (minval(abs(ritznum-expected(j))) > 1.0d-11) then
                print *, 'FAIL: missing nonnormal eigenvalue:', expected(j)
                status = status + 1
            end if
            call check_true_residual(3,j)
        end do
    end subroutine nonnormal_projection

    subroutine check_true_residual(n,j)
        integer, intent(in) :: n,j
        complex(kind=kind(1d0)) :: residual(n)
        double precision :: scale

        residual = matmul(test_matrix(1:n,1:n),eigvecs(:,j))-ritznum(j)*eigvecs(:,j)
        scale = maxval(sum(abs(test_matrix(1:n,1:n)),dim=2))*maxval(abs(eigvecs(:,j)))
        if (.not.all(ieee_is_finite(real(residual,kind(1d0)))) .or. &
            & .not.all(ieee_is_finite(aimag(residual)))) then
            print *, 'FAIL: non-finite true eigenvector residual'
            status = status + 1
        else if (maxval(abs(residual)) > 1.0d-12*scale) then
            print *, 'FAIL: true eigenvector residual:', maxval(abs(residual)), scale
            status = status + 1
        end if
    end subroutine check_true_residual

    subroutine matrix_operator(n,input,output)
        integer :: n
        complex(kind=kind(1d0)), dimension(n) :: input,output

        applications = applications + 1
        output = matmul(test_matrix(1:n,1:n),input)
    end subroutine matrix_operator

    subroutine rank_one_operator(n, input, output)
        integer :: n
        complex(kind=kind(1d0)), dimension(n) :: input, output

        applications = applications + 1
        output = (0.0d0, 0.0d0)
        output(1) = eigenvalue * input(1)
    end subroutine rank_one_operator

end program test_arnoldi_breakdown
