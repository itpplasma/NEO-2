program test_conservation_step
    ! Two-component model of the axisymmetric step: preconditioner P = I, filter
    ! vector e1, and operator A = diag(0.1, 1) for which e1 is not a null vector.
    ! Oracle: f = A^-1 q = (10, 2). The historical iterate filter settles at
    ! f1 = q1/1.1 with backward error 0.18; the update filter reaches f.
    ! Three Arnoldi steps close the two-dimensional Krylov space.
    use arnoldi_mod, only: iterator, f_init_arnoldi, conservation_step
    use collisionality_mod, only: num_spec
    use mpiprovider_module, only: mpro
    implicit none
    complex(kind=kind(1d0)), parameter :: diagonal(2) = [(0.1d0,0.0d0), (1.0d0,0.0d0)]
    complex(kind=kind(1d0)), parameter :: e1(2) = [(1.0d0,0.0d0), (0.0d0,0.0d0)]
    complex(kind=kind(1d0)) :: source(2), solution(2), exact(2)
    double precision :: residual
    integer :: status, iterations
    logical :: filter_update

    call mpro%init()
    num_spec = 1
    exact = [(10.0d0,0.0d0), (2.0d0,0.0d0)]
    source = diagonal*exact
    allocate(f_init_arnoldi(2))

    filter_update = .true.
    f_init_arnoldi = source
    solution = source
    call iterator(1, 2, 3, 1.0d-12, 100, solution, 0, apply_step, &
                  status, iterations, residual, original_residual)
    if (status /= 0 .or. maxval(abs(solution-exact)) > 1.0d-9) then
        print *, 'FAIL: update filter', status, iterations, residual, solution
        error stop
    end if
    call iterator(3, 2, 3, 1.0d-12, 100, solution, 0, apply_step)

    filter_update = .false.
    f_init_arnoldi = source
    solution = source
    call iterator(1, 2, 3, 1.0d-12, 100, solution, 0, apply_step, &
                  status, iterations, residual, original_residual)
    if (status /= 1 .or. residual < 0.1d0) then
        print *, 'FAIL: iterate filter bias not detected', status, residual
        error stop
    end if
    call iterator(3, 2, 3, 1.0d-12, 100, solution, 0, apply_step)

    call mpro%deinit()
    print *, 'All tests passed!'
contains
    subroutine apply_step(n, fold, fnew)
        use arnoldi_mod, only: fzero, mode
        integer :: n
        complex(kind=kind(1d0)) :: fold(n), fnew(n)
        fnew = -diagonal*fold
        if (mode /= 2) fnew = fnew + fzero
        call conservation_step(fold, fnew, filter_update, e1, e1, sum(e1*e1))
    end subroutine apply_step

    subroutine original_residual(n, candidate, value)
        integer, intent(in) :: n
        complex(kind=kind(1d0)), intent(in) :: candidate(n)
        double precision, intent(out) :: value
        complex(kind=kind(1d0)) :: lhs(n)
        lhs = diagonal*candidate
        value = sum(abs(source-lhs))/max(sum(abs(source))+sum(abs(lhs)), tiny(1.0d0))
    end subroutine original_residual
end program test_conservation_step
