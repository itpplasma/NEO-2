program test_iterator_convergence
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use arnoldi_mod, only: iterator, f_init_arnoldi
    use collisionality_mod, only: num_spec
    use mpiprovider_module, only: mpro
    implicit none
    complex(kind=kind(1d0)) :: source(2), solution(2), diagonal(2), exact(2)
    double precision :: residual
    integer :: status, iterations, rank
    logical :: corrupt_residual, projected_iteration

    call mpro%init()
    num_spec = mpro%getNumProcs()
    rank = mpro%getRank()
    corrupt_residual = .false.
    projected_iteration = .false.
    diagonal = [(0.2d0,0.0d0), (0.3d0,0.0d0)]
    exact = [(2.0d0,1.0d0), (-1.0d0,0.5d0)]
    source = (1.0d0-diagonal)*exact
    solution = source
    call iterator(0, 2, 2, 1.0d-12, 100, solution, rank, apply_iteration, &
                  status, iterations, residual, original_residual)
    if (status /= 0 .or. maxval(abs(solution-exact)) > 1.0d-10) error stop 'direct oracle'
    if (residual > 1.0d-12 .or. iterations <= 1 .or. iterations >= 100) error stop 'direct diagnostics'

    ! A deliberately inadequate budget must fail despite returning an iterate.
    diagonal = (0.99d0, 0.0d0)
    source = (1.0d0-diagonal)*exact
    solution = source
    call iterator(0, 2, 2, 1.0d-12, 1, solution, rank, apply_iteration, &
                  status, iterations, residual, original_residual)
    if (status /= 1 .or. iterations /= 1 .or. residual <= 1.0d-12) error stop 'cap accepted'

    source = (0.0d0,0.0d0)
    solution = source
    call iterator(0, 2, 2, 1.0d-12, 2, solution, rank, apply_iteration, &
                  status, iterations, residual, original_residual)
    if (status /= 0 .or. residual /= 0.0d0) error stop 'zero source'

    if (num_spec == 1) then
        ! Exact eigenvector subtraction: independent solution of (I-A)x=q.
        diagonal = [(2.0d0,0.0d0), (0.25d0,0.0d0)]
        source = (1.0d0-diagonal)*exact
        allocate(f_init_arnoldi(2))
        f_init_arnoldi = source
        solution = source
        call iterator(1, 2, 2, 1.0d-12, 100, solution, rank, apply_iteration, &
                      status, iterations, residual, original_residual)
        if (status /= 0 .or. maxval(abs(solution-exact)) > 1.0d-10) error stop 'stabilized oracle'
        call iterator(3, 2, 2, 1.0d-12, 100, solution, rank, apply_iteration)

        ! Singular but compatible equation: the unit eigenvalue of the
        ! iteration is a null-space mode and must not be subtracted.
        diagonal = [(1.0d0,0.0d0), (0.25d0,0.0d0)]
        source = [(0.0d0,0.0d0), (0.75d0,0.0d0)*exact(2)]
        f_init_arnoldi = [(1.0d0,0.0d0), (1.0d0,0.0d0)]
        solution = source
        call iterator(1, 2, 4, 1.0d-12, 100, solution, rank, apply_iteration, &
                      status, iterations, residual, original_residual)
        if (status /= 0 .or. abs(solution(2)-exact(2)) > 1.0d-10) error stop 'neutral mode oracle'
        call iterator(3, 2, 2, 1.0d-12, 100, solution, rank, apply_iteration)
        deallocate(f_init_arnoldi)
    end if

    ! A projected iteration can settle while violating the original equation.
    ! For L=diag(0,1), q=(1,0), the projected fixed point is f=(1,0),
    ! but the first original equation is 0=1 and its backward error stays one.
    projected_iteration = .true.
    diagonal = [(1.0d0,0.0d0), (0.0d0,0.0d0)]
    source = [(1.0d0,0.0d0), (0.0d0,0.0d0)]
    solution = source
    call iterator(0, 2, 2, 1.0d-12, 60, solution, rank, apply_iteration, &
                  status, iterations, residual, original_residual)
    if (status /= 1 .or. residual < 0.99d0) error stop 'projection hid equation defect'
    if (abs(solution(1)-1.0d0) > 1.0d-12) error stop 'projected fixed point oracle'
    projected_iteration = .false.

    solution = source
    call iterator(1, 2, 1, 1.0d-12, 1, solution, rank, apply_iteration, &
                  status, iterations, residual, original_residual)
    if (status /= 4 .or. iterations /= 0) error stop 'invalid Arnoldi budget'

    ! A single rank's invalid residual must reject the solve on every rank.
    corrupt_residual = .true.
    source = (1.0d0,0.0d0)
    solution = source
    call iterator(0, 2, 2, 1.0d-12, 1, solution, rank, apply_iteration, &
                  status, iterations, residual, original_residual)
    if (status == 0 .or. residual <= 1.0d-12) error stop 'rank failure accepted'

    call mpro%deinit()
    print *, 'All tests passed!'
contains
    subroutine apply_iteration(n, input, output)
        use arnoldi_mod, only: fzero, mode
        integer :: n
        complex(kind=kind(1d0)) :: input(n), output(n)
        output = diagonal*input
        if (mode /= 2) output = output + fzero
        if (projected_iteration) output(1) = 0.5d0*output(1)
    end subroutine apply_iteration

    subroutine original_residual(n, candidate, value)
        integer, intent(in) :: n
        complex(kind=kind(1d0)), intent(in) :: candidate(n)
        double precision, intent(out) :: value
        complex(kind=kind(1d0)) :: lhs(n)
        lhs = (1.0d0-diagonal)*candidate
        value = sum(abs(source-lhs))/max(sum(abs(source))+sum(abs(lhs)), tiny(1.0d0))
        if (corrupt_residual .and. rank == num_spec-1) then
            value = ieee_value(0.0d0, ieee_quiet_nan)
        end if
    end subroutine original_residual
end program test_iterator_convergence
