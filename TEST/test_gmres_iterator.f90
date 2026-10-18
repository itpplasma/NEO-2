program test_gmres_iterator
    use arnoldi_mod, only: gmres_iterator
    use collisionality_mod, only: num_spec
    use mpiprovider_module, only: mpro
    implicit none
    integer, parameter :: n = 4
    complex(kind=kind(1d0)) :: amat(n, n), source(n), solution(n), exact(n)
    double precision :: residual, c, s
    integer :: status, iterations, rank
    logical :: projected_iteration

    call mpro%init()
    num_spec = mpro%getNumProcs()
    rank = mpro%getRank()
    projected_iteration = .false.

    ! Richardson-hostile map: a rotation pair with |lambda| = 0.999 and a real
    ! eigenvalue 1.5 outside the unit circle. Oracle: q = (I-A) x_exact.
    c = 0.999d0*cos(1.0d0)
    s = 0.999d0*sin(1.0d0)
    amat = (0.0d0, 0.0d0)
    amat(1, 1:2) = [c, -s]
    amat(2, 1:2) = [s, c]
    amat(3, 3) = 1.5d0
    amat(4, 4) = 0.2d0
    amat(1, 3) = 0.3d0
    exact = [(2.0d0, 0.0d0), (-1.0d0, 0.0d0), (0.5d0, 0.0d0), (3.0d0, 0.0d0)]*(1 + rank)
    source = exact - matmul(amat, exact)
    solution = source
    call gmres_iterator(n, 4, 1.0d-12, 50, solution, rank, apply_iteration, &
                        original_residual, status, iterations, residual, .false.)
    if (status /= 0) error stop 'FAIL gmres status'
    if (maxval(abs(solution - exact)) > 1.0d-9*maxval(abs(exact))) error stop 'FAIL gmres oracle'
    if (residual > 1.0d-12 .or. iterations < 1 .or. iterations > 50) error stop 'FAIL gmres diagnostics'

    ! The projected fixed point (1,0) of L=diag(0,1), q=(1,0) violates 0=1:
    ! the equation gate must refuse it.
    projected_iteration = .true.
    amat = (0.0d0, 0.0d0)
    amat(1, 1) = 1.0d0
    source = (0.0d0, 0.0d0)
    source(1) = 1.0d0
    solution = source
    call gmres_iterator(n, 2, 1.0d-12, 20, solution, rank, apply_iteration, &
                        original_residual, status, iterations, residual, .false.)
    if (status == 0 .or. residual < 0.99d0) error stop 'FAIL projection accepted'
    projected_iteration = .false.

    call gmres_iterator(n, 0, 1.0d-12, 20, solution, rank, apply_iteration, &
                        original_residual, status, iterations, residual, .false.)
    if (status /= 4 .or. iterations /= 0) error stop 'FAIL invalid controls'

    call mpro%deinit()
    print *, 'All tests passed!'
contains
    subroutine apply_iteration(m, input, output)
        use arnoldi_mod, only: fzero, mode
        integer :: m
        complex(kind=kind(1d0)) :: input(m), output(m)
        output = matmul(amat, input)
        if (mode /= 2) output = output + fzero
        if (projected_iteration) output(1) = 0.5d0*output(1)
    end subroutine apply_iteration

    subroutine original_residual(m, candidate, value)
        use mpiprovider_module, only: mpro
        integer, intent(in) :: m
        complex(kind=kind(1d0)), intent(in) :: candidate(m)
        double precision, intent(out) :: value
        complex(kind=kind(1d0)) :: lhs(m)
        double precision :: buf(0:num_spec-1)
        lhs = candidate - matmul(amat, candidate)
        buf = 0.0d0
        buf(rank) = sum(abs(source - lhs))/max(sum(abs(source)) + sum(abs(lhs)), tiny(1.0d0))
        call mpro%allgather_inplace(buf)
        value = maxval(buf)
    end subroutine original_residual
end program test_gmres_iterator
