program test_bordered_conservation
    ! Two species with four unknowns each. Density and energy kets K span
    ! the null space of L0 = (I-Pi) D (I-Pi), Pi = K (B'K)^-1 B'. A weak
    ! exchange eta*E couples the species, so K is only nearly null for
    ! L = L0 + eta*E, as in the coupled NEO-2-QL operator. The bordered
    ! solution of [L K; B' 0] [f; lambda] = [q; 0] is the oracle.
    use arnoldi_mod, only: iterator, bordered_projection
    use collisionality_mod, only: num_spec
    use mpiprovider_module, only: mpro
    implicit none
    integer, parameter :: dp = kind(1.0d0), nper = 4, nspec = 2, n = nper*nspec, m = 2*nspec
    real(dp), parameter :: eta = 1.0d-6, rho = 30.0d0, relerr = 1.0d-5
    real(dp) :: kets(n, m), bras(n, m), eye(n, n), proj(n, n), diag(n, n), coupling(n, n)
    real(dp) :: lmat(n, n), pinv(n, n), work(n, n), directions(n, m), exact(n)
    real(dp) :: bordered(n + m, n + m), oracle(n + m, 1), velocity(nper)
    complex(dp) :: source(n), solution(n)
    real(dp) :: residual
    integer :: status, iterations, i, j, info, pivots(n + m)
    logical :: use_bordered

    call mpro%init()
    num_spec = mpro%getNumProcs()

    velocity = [0.5d0, 1.0d0, 1.5d0, 2.0d0]
    kets = 0.0d0
    bras = 0.0d0
    do i = 0, nspec - 1
        kets(i*nper + 1:(i + 1)*nper, 2*i + 1) = exp(-velocity**2)
        kets(i*nper + 1:(i + 1)*nper, 2*i + 2) = (velocity**2 - 1.5d0)*exp(-velocity**2)
        bras(i*nper + 1:(i + 1)*nper, 2*i + 1) = 1.0d0
        bras(i*nper + 1:(i + 1)*nper, 2*i + 2) = velocity**2 - 1.4d0
    end do
    eye = 0.0d0
    diag = 0.0d0
    do i = 1, n
        eye(i, i) = 1.0d0
        diag(i, i) = 1.0d0 + real(i - 1, dp)/real(n - 1, dp)
    end do
    proj = eye - matmul(kets, matmul(inverse(matmul(transpose(bras), kets), m), transpose(bras)))
    coupling = 0.0d0
    do i = 1, nper
        do j = 1, nper
            coupling(i, nper + j) = 1.0d0 - 2.0d0*real(j - 1, dp)/real(nper - 1, dp)
            coupling(nper + i, j) = 1.0d0 - 2.0d0*real(i - 1, dp)/real(nper - 1, dp)
        end do
    end do
    lmat = matmul(proj, matmul(diag, proj)) + eta*coupling
    ! Preconditioner P = 2D + rho (I-Pi) Z Pi maps the kets far outside
    ! their own span, like the regularized NEO-2 matrix.
    do i = 1, n
        do j = 1, n
            work(i, j) = sin(1.7d0*i + 0.9d0*j**2)
        end do
    end do
    pinv = inverse(2.0d0*diag + rho*matmul(proj, matmul(work, eye - proj)), n)
    directions = matmul(pinv, kets)

    do i = 1, n
        exact(i) = 1.0d0 - 1.5d0*real(i - 1, dp)/real(n - 1, dp) + 0.3d0*sin(real(i - 1, dp))
    end do
    source = matmul(lmat, exact)

    bordered = 0.0d0
    bordered(1:n, 1:n) = lmat
    bordered(1:n, n + 1:) = kets
    bordered(n + 1:, 1:n) = transpose(bras)
    oracle = 0.0d0
    oracle(1:n, 1) = real(source, dp)
    call dgesv(n + m, 1, bordered, n + m, pivots, oracle, n + m, info)
    if (info /= 0) error stop 'bordered oracle singular'

    use_bordered = .true.
    solution = source
    call iterator(0, n, 2, relerr, 200, solution, mpro%getRank(), apply_iteration, &
        status, iterations, residual, original_residual)
    print *, 'bordered: status', status, ' iterations', iterations, ' residual', residual
    if (status /= 0 .or. residual > relerr) then
        print *, 'FAIL: bordered iteration did not converge'
        error stop 1
    end if
    ! The bordered step contracts by 0.5, so the stopped iterate is near 2^-iterations.
    if (maxval(abs(solution - oracle(1:n, 1))) > 1.0d-4*maxval(abs(oracle(1:n, 1)))) then
        print *, 'FAIL: bordered iteration differs from the bordered oracle'
        error stop 1
    end if
    if (maxval(abs(matmul(transpose(bras), solution))) > 1.0d-10) then
        print *, 'FAIL: constraints violated'
        error stop 1
    end if

    ! The previous damped filter along K settles at an original-equation
    ! residual of 2.0e-4 for this system and must not be accepted.
    use_bordered = .false.
    solution = source
    call iterator(0, n, 2, relerr, 200, solution, mpro%getRank(), apply_iteration, &
        status, iterations, residual, original_residual)
    print *, 'damped filter: status', status, ' iterations', iterations, ' residual', residual
    if (status == 0 .or. residual < 1.0d-4) then
        print *, 'FAIL: damped filter reference no longer stalls'
        error stop 1
    end if

    call mpro%deinit()
    print *, 'All tests passed!'

contains

    function inverse(a, k) result(ainv)
        integer, intent(in) :: k
        real(dp), intent(in) :: a(k, k)
        real(dp) :: ainv(k, k), copy(k, k)
        integer :: piv(k), err, l
        copy = a
        ainv = 0.0d0
        do l = 1, k
            ainv(l, l) = 1.0d0
        end do
        call dgesv(k, k, copy, k, piv, ainv, k, err)
        if (err /= 0) error stop 'singular test matrix'
    end function inverse

    subroutine apply_iteration(nn, input, output)
        use arnoldi_mod, only: fzero, mode
        integer :: nn
        complex(dp) :: input(nn), output(nn)
        complex(dp) :: rhs(nn)
        integer :: c
        rhs = -matmul(lmat, input)
        if (mode /= 2) rhs = rhs + fzero
        output = input + matmul(pinv, rhs)
        if (use_bordered) then
            call bordered_projection(nn, m, bras, directions, output)
        else
            do c = 1, m
                output = output - 0.5d0*kets(:, c)*sum(bras(:, c)*output)/sum(bras(:, c)*kets(:, c))
            end do
        end if
    end subroutine apply_iteration

    subroutine original_residual(nn, candidate, value)
        integer, intent(in) :: nn
        complex(dp), intent(in) :: candidate(nn)
        real(dp), intent(out) :: value
        complex(dp) :: lhs(nn)
        lhs = matmul(lmat, candidate)
        value = sum(abs(source - lhs))/(sum(abs(source)) + sum(abs(lhs)))
    end subroutine original_residual
end program test_bordered_conservation
