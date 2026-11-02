program test_bordered_residual
    ! Two species with four unknowns each. L0 = (I-Pi) D (I-Pi) has the
    ! density and energy kets of each species as null vectors. An energy
    ! exchange of order one, nu*(w_1 - w_2)(u_1 - u_2), couples the species as
    ! D-T collisions do (u_a: energy amplitude of species a, w_a: its
    ! density-free energy ket). Holding each species' density and energy at
    ! zero then needs a multiplier: the bordered fixed point satisfies
    ! L f + K lambda = q with lambda /= 0, so the residual of L f = q stalls
    ! while the residual of the bordered system converges.
    use arnoldi_mod, only: iterator, bordered_projection, multiplier_force
    use collisionality_mod, only: num_spec
    use mpiprovider_module, only: mpro
    implicit none
    integer, parameter :: dp = kind(1.0d0), nper = 4, nspec = 2, n = nper*nspec, m = 2*nspec
    real(dp), parameter :: nu = 0.3d0, rho = 30.0d0, relerr = 1.0d-8
    real(dp) :: kets(n, m), bras(n, m), eye(n, n), proj(n, n), diag(n, n)
    real(dp) :: lmat(n, n), pinv(n, n), work(n, n), velocity(nper), exact(n)
    real(dp) :: directions(n, m), dual(n, m), free(n, nspec)
    real(dp) :: bordered(n + m, n + m), oracle(n + m, 1), residual
    complex(dp) :: source(n), solution(n), kforce(n)
    integer :: status, iterations, i, j, a, b, info, pivots(n + m)
    logical :: use_bordered_residual

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
    lmat = matmul(proj, matmul(diag, proj))
    dual = matmul(bras, transpose(inverse(matmul(transpose(kets), bras), m)))
    do a = 1, nspec
        free(:, a) = kets(:, 2*a) - kets(:, 2*a - 1)* &
            sum(bras(:, 2*a - 1)*kets(:, 2*a))/sum(bras(:, 2*a - 1)*kets(:, 2*a - 1))
    end do
    do a = 1, nspec
        do b = 1, nspec
            lmat = lmat + nu*sign(1.0d0, 0.5d0 - abs(a - b))*outer(free(:, a), dual(:, 2*b))
        end do
    end do
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
    exact = matmul(proj, exact) + 0.2d0*(free(:, 1) - free(:, 2))
    source = matmul(lmat, exact)

    bordered = 0.0d0
    bordered(1:n, 1:n) = lmat
    bordered(1:n, n + 1:) = kets
    bordered(n + 1:, 1:n) = transpose(bras)
    oracle = 0.0d0
    oracle(1:n, 1) = real(source, dp)
    call dgesv(n + m, 1, bordered, n + m, pivots, oracle, n + m, info)
    if (info /= 0) error stop 'bordered oracle singular'
    print *, 'bordered multiplier lambda:', oracle(n + 1:, 1)

    ! The residual of L f = q stops at the multiplier force K lambda.
    use_bordered_residual = .false.
    solution = source
    call iterator(0, n, 4, relerr, 300, solution, mpro%getRank(), apply_iteration, &
        status, iterations, residual, check_residual)
    print *, 'residual of L f = q: status', status, ' residual', residual
    if (status == 0 .or. residual < 1.0d-4) then
        print *, 'FAIL: the reference no longer stalls at the multiplier force'
        error stop 1
    end if

    ! The residual of the bordered system converges to the bordered solution.
    use_bordered_residual = .true.
    solution = source
    call iterator(0, n, 4, relerr, 300, solution, mpro%getRank(), apply_iteration, &
        status, iterations, residual, check_residual)
    print *, 'bordered residual: status', status, ' iterations', iterations, ' residual', residual
    if (status /= 0 .or. residual > relerr) then
        print *, 'FAIL: bordered residual did not converge'
        error stop 1
    end if
    if (maxval(abs(solution - oracle(1:n, 1))) > 1.0d-6*maxval(abs(oracle(1:n, 1)))) then
        print *, 'FAIL: solution differs from the bordered oracle'
        error stop 1
    end if
    call multiplier_force(n, m, kets, source - matmul(lmat, solution), kforce, identity_sum)
    if (maxval(abs(kforce - matmul(kets, oracle(n + 1:, 1)))) > &
        1.0d-6*maxval(abs(matmul(kets, oracle(n + 1:, 1))))) then
        print *, 'FAIL: multiplier force differs from K lambda of the oracle'
        error stop 1
    end if

    call mpro%deinit()
    print *, 'All tests passed!'

contains

    function outer(u, v) result(w)
        real(dp), intent(in) :: u(:), v(:)
        real(dp) :: w(size(u), size(v))
        integer :: k
        do k = 1, size(v)
            w(:, k) = u*v(k)
        end do
    end function outer

    function inverse(mat, k) result(minv)
        integer, intent(in) :: k
        real(dp), intent(in) :: mat(k, k)
        real(dp) :: minv(k, k), copy(k, k)
        integer :: piv(k), err, l
        copy = mat
        minv = 0.0d0
        do l = 1, k
            minv(l, l) = 1.0d0
        end do
        call dgesv(k, k, copy, k, piv, minv, k, err)
        if (err /= 0) error stop 'singular test matrix'
    end function inverse

    subroutine identity_sum(k, x)
        integer, intent(in) :: k
        complex(dp), intent(inout) :: x(k)
    end subroutine identity_sum

    subroutine apply_iteration(nn, input, output)
        use arnoldi_mod, only: fzero, mode
        integer :: nn
        complex(dp) :: input(nn), output(nn)
        complex(dp) :: rhs(nn)
        rhs = -matmul(lmat, input)
        if (mode /= 2) rhs = rhs + fzero
        output = input + matmul(pinv, rhs)
        call bordered_projection(nn, m, bras, directions, output)
    end subroutine apply_iteration

    subroutine check_residual(nn, candidate, value)
        integer, intent(in) :: nn
        complex(dp), intent(in) :: candidate(nn)
        real(dp), intent(out) :: value
        complex(dp) :: lhs(nn), r(nn), force(nn)
        lhs = matmul(lmat, candidate)
        r = source - lhs
        if (use_bordered_residual) then
            call multiplier_force(nn, m, kets, r, force, identity_sum)
            r = r - force
        end if
        value = sum(abs(r))/(sum(abs(source)) + sum(abs(lhs)))
    end subroutine check_residual
end program test_bordered_residual
