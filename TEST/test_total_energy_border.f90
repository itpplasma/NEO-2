program test_total_energy_border
    ! Two species with four unknowns each. L0 = (I-Pi) D (I-Pi) has the
    ! density and energy kets of each species as null vectors and their bras
    ! as left null vectors. An energy exchange of order one,
    ! nu*(w_1 - w_2)(u_1 - u_2), couples the species as D-T collisions do:
    ! u_a is the energy amplitude of species a (dual bra, blind to density)
    ! and w_a its density-free energy ket, so densities stay conserved and
    ! only the common energy ket K_E1 + K_E2 stays null. The exact solution
    ! has zero densities and zero total energy but opposite species energies.
    use arnoldi_mod, only: iterator, bordered_projection, &
        reduced_bordered_projection, energy_null_pair
    use collisionality_mod, only: num_spec
    use mpiprovider_module, only: mpro
    implicit none
    integer, parameter :: dp = kind(1.0d0), nper = 4, nspec = 2, n = nper*nspec
    real(dp), parameter :: nu = 0.3d0, rho = 30.0d0, relerr = 1.0d-8
    real(dp) :: kets(n, 2*nspec), bras(n, 2*nspec), eye(n, n), proj(n, n), diag(n, n)
    real(dp) :: lmat(n, n), pinv(n, n), work(n, n), velocity(nper), exact(n)
    real(dp) :: local_bras(n, 2*nspec), local_directions(n, 2*nspec)
    real(dp) :: shared_bras(n, nspec + 1), shared_directions(n, nspec + 1)
    real(dp) :: exchange(nspec, nspec), alpha(nspec), beta(nspec), ratio(nspec)
    real(dp) :: amplitude(nspec), residual, dual(n, 2*nspec), free(n, nspec)
    complex(dp) :: source(n), solution(n)
    integer :: status, iterations, i, j, a, b, energy_column(nspec)
    logical :: use_shared

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
        energy_column(i + 1) = 2*i + 2
    end do
    eye = 0.0d0
    diag = 0.0d0
    do i = 1, n
        eye(i, i) = 1.0d0
        diag(i, i) = 1.0d0 + real(i - 1, dp)/real(n - 1, dp)
    end do
    proj = eye - matmul(kets, matmul(inverse(matmul(transpose(bras), kets), 2*nspec), &
        transpose(bras)))
    lmat = matmul(proj, matmul(diag, proj))
    dual = matmul(bras, transpose(inverse(matmul(transpose(kets), bras), 2*nspec)))
    do a = 1, nspec
        free(:, a) = kets(:, energy_column(a)) - kets(:, 2*a - 1)* &
            sum(bras(:, 2*a - 1)*kets(:, energy_column(a)))/sum(bras(:, 2*a - 1)*kets(:, 2*a - 1))
    end do
    do a = 1, nspec
        do b = 1, nspec
            lmat = lmat + nu*sign(1.0d0, 0.5d0 - abs(a - b))* &
                outer(free(:, a), dual(:, energy_column(b)))
        end do
    end do
    do i = 1, n
        do j = 1, n
            work(i, j) = sin(1.7d0*i + 0.9d0*j**2)
        end do
    end do
    pinv = inverse(2.0d0*diag + rho*matmul(proj, matmul(work, eye - proj)), n)

    ! Exact solution: no density, opposite species energies, zero total.
    do i = 1, n
        exact(i) = 1.0d0 - 1.5d0*real(i - 1, dp)/real(n - 1, dp) + 0.3d0*sin(real(i - 1, dp))
    end do
    exact = matmul(proj, exact) + 0.2d0*(free(:, 1) - free(:, 2))
    source = matmul(lmat, exact)

    ! Holding the energy of each species is inconsistent with this exact
    ! solution; the bordered fixed point carries a multiplier residual.
    local_bras = bras
    local_directions = matmul(pinv, kets)
    use_shared = .false.
    solution = source
    call iterator(0, n, 4, relerr, 300, solution, mpro%getRank(), apply_iteration, &
        status, iterations, residual, original_residual)
    print *, 'per-species energy: status', status, ' residual', residual
    if (status == 0 .or. residual < 1.0d-4) then
        print *, 'FAIL: per-species energy constraints no longer leave a multiplier residual'
        error stop 1
    end if

    ! Total energy: weights from the exchange matrix of the assembled operator.
    do a = 1, nspec
        do b = 1, nspec
            exchange(a, b) = sum(bras(:, energy_column(a))*matmul(lmat, kets(:, energy_column(b))))
        end do
    end do
    call energy_null_pair(nspec, exchange, alpha, beta, ratio)
    print *, 'exchange singular values / largest:', ratio
    if (ratio(nspec) > 1.0d-12) then
        print *, 'FAIL: exchange matrix has no null pair'
        error stop 1
    end if
    shared_bras = 0.0d0
    shared_directions = 0.0d0
    do a = 1, nspec
        shared_bras(:, a) = bras(:, 2*a - 1)
        shared_directions(:, a) = kets(:, 2*a - 1)
        shared_bras(:, nspec + 1) = shared_bras(:, nspec + 1) + beta(a)*bras(:, energy_column(a))
        shared_directions(:, nspec + 1) = shared_directions(:, nspec + 1) + &
            alpha(a)*kets(:, energy_column(a))
    end do
    shared_directions = matmul(pinv, shared_directions)
    use_shared = .true.
    solution = source
    call iterator(0, n, 4, relerr, 300, solution, mpro%getRank(), apply_iteration, &
        status, iterations, residual, original_residual)
    print *, 'total energy: status', status, ' iterations', iterations, ' residual', residual
    if (status /= 0 .or. residual > relerr) then
        print *, 'FAIL: total-energy border did not converge'
        error stop 1
    end if
    if (maxval(abs(solution - exact)) > 1.0d-6*maxval(abs(exact))) then
        print *, 'FAIL: total-energy border differs from the exact solution'
        error stop 1
    end if
    do a = 1, nspec
        amplitude(a) = abs(sum(bras(:, energy_column(a))*solution))
    end do
    if (minval(amplitude) < 1.0d-3) then
        print *, 'FAIL: species energies were held at zero'
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

    function inverse(m, k) result(minv)
        integer, intent(in) :: k
        real(dp), intent(in) :: m(k, k)
        real(dp) :: minv(k, k), copy(k, k)
        integer :: piv(k), err, l
        copy = m
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
        if (use_shared) then
            call reduced_bordered_projection(nn, nspec + 1, shared_bras, shared_directions, &
                output, identity_sum)
        else
            call bordered_projection(nn, 2*nspec, local_bras, local_directions, output)
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
end program test_total_energy_border
