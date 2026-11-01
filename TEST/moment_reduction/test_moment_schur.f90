module dense_schur_fixture
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use moment_schur_solver
    use moment_test_support
    implicit none
    integer, parameter :: dp = real64
    complex(dp), allocatable :: p(:, :), d(:, :), w(:, :), lmat(:, :), kmat(:, :), bmat(:, :), lu(:, :)
    integer, allocatable :: piv(:)
    real(dp) :: inverse_scale = 1.0_dp
    integer :: inverse_failure = 0
    interface
        subroutine zgetrf(m, n, a, lda, ipiv, info)
            import dp
            integer, intent(in) :: m, n, lda
            complex(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine
        subroutine zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            import dp
            character, intent(in) :: trans
            integer, intent(in) :: n, nrhs, lda, ldb, ipiv(*)
            complex(dp), intent(in) :: a(lda, *)
            complex(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine
        subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import dp
            integer, intent(in) :: n, nrhs, lda, ldb
            complex(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine
        subroutine zgeev(jobvl, jobvr, n, a, lda, eig, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
            import dp
            character, intent(in) :: jobvl, jobvr
            integer, intent(in) :: n, lda, ldvl, ldvr, lwork
            complex(dp), intent(inout) :: a(lda, *)
            complex(dp), intent(out) :: eig(*), vl(ldvl, *), vr(ldvr, *), work(*)
            real(dp), intent(out) :: rwork(*)
            integer, intent(out) :: info
        end subroutine
    end interface
contains
    subroutine setup(op, n, r, c)
        type(moment_schur_operator_t), intent(out) :: op
        integer, intent(in) :: n, r, c
        integer :: i
        if (allocated(p)) deallocate (p, d, w, lmat, kmat, bmat, lu, piv)
        allocate (p(n, n), d(n, r), w(r, n), lmat(n, n), kmat(n, c), bmat(c, n), lu(n, n), piv(n))
        p = 0; d = 0; w = 0; lmat = 0; kmat = 0; bmat = 0
        do i = 1, n
            p(i, i) = 1
        end do
        inverse_scale = 1
        inverse_failure = 0
        op%nfield = n; op%nmoment = r; op%nconstraint = c
        op%inverse_p => pinv
        op%feedback => feedback
        op%moments => moments
        op%original => original
        op%force => force
        op%constraints => constraints
    end subroutine
    subroutine factor()
        integer :: info, n
        n = size(p, 1)
        lu = p
        call zgetrf(n, n, lu, n, piv, info)
        call check(info == 0, 'fixture preconditioner is nonsingular')
        lmat = p - matmul(d, w)
    end subroutine
    subroutine pinv(x, y, info)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: info
        complex(dp) :: rhs(size(x), 1)
        if (inverse_failure == 1) then
            info = 17
            return
        end if
        if (inverse_failure == 2) then
            y = cmplx(ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
            info = 0
            return
        end if
        rhs(:, 1) = x
        call zgetrs('N', size(x), 1, lu, size(x), piv, rhs, size(x), info)
        y = inverse_scale*rhs(:, 1)
    end subroutine
    subroutine feedback(x, y, info)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: info
        y = matmul(d, x); info = 0
    end subroutine
    subroutine moments(x, y, info)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: info
        y = matmul(w, x); info = 0
    end subroutine
    subroutine original(x, y, info)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: info
        y = matmul(lmat, x); info = 0
    end subroutine
    subroutine force(x, y, info)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: info
        y = matmul(kmat, x); info = 0
    end subroutine
    subroutine constraints(x, y, info)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: info
        y = matmul(bmat, x); info = 0
    end subroutine
    function direct_solve(a, rhs) result(x)
        complex(dp), intent(in) :: a(:, :), rhs(:)
        complex(dp) :: x(size(rhs)), ac(size(rhs), size(rhs)), bc(size(rhs), 1)
        integer :: ipiv(size(rhs)), info
        ac = a; bc(:, 1) = rhs
        call zgesv(size(rhs), 1, ac, size(rhs), ipiv, bc, size(rhs), info)
        call check(info == 0, 'independent full LAPACK solve')
        x = bc(:, 1)
    end function
    function eigenvalues(a) result(eig)
        complex(dp), intent(in) :: a(:, :)
        complex(dp) :: eig(size(a, 1)), ac(size(a, 1), size(a, 1)), vl(1, 1), vr(1, 1), work(8*size(a, 1))
        real(dp) :: rwork(2*size(a, 1))
        integer :: info, n
        n = size(a, 1); ac = a
        call zgeev('N', 'N', n, ac, n, eig, vl, 1, vr, 1, work, size(work), rwork, info)
        call check(info == 0, 'independent full spectrum')
    end function
end module

program test_moment_schur
    use dense_schur_fixture
    implicit none
    type(moment_schur_operator_t) :: op
    type(moment_schur_options_t) :: opts
    type(moment_schur_report_t) :: rep
    complex(dp), allocatable :: q(:), f(:), truth(:), g(:), lam(:), reference(:), xmat(:, :), full_eig(:), small_eig(:)
    complex(dp), allocatable :: border(:, :), rhsborder(:)
    complex(dp) :: none(0), q2(2), f2(2), lam2(2), g2(2), saved
    real(dp) :: spectral_error, err, best
    integer :: i, j, n, r, c, info, match
    logical, allocatable :: taken(:)

    ! General complex, nonnormal, nonorthogonal feedback; singular mixing allowed.
    n = 12; r = 4; c = 0
    call setup(op, n, r, c)
    allocate (q(n), f(n), truth(n), reference(n), xmat(n, r), full_eig(n), small_eig(r), taken(n))
    do i = 1, n
        p(i, i) = cmplx(2 + i*0.1_dp, 0.2_dp, dp)
        truth(i) = cmplx(sin(real(i, dp)), cos(real(2*i, dp)), dp)
        do j = 1, r
            d(i, j) = cmplx(sin(real(i + 3*j, dp)), cos(real(2*i - j, dp)), dp)*0.12_dp
            w(j, i) = cmplx(cos(real(3*i + j, dp)), sin(real(i - j, dp)), dp)*0.09_dp
        end do
    end do
    d(:, r) = 0 ! C need not be invertible.
    call factor()
    q = matmul(lmat, truth)
    opts%rtol = 1.0e-12_dp
    call solve_moment_schur(op, q, none, f, none, opts, rep)
    call check(rep%status == schur_ok, 'complex reduced solve')
    reference = direct_solve(lmat, q)
    err = maxval(abs(f - reference))
    call check(err < 2.0e-11_dp, 'reconstruct all fields vs full solve')
    call check(maxval(abs(f - truth)) < 2.0e-11_dp, 'manufactured field')
    do j = 1, r
        call pinv(d(:, j), xmat(:, j), info)
    end do
    full_eig = eigenvalues(matmul(xmat, w))
    small_eig = eigenvalues(matmul(w, xmat))
    taken = .false.; spectral_error = 0
    do j = 1, r
        best = huge(1.0_dp); match = 0
        do i = 1, n
            if (taken(i)) cycle
            if (abs(full_eig(i) - small_eig(j)) < best) then
                best = abs(full_eig(i) - small_eig(j)); match = i
            end if
        end do
        call check(match > 0, 'spectral multiplicity matching')
        taken(match) = .true.
        spectral_error = max(spectral_error, best)
    end do
    call check(spectral_error < 1.0e-11_dp, 'nonzero spectra and multiplicities agree')
    call check(maxval(abs(full_eig), mask=.not. taken) < 1.0e-11_dp, 'only extra zero eigenvalues')
    print '(a,es14.6,a,es14.6)', 'full_vs_reduced_error=', err, ' spectrum_error=', spectral_error

    ! Explicit residual must reject an inaccurate sparse inverse.
    inverse_scale = 1.001_dp
    call solve_moment_schur(op, q, none, f, none, opts, rep)
    call check(rep%status == schur_reconstruction, 'inaccurate inverse cannot pass via reduced residual')
    call check(rep%original_residual > 1.0e-5_dp, 'original residual exposed')
    inverse_scale = 1
    inverse_failure = 1
    call solve_moment_schur(op, q, none, f, none, opts, rep)
    call check(rep%status == schur_callback, 'propagate inverse failure')
    inverse_failure = 2
    call solve_moment_schur(op, q, none, f, none, opts, rep)
    call check(rep%status == schur_nonfinite, 'reject nonfinite inverse output')
    inverse_failure = 0
    saved = q(1); q(1) = cmplx(ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
    call solve_moment_schur(op, q, none, f, none, opts, rep)
    call check(rep%status == schur_nonfinite, 'reject nonfinite RHS')
    q(1) = saved
    call solve_moment_schur(op, q(1:n - 1), none, f, none, opts, rep)
    call check(rep%status == schur_invalid, 'reject wrong RHS dimension')
    opts%rtol = -1
    call solve_moment_schur(op, q, none, f, none, opts, rep)
    call check(rep%status == schur_invalid, 'reject negative tolerance')
    opts%rtol = 1.0e-12_dp
    q = 0
    call solve_moment_schur(op, q, none, f, none, opts, rep)
    call check(rep%status == schur_ok .and. rep%iterations == 0, 'zero RHS without Arnoldi normalization')
    call check(maxval(abs(f)) < tiny(1.0_dp), 'zero RHS gives exactly zero')
    q = 1.0e-100_dp*matmul(lmat, truth)
    call solve_moment_schur(op, q, none, f, none, opts, rep)
    call check(rep%status == schur_ok .and. rep%iterations > 0, 'tiny nonzero RHS is not silently zeroed')
    call check(maxval(abs(f/1.0e-100_dp - truth)) < 2.0e-11_dp, 'scale-aware solve')
    deallocate (q, f, truth, reference, xmat, full_eig, small_eig, taken)

    ! Bad Richardson eigenvalues do not obstruct the reduced linear solve.
    call setup(op, 4, 4, 0)
    do i = 1, 4
        w(i, i) = 1
    end do
    d(1, 1) = 2; d(2, 2) = 3; d(3, 3) = -1.5_dp; d(4, 4) = 0.2_dp
    call factor()
    allocate (q(4), f(4), truth(4))
    truth = [(cmplx(real(i, dp), -0.5_dp*i, dp), i=1,4)]
    q = matmul(lmat, truth)
    call solve_moment_schur(op, q, none, f, none, opts, rep)
    call check(rep%status == schur_ok .and. maxval(abs(f - truth)) < 1.0e-10_dp, 'unstable Richardson feedback')
    opts%maxiter = 1
    call solve_moment_schur(op, q, none, f, none, opts, rep)
    call check(rep%status == schur_limit, 'iteration budget exhaustion rejected')
    opts%maxiter = 300; opts%restart = 1
    d = 0
    do i = 1, 4
        d(i, i) = 0.1_dp*i
    end do
    call factor()
    q = matmul(lmat, truth)
    call solve_moment_schur(op, q, none, f, none, opts, rep)
    call check(rep%status == schur_ok .and. rep%iterations > 1, 'actual restarted GMRES path')
    opts%restart = 30
    d = p
    call factor()
    q = 1
    call solve_moment_schur(op, q, none, f, none, opts, rep)
    call check(rep%status == schur_breakdown, 'incompatible singular system fails without division by zero')
    deallocate (q, f, truth)

    ! Include the regularizer in D W: P=M+R -> L=P-(C+R), not P-C.
    call setup(op, 2, 2, 0)
    p(1, 1) = 1.7_dp; p(1, 2) = 0.3_dp; p(2, 2) = 2
    d(1, 1) = 0.7_dp; d(1, 2) = 0.3_dp
    w(1, 1) = 1; w(2, 2) = 1
    call factor()
    q2 = [cmplx(1.0_dp, 0.2_dp, dp), cmplx(4.0_dp, -0.4_dp, dp)]
    call solve_moment_schur(op, q2, none, f2, none, opts, rep)
    call check(rep%status == schur_ok, 'regularizer-aware exact reduction')
    call check(maxval(abs(f2 - direct_solve(lmat, q2))) < 1.0e-11_dp, 'regularizer recovery')
    d = 0 ! Keep original L; incorrectly omit R only from the reduction.
    call solve_moment_schur(op, q2, none, f2, none, opts, rep)
    call check(rep%status == schur_reconstruction, 'omitting regularizer is detected by original residual')

    ! Exact null mode: solve a bordered system without dividing by 1-lambda.
    call setup(op, 2, 1, 1)
    d(1, 1) = 1; w(1, 1) = 1; kmat(1, 1) = 1; bmat(1, 1) = 1
    call factor()
    q2 = [cmplx(0.0_dp, 0.0_dp, dp), cmplx(3.0_dp, -1.0_dp, dp)]
    g2 = 0
    call solve_moment_schur(op, q2, g2(1:1), f2, lam2(1:1), opts, rep)
    call check(rep%status == schur_ok, 'compatible constrained nullspace')
    call check(abs(f2(1)) < 1.0e-12_dp .and. abs(f2(2) - q2(2)) < 1.0e-11_dp, 'null coordinate fixed')
    g2(1) = (2.0_dp, -0.5_dp)
    call solve_moment_schur(op, q2, g2(1:1), f2, lam2(1:1), opts, rep)
    call check(rep%status == schur_ok .and. abs(f2(1) - g2(1)) < 1.0e-11_dp, 'nonzero constraint data')
    q2(1) = 1; g2 = 0
    call solve_moment_schur(op, q2, g2(1:1), f2, lam2(1:1), opts, rep)
    call check(rep%status == schur_original_residual, 'bordered success is not original-equation success')
    call check(rep%bordered_residual < 1.0e-11_dp .and. rep%original_residual > 0.9_dp, 'multiplier forcing reported')

    ! Two exchanging species: total energy invariant, difference is NOT a gauge.
    call setup(op, 2, 2, 1)
    d(1, 2) = 1; d(2, 1) = 1; w(1, 1) = 1; w(2, 2) = 1
    kmat(:, 1) = 1; bmat(1, :) = 1
    call factor()
    q2 = [cmplx(1.0_dp, 0.0_dp, dp), cmplx(-1.0_dp, 0.0_dp, dp)]
    g2 = 0
    call solve_moment_schur(op, q2, g2(1:1), f2, lam2(1:1), opts, rep)
    call check(rep%status == schur_ok, 'only total conserved energy constrained')
    call check(maxval(abs(f2 - 0.5_dp*q2)) < 1.0e-11_dp, 'retain physical interspecies energy exchange')
    allocate (border(3, 3), rhsborder(3))
    border(1:2, 1:2) = lmat; border(1:2, 3) = kmat(:, 1)
    border(3, 1:2) = bmat(1, :); border(3, 3) = 0
    rhsborder(1:2) = q2; rhsborder(3) = 0
    reference = direct_solve(border, rhsborder)
    call check(maxval(abs(reference(1:2) - f2)) < 1.0e-11_dp, 'independent full bordered solve')
    deallocate (reference, border, rhsborder)
    call setup(op, 2, 2, 2)
    d(1, 2) = 1; d(2, 1) = 1; w(1, 1) = 1; w(2, 2) = 1
    kmat(1, 1) = 1; kmat(2, 2) = 1; bmat = kmat
    call factor()
    call solve_moment_schur(op, q2, g2, f2, lam2, opts, rep)
    call check(rep%status == schur_original_residual, 'overconstraining separate species energies rejected')
    call check(maxval(abs(lam2 - q2)) < 1.0e-11_dp, 'overconstraint multipliers exposed')

    ! No feedback: reduction has size zero and still checks the original equation.
    call setup(op, 2, 0, 0)
    call factor()
    call solve_moment_schur(op, q2, none, f2, none, opts, rep)
    call check(rep%status == schur_ok .and. rep%iterations == 0, 'empty feedback direct solve')
    call check(maxval(abs(f2 - q2)) < 1.0e-12_dp, 'empty feedback reconstruction')
    call finished('QL moment Schur solver')
end program
