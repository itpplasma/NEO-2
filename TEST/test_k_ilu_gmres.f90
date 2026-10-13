program test_k_ilu_gmres
    use fortsparse, only: dp, sparse_solver_t, sparse_factor, sparse_solve, &
        sparse_free, fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_SUPERLU_ILU
    use fixed_point_gmres_mod, only: fixed_point_gmres_t
    implicit none
    integer, parameter :: n = 4
    type(sparse_solver_t) :: preconditioner
    type(fortsparse_status_t) :: status
    type(fixed_point_gmres_t) :: gm
    real(dp) :: a(n,n), c(n,n), p(n,n), full(n,n), q(n), f0(n), y(n)
    real(dp) :: want(n), reference(n), wrong(n), residual(n), pv(n)
    integer :: ipiv(n), info, j

    ! Deliberately P /= A: replacing A by P would solve the wrong equation.
    a = reshape([4.0_dp, -1.0_dp, 0.0_dp, 0.5_dp, &
        2.0_dp, 5.0_dp, -0.5_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 3.0_dp, -1.0_dp, &
        -0.5_dp, 0.0_dp, 2.0_dp, 6.0_dp], [n,n])
    c = 0.05_dp
    p = 0.0_dp
    do j = 1, n
        p(j,j) = real(j + 1, dp)
    end do
    pv = [2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    want = [1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp]
    q = matmul(a - c, want)
    full = a - c
    reference = q
    call dgesv(n, 1, full, n, ipiv, reference, n, info)
    if (info /= 0) error stop 'LAPACK reference failed'
    full = p - c
    wrong = q
    call dgesv(n, 1, full, n, ipiv, wrong, n, info)
    if (maxval(abs(wrong - want)) < 0.1_dp) error stop 'Weak wrong-system oracle'

    preconditioner%backend_id = FORTSPARSE_BACKEND_SUPERLU_ILU
    call sparse_factor(preconditioner, n, n, n, [1,2,3,4,5], [1,2,3,4], pv, status)
    if (.not. status_ok(status)) error stop 'ILU factor failed'
    call sparse_solve(preconditioner, q, f0, status)
    if (.not. status_ok(status)) error stop 'ILU initial source failed'
    call gm%start(f0, 1.0e-12_dp, 20)
    do while (gm%needs_apply())
        y = matmul(a, gm%vin) - matmul(c, gm%vin)
        call sparse_solve(preconditioner, y, status)
        if (.not. status_ok(status)) error stop 'ILU application failed'
        call gm%put(gm%vin - y)
    end do
    if (.not. gm%converged) error stop 'GMRES did not converge'
    if (maxval(abs(gm%x - reference)) > 1.0e-10_dp) &
        error stop 'ILU GMRES differs from LAPACK full-equation oracle'
    residual = q - matmul(a - c, gm%x)
    if (sum(abs(residual))/sum(abs(q)) > 1.0e-12_dp) &
        error stop 'Original full equation residual failed'
    call sparse_solve(preconditioner, residual, status)
    if (sum(abs(residual))/sum(abs(gm%x)) > 1.0e-12_dp) &
        error stop 'Independent preconditioned residual failed'
    print *, 'Full-equation ILU GMRES PASS, applications:', gm%napply
    call sparse_free(preconditioner)
end program
