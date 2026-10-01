!> Tests fixed_point_gmres_mod on x = f0 + M x with a slowly contracting,
!> non-normal dense M (Richardson contraction ~0.97, like the PAR integral
!> part). Oracle: direct LAPACK solve of (I - M) x = f0.
program test_fixed_point_gmres
  use fixed_point_gmres_mod, only : fixed_point_gmres_t
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: n = 300
  real(dp) :: m(n, n), a(n, n), f0(n), xref(n), xr(n), xnew(n)
  integer  :: ipiv(n), info, i, j, nrich
  integer(8) :: seed
  type(fixed_point_gmres_t) :: gm
  logical :: ok

  ok = .true.
  seed = 12345_8

  ! M = diag(lambda) + strictly upper non-normal part, lambda in [-0.4, 0.97].
  m = 0.0_dp
  do i = 1, n
     m(i, i) = -0.4_dp + 1.37_dp * real(i - 1, dp) / real(n - 1, dp)
     do j = i + 1, n
        m(i, j) = 0.3_dp * (lcg(seed) - 0.5_dp) / sqrt(real(n, dp))
     end do
  end do
  do i = 1, n
     f0(i) = 1.0_dp + lcg(seed)
  end do

  a = -m
  do i = 1, n
     a(i, i) = a(i, i) + 1.0_dp
  end do
  xref = f0
  call dgesv(n, 1, a, n, ipiv, xref, n, info)
  call check(info == 0, 'LAPACK reference solve')

  ! Richardson with the legacy stopping rule, for the solve count.
  xr = f0
  do nrich = 1, 100000
     xnew = f0 + matmul(m, xr)
     if (sum(abs(xnew - xr)) < 1.0e-10_dp * sum(abs(xr))) exit
     xr = xnew
  end do

  ! 1. Converges to the direct solution with far fewer applications.
  call run(gm, 1.0e-10_dp, 1000, 30)
  call check(gm%converged, 'GMRES converged')
  call check(maxval(abs(gm%x - xref)) / maxval(abs(xref)) < 1.0e-8_dp, &
       'GMRES agrees with direct solve to 1e-8')
  call check(gm%resid_rel < 1.0e-10_dp, 'reported residual below eps')
  call check(4 * gm%napply < nrich, 'GMRES needs < 1/4 of Richardson solves')
  print '(a,i0,a,i0)', ' applications: GMRES ', gm%napply, ', Richardson ', nrich

  ! 2. Short restart still converges (more applications).
  call run(gm, 1.0e-10_dp, 5000, 5)
  call check(gm%converged, 'GMRES(5) converged')
  call check(maxval(abs(gm%x - xref)) / maxval(abs(xref)) < 1.0e-8_dp, &
       'GMRES(5) agrees with direct solve')

  ! 3. Budget exhaustion is reported, never silent.
  call run(gm, 1.0e-14_dp, 7, 30)
  call check(.not. gm%converged, 'budget exhaustion flagged')
  call check(gm%napply == 7, 'budget respected exactly')

  ! 4. Zero right-hand side converges immediately to zero.
  call gm%start([(0.0_dp, i = 1, n)], 1.0e-10_dp, 10)
  do while (gm%needs_apply())
     call gm%put(matmul(m, gm%vin))
  end do
  call check(gm%converged .and. gm%napply == 1 .and. all(gm%x == 0.0_dp), &
       'zero source')

  ! 5. Strongly contracting M (few Richardson steps, like unit-vector
  !    propagator sources): GMRES must not need more applications.
  m = 0.05_dp * m
  xr = f0
  do nrich = 1, 100000
     xnew = f0 + matmul(m, xr)
     if (sum(abs(xnew - xr)) < 1.0e-10_dp * sum(abs(xr))) exit
     xr = xnew
  end do
  call run(gm, 1.0e-10_dp, 1000, 30)
  call check(gm%converged .and. gm%napply <= nrich, &
       'GMRES applications <= Richardson for fast contraction')
  print '(a,i0,a,i0)', ' applications: GMRES ', gm%napply, ', Richardson ', nrich

  if (ok) then
     print *, 'All tests passed!'
  else
     print *, 'FAIL'
     error stop 1
  end if

contains

  subroutine run(gm, eps, maxapply, nrestart)
    type(fixed_point_gmres_t), intent(inout) :: gm
    real(dp), intent(in) :: eps
    integer,  intent(in) :: maxapply, nrestart

    call gm%start(f0, eps, maxapply, nrestart)
    do while (gm%needs_apply())
       call gm%put(matmul(m, gm%vin))
    end do
  end subroutine run

  subroutine check(cond, what)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: what
    if (cond) then
       print '(2a)', ' PASS: ', what
    else
       print '(2a)', ' FAIL: ', what
       ok = .false.
    end if
  end subroutine check

  real(dp) function lcg(s)
    integer(8), intent(inout) :: s
    s = modulo(16807_8 * s, 2147483647_8)
    lcg = real(s, dp) / 2147483647.0_dp
  end function lcg

end program test_fixed_point_gmres
