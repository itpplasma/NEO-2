!> Test of arnoldi_mod::gmres_iterator on a dense fixed-point problem
!> f = A f + q whose iteration matrix A has the eigenvalue pattern of the
!> NEO-2-QL non-axisymmetric pass (LAG=6 AUG 30835): eigenvalues on and
!> outside the unit circle (1.0047-0.0062i, -1.008, 1.0072) and very close
!> to 1 (1.0000036+3e-5i), where plain or deflated Richardson iteration
!> diverges. Oracle: dense LU solve of (I - A) f = q with LAPACK zgesv.
module gmres_test_operator
  implicit none
  integer, parameter :: dp = kind(1d0)
  complex(dp), allocatable :: amat(:,:)
  double precision :: noise = 0d0
  logical :: corrupt_operator = .false.
  integer(8) :: nseed = 777_8
contains
  subroutine next_iteration(n, fold, fnew)
    use arnoldi_mod, only : fzero, mode
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    integer :: n
    complex(dp), dimension(n) :: fold, fnew
    integer :: i
    fnew = matmul(amat, fold)
    if (mode .ne. 2) fnew = fnew + fzero
    if (corrupt_operator) then
      if (mode == 2) fnew(1) = cmplx(ieee_value(0d0, ieee_quiet_nan), 0d0, dp)
    end if
    ! optional rounding-like noise, mimicking inexact sparse solves
    if (noise > 0d0) then
      do i = 1, n
        nseed = mod(16807_8*nseed, 2147483647_8)
        fnew(i) = fnew(i) + noise*sqrt(sum(abs(fold)**2)/n) &
          *(dble(nseed)/2147483647d0 - 0.5d0)
      end do
    end if
  end subroutine next_iteration
end module gmres_test_operator

program test_gmres_iterator
  use gmres_test_operator
  use arnoldi_mod, only : gmres_iterator
  use mpiprovider_module, only : mpro
  use collisionality_mod, only : num_spec
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, &
    ieee_positive_inf, ieee_quiet_nan
  implicit none

  integer, parameter :: n = 200
  complex(dp), parameter :: special(6) = [ (1.0047016532d0, -6.16d-3), &
    (1.0000036189d0, 3.0d-5), (0.98769d0, 0d0), (-0.74147d0, -2.49d-2), &
    (-1.0079527d0, 0d0), (1.00723d0, 0d0) ]
  complex(dp) :: smat(n,n), sinv(n,n), lam(n), q(n), fref(n), f(n), work(n,n)
  integer :: ipiv(n), info, i, j
  integer(8) :: seed
  logical :: ok
  logical :: converged
  double precision :: invalid_relerr(5), invalid_target(5)

  call mpro%init()
  num_spec = 1
  seed = 12345_8

  ! A = S diag(lam) S^-1, S = I + 0.3 R (non-normal)
  do j = 1, n
    do i = 1, n
      smat(i, j) = 0.3d0*cmplx(urand() - 0.5d0, urand() - 0.5d0, dp)/sqrt(dble(n))
    end do
    smat(j, j) = smat(j, j) + 1d0
  end do
  lam(1:6) = special
  do i = 7, n
    lam(i) = 0.6d0*urand()*exp(cmplx(0d0, 6.283185307d0*urand(), dp))
  end do
  sinv = (0d0, 0d0)
  do i = 1, n
    sinv(i, i) = 1d0
  end do
  work = smat
  call zgesv(n, n, work, n, ipiv, sinv, n, info)
  if (info /= 0) stop 'FAIL: zgesv S'
  do j = 1, n
    work(:, j) = smat(:, j)*lam(j)
  end do
  allocate(amat(n, n))
  amat = matmul(work, sinv)

  do i = 1, n
    q(i) = cmplx(urand() - 0.5d0, urand() - 0.5d0, dp)
  end do

  ! Oracle: (I - A) fref = q
  work = -amat
  do i = 1, n
    work(i, i) = work(i, i) + 1d0
  end do
  fref = q
  call zgesv(n, 1, work, n, ipiv, fref, n, info)
  if (info /= 0) stop 'FAIL: zgesv I-A'

  ok = .true.
  ! GMRES(60) converges within one cycle; GMRES(45) exercises restarts.
  ! I - A has eigenvalues of modulus ~3e-5, so the error bound for a
  ! residual of 1e-10 is ~1e-5 (cond ~ 1e5).
  call check(60, 300, 1d-10, 1d-10, 1d-5)
  call check(45, 600, 1d-10, 1d-10, 1d-5)
  ! Noisy operator: the target 1e-16 is below the noise floor, so the
  ! iteration has to stop by stagnation, still reaching relerr = 1e-9.
  noise = 1d-12
  call check(60, 300, 1d-9, 1d-16, 1d-4)

  ! A noiseless restarted problem need not halve its residual every cycle.
  ! (I-A) = diag(1,10) has the exact solution (1,0.1).
  noise = 0d0
  amat = (0d0, 0d0)
  amat(2, 2) = -9d0
  q = (0d0, 0d0)
  q(1:2) = 1d0
  fref = q
  fref(2) = 0.1d0
  call check(1, 100, 1d-10, 1d-10, 1d-9)

  ! The fixed-point step N(x) can amplify a residual. Return checked x.
  q(2) = 1d-8
  fref(2) = 1d-9
  call check(1, 1, 1d-7, 1d-10, 1d-7)

  ! A=I and nonzero q is inconsistent: report failure and keep a finite x.
  amat = (0d0, 0d0)
  do i = 1, n
    amat(i, i) = 1d0
  end do
  q(2) = 0d0
  f = q
  call gmres_iterator(n, 1, 1d-7, 1d-10, 5, f, 0, next_iteration, converged)
  if (converged) ok = .false.
  if (.not. all(ieee_is_finite(dble(f)))) ok = .false.
  if (.not. all(ieee_is_finite(aimag(f)))) ok = .false.

  ! A finite source followed by a nonfinite operator must fail explicitly.
  amat = (0d0, 0d0)
  corrupt_operator = .true.
  f = q
  call gmres_iterator(n, 1, 1d-7, 1d-10, 5, f, 0, next_iteration, converged)
  call check_failure('nonfinite operator')
  corrupt_operator = .false.

  ! Invalid tolerances cannot turn a failure sentinel into convergence.
  invalid_relerr = 1d-7
  invalid_target = 1d-10
  invalid_relerr(1) = ieee_value(0d0, ieee_positive_inf)
  invalid_relerr(2) = ieee_value(0d0, ieee_quiet_nan)
  invalid_target(3) = ieee_value(0d0, ieee_positive_inf)
  invalid_target(4) = -1d-10
  invalid_relerr(5) = 0d0
  do i = 1, size(invalid_relerr)
    f = q
    call gmres_iterator(n, 1, invalid_relerr(i), invalid_target(i), 5, &
      f, 0, next_iteration, converged)
    call check_failure('invalid tolerance')
  end do

  if (ok) then
    print *, 'All tests passed!'
  else
    print *, 'FAIL'
  end if
  call mpro%deinit(.false.)
  if (.not. ok) error stop 'FAIL: gmres_iterator regression'

contains

  subroutine check_failure(label)
    character(len=*), intent(in) :: label
    logical :: valid
    valid = .not. converged
    if (.not. all(ieee_is_finite(dble(f)))) valid = .false.
    if (.not. all(ieee_is_finite(aimag(f)))) valid = .false.
    print *, label, ': finite unconverged result = ', valid
    if (.not. valid) ok = .false.
  end subroutine check_failure

  subroutine check(mrestart, itermax, relerr, reltarget, maxerr)
    integer, intent(in) :: mrestart, itermax
    double precision, intent(in) :: relerr, reltarget, maxerr
    double precision :: err, residual
    logical :: achieved
    f = q
    call gmres_iterator(n, mrestart, relerr, reltarget, itermax, f, 0, &
      next_iteration, achieved)
    err = sqrt(sum(abs(f - fref)**2)/sum(abs(fref)**2))
    print '(a,i4,a,es10.3)', ' restart ', mrestart, ': relative error vs zgesv ', err
    if (.not. (err <= maxerr)) ok = .false.
    residual = sqrt(sum(abs(q + matmul(amat, f) - f)**2)) &
      /max(sqrt(sum(abs(f)**2)), tiny(1d0))
    if (.not. achieved) ok = .false.
    if (.not. (residual <= relerr)) ok = .false.
  end subroutine check

  double precision function urand()
    ! Park-Miller minimal standard generator (no integer overflow)
    seed = mod(16807_8*seed, 2147483647_8)
    urand = dble(seed)/2147483647d0
  end function urand

end program test_gmres_iterator
