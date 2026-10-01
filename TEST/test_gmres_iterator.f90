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
  integer(8) :: nseed = 777_8
contains
  subroutine next_iteration(n, fold, fnew)
    use arnoldi_mod, only : fzero, mode
    integer :: n
    complex(dp), dimension(n) :: fold, fnew
    integer :: i
    fnew = matmul(amat, fold)
    if (mode .ne. 2) fnew = fnew + fzero
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
  implicit none

  integer, parameter :: n = 200
  complex(dp), parameter :: special(6) = [ (1.0047016532d0, -6.16d-3), &
    (1.0000036189d0, 3.0d-5), (0.98769d0, 0d0), (-0.74147d0, -2.49d-2), &
    (-1.0079527d0, 0d0), (1.00723d0, 0d0) ]
  complex(dp) :: smat(n,n), sinv(n,n), lam(n), q(n), fref(n), f(n), work(n,n)
  integer :: ipiv(n), info, i, j
  integer(8) :: seed
  logical :: ok

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

  if (ok) then
    print *, 'All tests passed!'
  else
    print *, 'FAIL'
  end if
  call mpro%deinit(.false.)

contains

  subroutine check(mrestart, itermax, relerr, reltarget, maxerr)
    integer, intent(in) :: mrestart, itermax
    double precision, intent(in) :: relerr, reltarget, maxerr
    double precision :: err
    f = q
    call gmres_iterator(n, mrestart, relerr, reltarget, itermax, f, 0, next_iteration)
    err = sqrt(sum(abs(f - fref)**2)/sum(abs(fref)**2))
    print '(a,i4,a,es10.3)', ' restart ', mrestart, ': relative error vs zgesv ', err
    if (.not. (err <= maxerr)) ok = .false.
  end subroutine check

  double precision function urand()
    ! Park-Miller minimal standard generator (no integer overflow)
    seed = mod(16807_8*seed, 2147483647_8)
    urand = dble(seed)/2147483647d0
  end function urand

end program test_gmres_iterator
