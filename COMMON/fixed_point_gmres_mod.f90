!> Restarted GMRES for the fixed-point problem x = f0 + M x, i.e.
!> (I - M) x = f0, with reverse communication for the operator M.
!>
!> The stopping test is the one of the Richardson (fixed-point) iteration
!> x_new = f0 + M x it replaces: sum(abs(r)) < eps * sum(abs(x)) with
!> r = f0 + M x - x. It is evaluated in the 1-norm on the GMRES residual
!> after every Arnoldi step (vector work only, no operator application).
!> A projected convergence estimate always requests a true residual check
!> before convergence is reported. That check counts against maxapply.
!> The converged result is the iterate whose true residual was checked.
!>
!> Usage:
!>   call gm%start(f0, eps, maxapply)
!>   do while (gm%needs_apply())
!>     y = M gm%vin                 ! caller applies the operator
!>     call gm%put(y)
!>   end do
!>   x = gm%x ; if (.not. gm%converged) warn
!>
!> maxapply bounds the number of operator applications.
module fixed_point_gmres_mod
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter, public :: fixed_point_gmres_default_restart = 30

  integer, parameter :: PHASE_DONE = 0, PHASE_RESID = 1, PHASE_ARNOLDI = 2

  type, public :: fixed_point_gmres_t
     integer  :: n = 0, nrestart = 0, maxapply = 0
     integer  :: napply = 0      !< operator applications so far
     integer  :: j = 0           !< Arnoldi step within the current cycle
     integer  :: phase = PHASE_DONE
     logical  :: converged = .false.
     logical  :: failed = .false.
     real(dp) :: eps = 0.0_dp
     real(dp) :: resid_rel = huge(1.0_dp) !< last sum|r| / sum|x|
     real(dp), allocatable :: f0(:), x(:), vin(:)
     real(dp), allocatable :: v(:,:), h(:,:), cs(:), sn(:), g(:), t(:)
   contains
     procedure :: start => fpg_start
     procedure :: needs_apply => fpg_needs_apply
     procedure :: put => fpg_put
     procedure :: free => fpg_free
  end type fixed_point_gmres_t

contains

  subroutine fpg_start(self, f0, eps, maxapply, nrestart)
    class(fixed_point_gmres_t), intent(inout) :: self
    real(dp), intent(in) :: f0(:)
    real(dp), intent(in) :: eps
    integer,  intent(in) :: maxapply
    integer,  intent(in), optional :: nrestart

    integer :: n, m

    n = size(f0)
    m = fixed_point_gmres_default_restart
    if (present(nrestart)) m = nrestart
    m = max(1, min(m, maxapply, n))

    if (self%n /= n .or. self%nrestart /= m) then
       call self%free()
       allocate(self%f0(n), self%x(n), self%vin(n), self%t(n))
       allocate(self%v(n, m+1), self%h(m+1, m), self%cs(m), self%sn(m), &
            self%g(m+1))
       self%n = n
       self%nrestart = m
    end if

    self%eps = eps
    self%maxapply = maxapply
    self%napply = 0
    self%j = 0
    self%converged = .false.
    self%failed = .false.
    self%resid_rel = huge(1.0_dp)
    self%f0 = f0
    self%x = f0
    self%vin = self%x
    self%phase = PHASE_RESID
    if (maxapply < 1) self%phase = PHASE_DONE
    if (.not. all(ieee_is_finite(f0)) .or. .not. ieee_is_finite(eps)) then
       self%x = 0.0_dp
       call mark_failure(self)
    else if (eps < 0.0_dp) then
       call mark_failure(self)
    end if
  end subroutine fpg_start

  logical function fpg_needs_apply(self)
    class(fixed_point_gmres_t), intent(in) :: self
    fpg_needs_apply = self%phase /= PHASE_DONE
  end function fpg_needs_apply

  !> Hand in y = M vin and advance to the next request.
  subroutine fpg_put(self, y)
    class(fixed_point_gmres_t), intent(inout) :: self
    real(dp), intent(in) :: y(:)

    self%napply = self%napply + 1
    if (.not. all(ieee_is_finite(y))) then
       call mark_failure(self)
       return
    end if
    select case (self%phase)
    case (PHASE_RESID)
       call true_residual_step(self, y)
    case (PHASE_ARNOLDI)
       call arnoldi_step(self, y)
    end select
  end subroutine fpg_put

  subroutine true_residual_step(self, y)
    type(fixed_point_gmres_t), intent(inout) :: self
    real(dp), intent(in) :: y(:)

    real(dp) :: rnorm1, xnorm1, beta

    ! Residual r = f0 + M x - x, kept in t.
    self%t = self%f0 + y - self%x
    if (.not. all(ieee_is_finite(self%t))) then
       call mark_failure(self)
       return
    end if
    rnorm1 = sum(abs(self%t))
    xnorm1 = sum(abs(self%x))
    if (.not. ieee_is_finite(rnorm1) .or. .not. ieee_is_finite(xnorm1)) then
       call mark_failure(self)
       return
    end if
    self%resid_rel = rnorm1 / max(xnorm1, tiny(1.0_dp))
    self%converged = rnorm1 <= self%eps * xnorm1

    if (self%converged .or. self%napply >= self%maxapply) then
       if (.not. self%converged) then
          self%vin = self%x + self%t
          if (.not. all(ieee_is_finite(self%vin))) then
             call mark_failure(self)
             return
          end if
          self%x = self%vin
       end if
       self%phase = PHASE_DONE
       return
    end if

    beta = norm2(self%t)
    self%v(:, 1) = self%t / beta
    self%g = 0.0_dp
    self%g(1) = beta
    self%h = 0.0_dp
    self%j = 0
    self%vin = self%v(:, 1)
    self%phase = PHASE_ARNOLDI
  end subroutine true_residual_step

  subroutine arnoldi_step(self, y)
    type(fixed_point_gmres_t), intent(inout) :: self
    real(dp), intent(in) :: y(:)

    integer  :: i, j
    real(dp) :: hij, denom, temp, wnorm
    logical  :: breakdown, estimated

    self%j = self%j + 1
    j = self%j

    ! w = (I - M) v_j, orthogonalised by modified Gram-Schmidt into t.
    self%t = self%v(:, j) - y
    if (.not. all(ieee_is_finite(self%t))) then
       call mark_failure(self)
       return
    end if
    wnorm = norm2(self%t)
    do i = 1, j
       hij = dot_product(self%t, self%v(:, i))
       self%h(i, j) = hij
       self%t = self%t - hij * self%v(:, i)
    end do
    self%h(j+1, j) = norm2(self%t)
    breakdown = self%h(j+1, j) <= epsilon(1.0_dp) * wnorm
    self%v(:, j+1) = 0.0_dp
    if (.not. breakdown) self%v(:, j+1) = self%t / self%h(j+1, j)

    ! Givens rotations: keep h(1:j,1:j) upper triangular.
    do i = 1, j - 1
       temp = self%cs(i) * self%h(i, j) + self%sn(i) * self%h(i+1, j)
       self%h(i+1, j) = -self%sn(i) * self%h(i, j) + self%cs(i) * self%h(i+1, j)
       self%h(i, j) = temp
    end do
    denom = hypot(self%h(j, j), self%h(j+1, j))
    if (denom > 0.0_dp) then
       self%cs(j) = self%h(j, j) / denom
       self%sn(j) = self%h(j+1, j) / denom
    else
       self%cs(j) = 1.0_dp
       self%sn(j) = 0.0_dp
    end if
    self%h(j, j) = denom
    self%h(j+1, j) = 0.0_dp
    self%g(j+1) = -self%sn(j) * self%g(j)
    self%g(j) = self%cs(j) * self%g(j)

    estimated = gmres_resid_small(self)
    if (self%failed) return
    if (estimated) then
       ! The estimate can drift for an inexact operator application.
       ! Verify the actual iterate without exceeding the solve budget.
       self%x = self%t
       if (self%napply >= self%maxapply) then
          self%phase = PHASE_DONE
       else
          self%vin = self%x
          self%phase = PHASE_RESID
       end if
    else if (self%napply >= self%maxapply) then
       ! Budget exhausted: return the GMRES iterate, flagged unconverged.
       self%x = self%t
       self%phase = PHASE_DONE
    else if (breakdown .or. j == self%nrestart) then
       self%x = self%t
       self%vin = self%x
       self%phase = PHASE_RESID
    else
       self%vin = self%v(:, j+1)
    end if
  end subroutine arnoldi_step

  !> Fixed-point test on the current GMRES iterate in the 1-norm.
  !> The GMRES residual is r_j = V_{j+1} Q_j^T (g(j+1) e_{j+1}).
  !> On return t holds the current iterate x_j.
  logical function gmres_resid_small(self)
    type(fixed_point_gmres_t), intent(inout) :: self

    real(dp) :: z(self%j+1), temp, rnorm1, xnorm1
    integer  :: i, j

    j = self%j
    z = 0.0_dp
    z(j+1) = self%g(j+1)
    do i = j, 1, -1
       temp = self%cs(i) * z(i) - self%sn(i) * z(i+1)
       z(i+1) = self%sn(i) * z(i) + self%cs(i) * z(i+1)
       z(i) = temp
    end do
    self%t = matmul(self%v(:, 1:j+1), z)
    rnorm1 = sum(abs(self%t))

    self%t = self%x
    call update_solution(self, self%t)
    xnorm1 = sum(abs(self%t))
    if (.not. all(ieee_is_finite(self%t)) .or. &
         .not. ieee_is_finite(rnorm1) .or. .not. ieee_is_finite(xnorm1)) then
       gmres_resid_small = .false.
       call mark_failure(self)
       return
    end if

    self%resid_rel = rnorm1 / max(xnorm1, tiny(1.0_dp))
    gmres_resid_small = rnorm1 <= self%eps * xnorm1
  end function gmres_resid_small

  !> xout = xout + V_j H_j^{-1} g_j (xout holds the cycle start on entry).
  subroutine update_solution(self, xout)
    type(fixed_point_gmres_t), intent(in) :: self
    real(dp), intent(inout) :: xout(:)

    real(dp) :: c(self%j)
    integer  :: i, j

    j = self%j
    do i = j, 1, -1
       c(i) = self%g(i) - dot_product(self%h(i, i+1:j), c(i+1:j))
       if (self%h(i, i) /= 0.0_dp) then
          c(i) = c(i) / self%h(i, i)
       else
          c(i) = 0.0_dp
       end if
    end do
    xout = xout + matmul(self%v(:, 1:j), c)
  end subroutine update_solution

  subroutine mark_failure(self)
    type(fixed_point_gmres_t), intent(inout) :: self

    self%converged = .false.
    self%failed = .true.
    self%resid_rel = huge(1.0_dp)
    self%phase = PHASE_DONE
  end subroutine mark_failure

  subroutine fpg_free(self)
    class(fixed_point_gmres_t), intent(inout) :: self
    if (allocated(self%f0)) deallocate(self%f0, self%x, self%vin, self%t)
    if (allocated(self%v)) deallocate(self%v, self%h, self%cs, self%sn, self%g)
    self%n = 0
    self%nrestart = 0
    self%phase = PHASE_DONE
  end subroutine fpg_free

end module fixed_point_gmres_mod
