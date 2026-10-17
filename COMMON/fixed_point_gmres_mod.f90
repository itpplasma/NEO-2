!> Restarted GMRES(m) for the fixed point x = b + T x, i.e. (I - T) x = b,
!> with reverse communication for the operator T:
!>
!>   call gm%start(b, eps, maxapply, m)
!>   do while (gm%needs_apply())
!>     call gm%put(T gm%vin)
!>   end do
!>
!> The initial guess is b. Each cycle starts from the true residual
!> r = b + T x - x and ends after m Arnoldi steps, on breakdown, or when the
!> GMRES residual passes the Richardson test sum|r| <= eps sum|x|. Converged
!> is set only by a true residual. maxapply bounds all applications of T.
module fixed_point_gmres_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  integer, parameter :: dp = kind(1.0d0)

  type, public :: fixed_point_gmres_t
     integer  :: m = 0, maxapply = 0, napply = 0, j = 0
     logical  :: arnoldi = .false., done = .true., converged = .false.
     real(dp) :: eps = 0.0_dp
     real(dp), allocatable :: b(:), x(:), vin(:), w(:)
     real(dp), allocatable :: v(:,:), h(:,:), cs(:), sn(:), g(:)
   contains
     procedure :: start => fpg_start
     procedure :: needs_apply => fpg_needs_apply
     procedure :: put => fpg_put
  end type fixed_point_gmres_t

contains

  subroutine fpg_start(self, b, eps, maxapply, m)
    class(fixed_point_gmres_t), intent(inout) :: self
    real(dp), intent(in) :: b(:), eps
    integer, intent(in) :: maxapply, m

    integer :: n

    n = size(b)
    self%m = max(1, min(m, maxapply, n))
    if (allocated(self%v)) deallocate(self%v, self%h, self%cs, self%sn, self%g)
    allocate(self%v(n, self%m+1), self%h(self%m+1, self%m), self%cs(self%m), &
             self%sn(self%m), self%g(self%m+1))
    self%b = b
    self%x = b
    self%vin = b
    self%w = b
    self%eps = eps
    self%maxapply = maxapply
    self%napply = 0
    self%arnoldi = .false.
    self%converged = .false.
    self%done = maxapply < 1
  end subroutine fpg_start

  logical function fpg_needs_apply(self)
    class(fixed_point_gmres_t), intent(in) :: self
    fpg_needs_apply = .not. self%done
  end function fpg_needs_apply

  !> Hand in y = T vin and advance to the next request.
  subroutine fpg_put(self, y)
    class(fixed_point_gmres_t), intent(inout) :: self
    real(dp), intent(in) :: y(:)

    self%napply = self%napply + 1
    if (self%arnoldi) then
       call arnoldi_step(self, y)
    else
       call restart(self, y)
    end if
  end subroutine fpg_put

  !> y = T x: test the true residual and start a new cycle from it.
  subroutine restart(self, y)
    type(fixed_point_gmres_t), intent(inout) :: self
    real(dp), intent(in) :: y(:)

    real(dp) :: beta

    self%w = self%b + y - self%x
    self%converged = sum(abs(self%w)) <= self%eps * sum(abs(self%x))
    beta = norm2(self%w)
    self%done = self%converged .or. self%napply >= self%maxapply &
         .or. .not. ieee_is_finite(beta)
    if (self%done) return

    self%v(:, 1) = self%w / beta
    self%g = 0.0_dp
    self%g(1) = beta
    self%j = 0
    self%vin = self%v(:, 1)
    self%arnoldi = .true.
  end subroutine restart

  !> y = T v_j: extend the Arnoldi basis and the Givens QR of H.
  subroutine arnoldi_step(self, y)
    type(fixed_point_gmres_t), intent(inout) :: self
    real(dp), intent(in) :: y(:)

    integer  :: i, j
    real(dp) :: temp, wnorm, z(self%m+1)
    logical  :: breakdown, small

    self%j = self%j + 1
    j = self%j

    ! w = (I - T) v_j, modified Gram-Schmidt.
    self%w = self%v(:, j) - y
    wnorm = norm2(self%w)
    do i = 1, j
       self%h(i, j) = dot_product(self%w, self%v(:, i))
       self%w = self%w - self%h(i, j) * self%v(:, i)
    end do
    self%h(j+1, j) = norm2(self%w)
    breakdown = self%h(j+1, j) <= epsilon(1.0_dp) * wnorm
    if (.not. breakdown) self%v(:, j+1) = self%w / self%h(j+1, j)

    do i = 1, j - 1
       temp = self%cs(i) * self%h(i, j) + self%sn(i) * self%h(i+1, j)
       self%h(i+1, j) = -self%sn(i) * self%h(i, j) + self%cs(i) * self%h(i+1, j)
       self%h(i, j) = temp
    end do
    temp = hypot(self%h(j, j), self%h(j+1, j))
    self%cs(j) = self%h(j, j) / temp
    self%sn(j) = self%h(j+1, j) / temp
    self%h(j, j) = temp
    self%g(j+1) = -self%sn(j) * self%g(j)
    self%g(j) = self%cs(j) * self%g(j)

    ! GMRES residual V_{j+1} Q_j^T g(j+1) e_{j+1} and iterate x + V_j R^{-1} g.
    z(1:j+1) = 0.0_dp
    z(j+1) = self%g(j+1)
    do i = j, 1, -1
       temp = self%cs(i) * z(i) - self%sn(i) * z(i+1)
       z(i+1) = self%sn(i) * z(i) + self%cs(i) * z(i+1)
       z(i) = temp
    end do
    self%w = matmul(self%v(:, 1:j+1), z(1:j+1))
    temp = sum(abs(self%w))
    do i = j, 1, -1
       z(i) = (self%g(i) - dot_product(self%h(i, i+1:j), z(i+1:j))) / self%h(i, i)
    end do
    self%w = self%x + matmul(self%v(:, 1:j), z(1:j))
    small = temp <= self%eps * sum(abs(self%w))

    if (small .or. breakdown .or. j == self%m .or. &
        self%napply >= self%maxapply .or. .not. ieee_is_finite(temp)) then
       self%x = self%w
       self%vin = self%x
       self%arnoldi = .false.
       self%done = self%napply >= self%maxapply
    else
       self%vin = self%v(:, j+1)
    end if
  end subroutine arnoldi_step

end module fixed_point_gmres_mod
