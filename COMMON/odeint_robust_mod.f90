module odeint_robust_mod
  ! Adaptive ODE step that does not fail on components that start at zero.
  !
  ! The error scale of odeint_allroutines is |y| + |h*dy/dx| + 1e-30. A
  ! component that is exactly zero with zero derivative at the start of the
  ! interval (e.g. the running integrals y(11:14) at the start of a field
  ! line) makes the error test fail down to the roundoff level, so the step
  ! size underflows and odeint returns without integrating. Whether this
  ! happens depends on the last bits of the error estimate, i.e. on compiler
  ! and platform. odeint_robust then redoes the interval with these
  ! components shifted by one, which gives them an absolute tolerance.
  use odeint_allroutines_sub, only: odeint_allroutines, odeint_has_failed
  implicit none
  private

  public :: odeint_robust

  abstract interface
    subroutine rhs_interface(x, y, dydx)
      double precision, intent(in)  :: x
      double precision, intent(in)  :: y(:)
      double precision, intent(out) :: dydx(:)
    end subroutine rhs_interface
  end interface

  procedure(rhs_interface), pointer :: rhs_unshifted => null()
  double precision, allocatable :: y_offset(:)
  !$omp threadprivate(rhs_unshifted, y_offset)

contains

  subroutine odeint_robust(y, x1, x2, eps, derivs, ierr)
    double precision, intent(inout) :: y(:)
    double precision, intent(in)    :: x1, x2, eps
    procedure(rhs_interface)        :: derivs
    integer, intent(out)            :: ierr

    double precision :: y_start(size(y))

    ierr = 0
    y_start = y
    call odeint_allroutines(y, size(y), x1, x2, eps, derivs)
    if (.not. odeint_has_failed()) return

    y_offset = merge(1.0d0, 0.0d0, y_start .eq. 0.0d0)
    rhs_unshifted => derivs
    y = y_start + y_offset
    call odeint_allroutines(y, size(y), x1, x2, eps, rhs_shifted)
    y = y - y_offset
    rhs_unshifted => null()
    if (odeint_has_failed()) ierr = 1
  end subroutine odeint_robust

  subroutine rhs_shifted(x, y_shifted, dydx)
    double precision, intent(in)  :: x
    double precision, intent(in)  :: y_shifted(:)
    double precision, intent(out) :: dydx(:)

    call rhs_unshifted(x, y_shifted - y_offset, dydx)
  end subroutine rhs_shifted

end module odeint_robust_mod
