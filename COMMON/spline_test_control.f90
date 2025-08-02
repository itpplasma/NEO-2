!> Module for controlling spline implementation behavior during testing
module spline_test_control
  use nrtype, only: DP
  implicit none
  
  !> Flag to disable fast path for testing purposes
  logical :: disable_fast_path = .false.
  
  !> Flag to enable fast splines (configurable via neo2.in)
  logical :: use_fast_splines = .false.
  
end module spline_test_control