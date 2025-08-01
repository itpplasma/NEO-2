!> Module for controlling spline implementation behavior during testing
module spline_test_control_mod
  use nrtype, only: DP
  implicit none
  
  !> Flag to disable fast path for testing purposes
  logical :: disable_fast_path = .false.
  
end module spline_test_control_mod