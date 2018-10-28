module parallelStorage_module

  use propagator_mod
  use magnetics_mod
  use nrtype, only : dp

  implicit none

  type :: parallelStorage
    ! Special for workunit creation:
    integer :: bin_split_mode
    real(kind=dp), allocatable :: eta_ori(:)
    TYPE(fieldperiod_struct), POINTER :: fieldperiod
    TYPE(fieldline_struct), POINTER :: fieldline

    real(kind=dp) :: timeSolver = 0
    real(kind=dp) :: timeJoiner = 0
    integer :: countSolver = 0
    integer :: countJoiner = 0

  end type

  type(parallelStorage) :: globalstorage

end module parallelStorage_module
