module parallelStorage_module

     use propagator_mod
     use magnetics_mod

     implicit none

     type :: parallelStorage
          ! Special for workunit creation:
          integer :: bin_split_mode
          REAL(kind=KIND(1.0d0)), allocatable :: eta_ori(:)
          TYPE(fieldperiod_struct), POINTER :: fieldperiod
          TYPE(fieldline_struct), POINTER :: fieldline

          double precision :: timeSolver = 0
          double precision :: timeJoiner = 0
          integer :: countSolver = 0
          integer :: countJoiner = 0
 
     end type

     type(parallelStorage) :: globalstorage

end module parallelStorage_module
