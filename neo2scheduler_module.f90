!> Module for the class Neo2Scheduler
module neo2scheduler_module

  use scheduler_module
  use wuGenericNeo2Workunit_module

  implicit none

  !> Neo2Scheduler, which is responsible for NEO-2 parallelization
  type, extends(scheduler) :: neo2scheduler

  contains
    procedure :: initMaster => initMaster_neo2scheduler

    procedure :: createWorkunit => createWorkunit_neo2scheduler
    procedure :: allocateSpecific => allocateSpecific_neo2scheduler
    procedure :: allocateSpecificMergeWU  => allocateSpecificMergeWU_neo2scheduler

    procedure :: deinit => deinit_neo2scheduler
  end type neo2scheduler

contains

  !> Init routine of neo2scheduler
  subroutine initMaster_neo2scheduler(this)
    class(neo2scheduler) :: this
    integer :: client_count
    integer :: fperiod_count
    integer :: current_client
    integer :: proptag_start_client, proptag_end_client
    integer :: current_workunit_client

    ! Calculate client count and fieldperiods count
    client_count = mpro%getNumProcs() - 1
    fperiod_count = parallel_storage%fieldline%ch_las%tag - parallel_storage%fieldline%ch_fir%tag + 1

    ! Part the number of fieldperiods to the number of clients
    call this%partNearlyFair(fperiod_count)

    ! Store some global variables
    parallel_storage%fieldperiod => parallel_storage%fieldline%ch_fir

    ! Iterate through the number of clients
    do current_client = 1, client_count

      ! Iterate over number of workunits, which the client will have
      do current_workunit_client = 1, this%workunits_per_client(current_client)

        ! Fill workunit with data
        proptag_start_client = parallel_storage%fieldperiod%ch_fir%tag
        proptag_end_client   = parallel_storage%fieldperiod%ch_las%tag

        !write (*,*) "Client", current_client, proptag_start_client, proptag_end_client

        ! Create the workunit
        call this%createWorkunit(current_client, proptag_start_client, proptag_end_client)

        ! Go to the next fieldperiod
        if (associated(parallel_storage%fieldperiod%next)) then
          parallel_storage%fieldperiod => parallel_storage%fieldperiod%next
        else
          write (*,*) "This was the last fieldperiod!"
        end if
      end do
    end do

  end subroutine initMaster_neo2scheduler

  !> Inherited from genericScheduler
  function allocateSpecific_neo2scheduler(this, wuType) result(res)
    class(neo2scheduler) :: this
    character(len=maxStrLen) :: wuType
    class(workunit), pointer :: res

    nullify(res)
    select case (wuType)
      case ("wuSolvePropagator")
        allocate(wuSolvePropagator :: res)
      case ("wuExternalJoin")
        allocate(wuExternalJoin :: res)
    end select

  end function allocateSpecific_neo2scheduler

  function allocateSpecificMergeWU_neo2scheduler(this) result(res)
    class(neo2scheduler) :: this
    class(wuMergeWorkunit), pointer :: res

    nullify(res)
    allocate(wuExternalJoin :: res)

  end function allocateSpecificMergeWU_neo2scheduler

  !> This function creates a wuSolvePropagator workunit and adds it to the unprocessed workunits list
  subroutine createWorkunit_neo2scheduler(this, client, proptag_start, proptag_end)
    class(neo2scheduler) :: this
    integer :: client
    integer :: proptag_start
    integer :: proptag_end

    integer :: ln
    integer :: rn
    type(wuSolvePropagator), pointer, save :: wu => null()

    ! Initialize the neighbors
    ln = -1
    rn = -1
    if (associated(wu)) then
      ln = wu%uid
      wu%rightNeighbor = mpro%storage%nextUID
    end if

    ! Allocate and fill the workunit
    allocate(wu)
    call wu%init()
    wu%client               = client
    wu%fracIndex            = proptag_start
    wu%proptag_start_client = proptag_start
    wu%proptag_end_client   = proptag_end

    ! Add the workunit to the unprocessed list
    call this%addWorkunit(wu)
    call wu%setNeighbors(ln, rn)

  end subroutine createWorkunit_neo2scheduler

  subroutine deinit_neo2scheduler(this)
    class(neo2scheduler) :: this

    !write (*,*) "Client", mpro%getRank(), " Solver: ", parallel_Storage%timeSolver / parallel_Storage%countSolver, &
    !                                                   parallel_Storage%countSolver, &
    !                                      " Joiner: ", parallel_Storage%timeJoiner / parallel_Storage%countJoiner, &
    !                                                   parallel_Storage%countJoiner

    call this%scheduler%deinit()
  end subroutine deinit_neo2scheduler

end module neo2scheduler_module
