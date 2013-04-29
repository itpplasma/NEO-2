!> Module for the class Neo2Scheduler
module neo2scheduler_module

  use scheduler_module
  use wuGenericNeo2Workunit_module

  implicit none

  !> Neo2Scheduler, which is responsible for
  type, extends(scheduler) :: neo2scheduler

  contains
    procedure :: initMaster => initMaster_neo2scheduler

    procedure :: createWorkunit => createWorkunit_neo2scheduler
    procedure :: allocateSpecific => allocateSpecific_neo2scheduler
    procedure :: checkIfClientDone => checkIfClientDone_neo2scheduler
    procedure :: loadBalancing => loadBalancing_neo2scheduler

    procedure :: deinit => deinit_neo2scheduler
  end type neo2scheduler

contains

  !> Define behaviour of loadbalancing
  subroutine loadBalancing_neo2Scheduler(this)
    class(neo2scheduler) :: this
    integer :: k

    class(workunit), pointer :: selectWU => null()

    ! If balancing is activated in config file
    if (this%balance) then

      ! Iterate through all clients
      do k = 1, mpro%getNumProcs()-1

        ! Check if the client has no workunit at the moment
        if (this%clientStats(k)%isReady) then

          ! Check if the client has nothing more to do
          if (this%clientStats(k)%isDone) then
            !write (*,*) "Client", k, "has nothing more to do."
            call mpro%storage%waitingWorkunits%rewind()

            ! Search unprocessed workunit
            do while (associated(mpro%storage%waitingWorkunits%currentElement))
              ! Tell the fortran compiler that the list contains workunits :-)
              selectWU => mpro%storage%waitingWorkunits%getCurrent()
              select type (q => selectWU)
                class is (wuSolvePropagator)

                ! Search next free workunit
                if ((q%client /= k) .and. (.not. this%clientStats(q%client)%isReady)) then
                  write (*,*) "------- Balancing", q%uid, "from", q%client, "to", k, " -------"

                  ! Set client of the workunit to the client, which is done
                  call q%setClient(k)
                  mpro%balanceCount = mpro%balanceCount + 1
                  exit
                end if
              end select

              call mpro%storage%waitingWorkunits%gotoNext()
            end do
          else
            ! If a client is ready, but not done yet, set it to done
            ! If it will be done the next iteration of scheduler, it will get a balanced workunit
            this%clientStats(k)%isDone = .true.
          end if
        end if
      end do
    end if

  end subroutine loadBalancing_neo2Scheduler

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

  subroutine deinit_neo2scheduler(this)
    class(neo2scheduler) :: this

    !write (*,*) "Client", mpro%getRank(), " Solver: ", parallel_Storage%timeSolver / parallel_Storage%countSolver, &
    !                                                   parallel_Storage%countSolver, &
    !                                      " Joiner: ", parallel_Storage%timeJoiner / parallel_Storage%countJoiner, &
    !                                                   parallel_Storage%countJoiner

    call this%scheduler%deinit()
  end subroutine deinit_neo2scheduler

  !> Check if workunits have to be created
  subroutine checkIfClientDone_neo2scheduler(this)
    class(neo2scheduler) :: this
    integer :: k
    type(wuDataRequester), pointer :: dr
    class(wuExternalJoin), pointer :: wu
    class(wuMergeWorkunit), pointer :: p
    integer :: temp

    class(workunit), pointer :: selectWU => null()
    class(workunit), pointer :: selectWU2 => null()

    do k = 1, mpro%getRank() -1
      this%clientStats(k)%localMerge = .false.
    end do

    ! In the first iteration, local merges will be done, in the second iteration merges over different clients will be done
    ! So, local merge has always higher priority
    do k = 1, 2

      call mpro%storage%processedWorkunits%rewind()
      !if ((mpro%waitingWorkunits%getCount() == 0 .and. mpro%pendingWorkunits%getCount() == 0) .or. this%readyToMerge) then

      do while (associated(mpro%storage%processedWorkunits%currentElement))

        selectWU => mpro%storage%processedWorkunits%getCurrent()
        select type (q => selectWU)
          class is (wuMergeWorkunit)

          !Do start with the second element
          call mpro%storage%processedWorkunits%gotoNext()

          !Check if element was not already involved in a merge-process
          if (.not. q%isMerged) then
            nullify(p)

            ! Get needed object to merge
            selectWU2 => mpro%storage%processedWorkunits%get(q%leftNeighbor, .false.)
            if (associated(selectWU2)) then

            select type (selectWU2)
              class is (wuMergeWorkunit)
              p => selectWU2
            end select

            end if

            ! If found (object was already processed and is in processedWorkunits list)
            if (associated(p)) then

              !write (*,*) "I would merge ", p%uid, q%uid, "on client", p%client

              !Check if clients have nothing to do
              if (this%clientStats(p%client)%isReady .and. this%clientStats(q%client)%isReady .and. &
                .not. this%clientStats(p%client)%isBlocked .and. .not. this%clientStats(q%client)%isBlocked) then

                ! Local merge has first priority !
                if ((k == 1) .and. (p%client == q%client)) then

                  this%clientStats(p%client)%localMerge = .true.

                  !wu => this%allocateMergeWU("wuExternalJoin")
                  allocate(wu)
                  call wu%init()
                  wu%fieldperiod1_uid = p%uid
                  wu%fieldperiod2_uid = q%uid
                  wu%client = p%client
                  call wu%setNeighbors(p%leftNeighbor, q%rightNeighbor)

                  wu%iend = 0
                  if (wu%leftNeighbor == wu%rightNeighbor) then
                    wu%iend = 1
                    write (*,*) "Scheduler detected last joining process, setting iend = 1 and calling final join"
                  end if

                  wu%fracIndex = 0
                  wu%resultUID = wu%uid

                  write (*,*) "Local merge on client", p%client, p%uid, q%uid

                  if (p%druid /= -1) then
                    call wu%neededWUs%add(p%druid)
                     ! write (*,*) "Setting dep"
                  end if

                  if (q%druid /= -1) then
                    call wu%neededWUs%add(q%druid)
                  end if

                  call mpro%storage%waitingWorkunits%add(wu)

                  call q%setMerged(.true.)
                  call this%repairNeighbors(p%uid, wu%leftNeighbor, wu%uid, wu%rightNeighbor)
                  call mpro%storage%processedWorkunits%rewind()


                  !call mpro%storage%waitingWorkunits%print()
                else
                  if ((k == 2)) then
                    if (q%oldClient == -1) then
                      nullify(dr)
                      allocate(dr)
                      call dr%init()
                      dr%client   = q%client
                      dr%dest     = p%client !Right to left
                      dr%whichUID = q%uid
                      dr%fracIndex = 0
                      call q%setOldClient(q%client)
                      !call dr%neededWUs%add(dr%whichUID)
                      call q%setClient(p%client)

                      !this%clientStats(dr%dest)%isReady = .false.
                      this%clientStats(dr%dest)%isBlocked = .true.
                      !write (*,*) "Blocking client", dr%dest

                      temp = dr%uid
                      call mpro%storage%waitingWorkunits%add(dr)
                      call p%setDrUID(temp)
                      call q%setDrUID(temp)

                      call myLog%logCreateDR(dr%uid, dr%client, dr%dest, dr%whichUID)
                       write (*,*) "Creating dataRequester", dr%uid,"for client", dr%client, "UID=", &
                            dr%whichUID, "Dest=", dr%dest

                    end if
                  end if
                end if

                 !exit
              else

              end if !Local merge ?

            end if ! Client are ready?

          else
            cycle
          end if ! workunit already merged

             ! Skip other workunits
          class default
          call mpro%storage%processedWorkunits%gotoNext()

        end select !Get only mergeWorkunits

      end do !Loop over processedWorkunits

    end do ! loop over k

  end subroutine checkIfClientDone_neo2scheduler

end module neo2scheduler_module
