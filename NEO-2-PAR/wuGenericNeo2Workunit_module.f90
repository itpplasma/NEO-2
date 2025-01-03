!> Module for Neo2Workunit
module wuGenericNeo2Workunit_module

  use scheduler_module
  use mpiprovider_module
  use parallelStorage_module

  !> Class wuGenericNeo2Workunit
  type, extends(wuMergeWorkunit) :: wuGenericNeo2Workunit
    type(propagator), pointer :: prop_res => NULL()
    contains
      procedure :: packPropagator => packPropagator_wuGenericNeo2Workunit         !< For packing a propagator
      procedure :: packBinarySplit => packBinarySplit_wuGenericNeo2Workunit       !< Helper function for packing the binarySplit
      procedure :: packQEPropagator => packQEPropagator_wuGenericNeo2Workunit     !< Helper function for packing the QEPropagator

      procedure :: unpackPropagator => unpackPropagator_wuGenericNeo2Workunit     !< For unpacking a propagator
      procedure :: unpackBinarySplit => unpackBinarySplit_wuGenericNeo2Workunit   !< Helper function for unpacking the binarySplit
      procedure :: unpackQEPropagator => unpackQEPropagator_wuGenericNeo2Workunit !< Helper function for unpacking the QEPropagator
  end type wuGenericNeo2Workunit

  !> Workunit for calling the routine propsolve_all
  type, extends(wuGenericNeo2Workunit) :: wuSolvePropagator
    integer :: proptag_start_client   !< Defines the tag of the first propagator in the fieldperiod
    integer :: proptag_end_client     !< Defines the tag of the last propagator in the fieldperiod
  contains
    procedure :: init   => init_wuSolvePropagator     !< Initialize the workunit
    procedure :: process => process_wuSolvePropagator !< Call propsolve_all with the start and end tag
    procedure :: pack => pack_wuSolvePropagator       !< Pack the result
    procedure :: unpack => unpack_wuSolvePropagator   !< Unpack a received result

    procedure :: print => print_wuSolvePropagator     !< Print debug information
    
  end type wuSolvePropagator
  
  !> Workunit for joining two fieldperiods
  type, extends(wuGenericNeo2Workunit) :: wuExternalJoin

    integer :: iend
    integer :: fieldperiod1_uid     !< First period to join
    integer :: fieldperiod2_uid     !< Second period to join

  contains

    procedure :: init   => init_wuExternalJoin      !< Init the workunit
    procedure :: pack   => pack_wuExternalJoin      !< Pack the result
    procedure :: unpack => unpack_wuExternalJoin    !< Unpack a received result
    procedure :: process => process_wuExternalJoin  !< Call the joining process
    procedure :: print   => print_wuExternalJoin    !< Print debug information

    procedure :: setMergeInfo => setMergeInfo_wuExternalJoin

  end type wuExternalJoin

  contains

  !> Init the workunit to set the type
  subroutine init_wuSolvePropagator(this)
    class(wuSolvePropagator) :: this

    call this%wuMergeWorkunit%init()
    this%type = "wuSolvePropagator"
    this%balance = .true.
  end subroutine init_wuSolvePropagator

  !> Solve propagators on the client
  subroutine process_wuSolvePropagator(this)
    class(wuSolvePropagator) :: this
    double precision :: stime

    write (*,*) "Client", mpro%getRank(), " calls propagator_solver() with", this%proptag_start_client, this%proptag_end_client

    stime = MPI_WTime()

    ! Call the Neo-2 algorithm for solving propagators from propagator.f90
    call propagator_solver(this%proptag_start_client, this%proptag_end_client, &
                           globalstorage%bin_split_mode, globalstorage%eta_ori, parallelMode = .true.)

    ! Do some performance recordings
    globalstorage%countSolver = globalstorage%countSolver + 1
    globalstorage%timeSolver = globalstorage%timeSolver + (MPI_WTime() - stime)

    ! Remember the pointer to the result
    if (associated(prop_c)) then

       if (associated(prop_c%prev)) then
       
          ! Relinking pointers for storing the propagator without copying its content
          this%prop_res => prop_c%prev
          this%prop_res%next => null()
          prop_c%prev => null()
          prop_r => prop_c

       else
          if (prop_reconstruct == 2) then
             write (*,*) "prop_c%prev is not associated because of prop_reconstruct==2"
          else
             write (*,*) "ERROR: prop_c%prev is not associated"
             stop 
          end if
       end if
    else
       write (*,*) "An error occurred in process_wuSolvePropagator, one of the propagators is null", &
            associated(this%prop_res), associated(prop_c)
      stop
    end if

    ! Free the source propagators
    call destruct_all_propagators()
    this%isProcessed = .true.
  end subroutine process_wuSolvePropagator

  !> Call helper functions to pack the result of the workunit
  subroutine pack_wuSolvePropagator(this)
    class(wuSolvePropagator) :: this

    call this%genericWorkunit%pack()
  
    associate (b => mpro%packbuffer)
      call b%add(this%proptag_start_client)
      call b%add(this%proptag_end_client)
    end associate
    call this%packPropagator()
  end subroutine pack_wuSolvePropagator

  !> Call helper functions to unpack the results
  subroutine unpack_wuSolvePropagator(this)
    class(wuSolvePropagator) :: this

    call this%genericWorkunit%unpack()

    associate (b => mpro%packbuffer)
      call b%get(this%proptag_start_client)
      call b%get(this%proptag_end_client)
    end associate
    call this%unpackPropagator()
  end subroutine unpack_wuSolvePropagator

  !> Print debug information
  subroutine print_wuSolvePropagator(this)
    class(wuSolvePropagator) :: this

    write (*,*) trim(this%type), this%uid, this%leftNeighbor, this%rightNeighbor, &
      this%proptag_start_client, this%proptag_end_client
  end subroutine print_wuSolvePropagator

  !> Init the workunit for external joining
  subroutine init_wuExternalJoin(this)
    class(wuExternalJoin) :: this

    call this%wuMergeWorkunit%init()
    this%type = "wuExternalJoin"
  end subroutine init_wuExternalJoin

  !> Call the external join routines
  subroutine process_wuExternalJoin(this)
    class(wuExternalJoin) :: this
    type(propagator), pointer :: prop1 => null()
    type(propagator), pointer :: prop2 => null()

    double precision :: stime
    class(workunit), pointer :: wu => null()


    ! Search the first propagator to join
    wu => mpro%storage%waitingWorkunits%get(this%fieldperiod1_uid)
    select type (wu)
      class is (wuGenericNeo2Workunit)
        prop1 => wu%prop_res
      class default
        write (*,*) "An error occurred in wuExternal_join, undefined work unit type"
        stop
    end select

    ! Search the second propagator to join
    wu => mpro%storage%waitingWorkunits%get(this%fieldperiod2_uid)
    select type (wu)
      class is (wuGenericNeo2Workunit)
        prop2 => wu%prop_res
      class default
        write (*,*) "An error occurred in wuExternal_join, undefined work unit type"
        stop
    end select

    write (*,*) "  Calling external join for fieldperiods: ", this%fieldperiod1_uid, this%fieldperiod2_uid

    ! Check if they have been found
    if (associated(prop1) .and. associated(prop2)) then

      ! Set global variables to define the two propagator to be joined
      prop_c_old => prop1
      prop_c_new => prop2

      stime = MPI_WTime()

      ! Call Neo-2 algorithm to join in propagator.f90
      call external_joining()

      ! Do some performance recordings
      globalstorage%timeJoiner = globalstorage%timeJoiner + (MPI_WTime() - stime)
      globalstorage%countJoiner = globalstorage%countJoiner + 1

      ! Relinking pointers for storing the propagator without copying its content
      this%prop_res => prop_c_old

      if (this%iend == 1) then

        ! Copy one propagator
        call assign_propagator_content(prop_c_new, prop_c_old)

        ! Call final join from propagator.f90
        call final_joining()
      end if

      this%isProcessed = .true.
    else
      write (*,*) "An error occurred in process_wuExternalJoin, some of the propagators are not associated", &
                  associated(prop1), associated(prop2)
      stop
    end if

  end subroutine process_wuExternalJoin

  !> Call helper functions to pack the result
  subroutine pack_wuExternalJoin(this)
    class(wuExternalJoin) :: this

    call this%genericWorkunit%pack()
    associate (b => mpro%packbuffer)
      call b%add(this%iend)
      call b%add(this%fieldperiod1_uid)
      call b%add(this%fieldperiod2_uid)
    end associate
    call this%packPropagator()

  end subroutine pack_wuExternalJoin

  !> Call helper functions to unpack the result
  subroutine unpack_wuExternalJoin(this)
    class(wuExternalJoin) :: this

    call this%genericWorkunit%unpack()
    associate (b => mpro%packbuffer)
      call b%get(this%iend)
      call b%get(this%fieldperiod1_uid)
      call b%get(this%fieldperiod2_uid)
    end associate
    call this%unpackPropagator()

  end subroutine unpack_wuExternalJoin

  !> Print debug information
  subroutine print_wuExternalJoin(this)
    class(wuExternalJoin) :: this

    write (*,*) trim(this%type), this%uid, this%fieldperiod1_uid, this%fieldperiod2_uid
  end subroutine print_wuExternalJoin

  subroutine setMergeInfo_wuExternalJoin(this, left_uid, right_uid)
      class(wuExternalJoin) :: this
      integer :: left_uid, right_uid

      this%fieldperiod1_uid = left_uid
      this%fieldperiod2_uid = right_uid
      this%iend = 0
      if (this%leftNeighbor == this%rightNeighbor) then
          this%iend = 1
          write (*,*) "Scheduler detected last joining process, setting iend = 1 and calling final join"
      end if


  end subroutine setMergeInfo_wuExternalJoin

  !> Helper function to pack a propagator
  subroutine packPropagator_wuGenericNeo2Workunit(this)
    class(wuGenericNeo2Workunit) :: this

    if (this%isProcessed) then

      associate (b => mpro%packBuffer)
        associate (prop => this%prop_res)

          call b%add(prop%tag)
          call b%add(prop%nr_joined)
          call b%add(prop%bin_split_mode)

          call this%packBinarySplit(prop%eta_bs_l)
          call this%packBinarySplit(prop%eta_bs_r)

          call b%add(prop%fieldpropagator_tag_s)
          call b%add(prop%fieldpropagator_tag_e)
          call b%add(prop%fieldperiod_tag_s)
          call b%add(prop%fieldperiod_tag_e)
          call b%add(prop%y)
          call b%add(prop%phi_l)
          call b%add(prop%phi_r)

          call this%packQEPropagator(prop%p)
        end associate
      end associate

    end if
  end subroutine packPropagator_wuGenericNeo2Workunit

  !> Helper function to unpack a propagator
  subroutine unpackPropagator_wuGenericNeo2Workunit(this)
    class(wuGenericNeo2Workunit) :: this

    if (this%isProcessed)  then

      allocate(this%prop_res)

      associate (b => mpro%packBuffer)
        associate (prop => this%prop_res)
          write (*,*) "Unpacking Propagator"
          call b%get(prop%tag)
          call b%get(prop%nr_joined)
          call b%get(prop%bin_split_mode)

          call this%unpackBinarySplit(prop%eta_bs_l)
          call this%unpackBinarySplit(prop%eta_bs_r)

          call b%get(prop%fieldpropagator_tag_s)
          call b%get(prop%fieldpropagator_tag_e)
          call b%get(prop%fieldperiod_tag_s)
          call b%get(prop%fieldperiod_tag_e)

          call b%get(prop%y)
          call b%get(prop%phi_l)
          call b%get(prop%phi_r)

          call this%unpackQEPropagator(prop%p)

        end associate
      end associate

    end if
  end subroutine unpackPropagator_wuGenericNeo2Workunit

  !> Pack QEPropagator
  subroutine packQEPropagator_wuGenericNeo2Workunit(this, p)
    class(wuGenericNeo2Workunit) :: this
    type(prop_qe) :: p

    associate (b => mpro%packBuffer)
      call b%add(p%npart)
      call b%add(p%npass_l)
      call b%add(p%npass_r)
      call b%add(p%nvelocity)

      call b%add(p%amat_p_p)
      call b%add(p%amat_m_m)
      call b%add(p%amat_p_m)
      call b%add(p%amat_m_p)

      call b%add(p%source_p)
      call b%add(p%source_m)

      call b%add(p%flux_p)
      call b%add(p%flux_m)

      call b%add(p%qflux)
      call b%add(p%cmat)

      call b%add(p%eta_l)
      call b%add(p%eta_r)

      call b%add(p%eta_boundary_l)
      call b%add(p%eta_boundary_r)
      call b%add(p%w)

    end associate
  end subroutine packQEPropagator_wuGenericNeo2Workunit

  !> Unpack QEPropagator
  subroutine unpackQEPropagator_wuGenericNeo2Workunit(this, p)
    class(wuGenericNeo2Workunit) :: this
    type(prop_qe) :: p

    associate (b => mpro%packBuffer)
      call b%get(p%npart)
      call b%get(p%npass_l)
      call b%get(p%npass_r)
      call b%get(p%nvelocity)

      call b%get(p%amat_p_p)
      call b%get(p%amat_m_m)
      call b%get(p%amat_p_m)
      call b%get(p%amat_m_p)

      call b%get(p%source_p)
      call b%get(p%source_m)

      call b%get(p%flux_p)
      call b%get(p%flux_m)

      call b%get(p%qflux)
      call b%get(p%cmat)

      call b%get(p%eta_l)
      call b%get(p%eta_r)

      call b%get(p%eta_boundary_l)
      call b%get(p%eta_boundary_r)

      call b%get(p%w)

    end associate
  end subroutine unpackQEPropagator_wuGenericNeo2Workunit

  !> Pack binarySplit
  subroutine packBinarySplit_wuGenericNeo2Workunit(this, bs)
    class(wuGenericNeo2Workunit) :: this
    type(binarysplit) :: bs
    integer :: k

    associate (b => mpro%packBuffer)

        call b%add(bs%n_ori)
        call b%add(bs%n_split)

        do k = 0, bs%n_ori
           call b%add(bs%x_ori_bin_sparse(k)%idxvec)
           call b%add(bs%x_ori_bin_sparse(k)%values)
           call b%add(bs%x_ori_bin_sparse(k)%len)
           call b%add(bs%x_ori_bin_sparse(k)%len_sparse)
        end do
        
        call b%add(bs%x_ori_poi)
        call b%add(bs%x_poi)
        call b%add(bs%x_split)
        call b%add(bs%x_pos)
        call b%add(bs%x)
        call b%add(bs%y)
        call b%add(bs%int)
        call b%add(bs%err)

    end associate
  end subroutine packBinarySplit_wuGenericNeo2Workunit

  !> Unpack binarySplit
  subroutine unpackBinarySplit_wuGenericNeo2Workunit(this, bs)
    class(wuGenericNeo2Workunit) :: this
    type(binarysplit) :: bs

    associate (b => mpro%packBuffer)
      call b%get(bs%n_ori)
      call b%get(bs%n_split)
      !call b%get(bs%x_ori_bin)

      if (allocated(bs%x_ori_bin_sparse)) deallocate(bs%x_ori_bin_sparse)
      allocate(bs%x_ori_bin_sparse(0:bs%n_ori))
      do k = 0, bs%n_ori
         call b%get(bs%x_ori_bin_sparse(k)%idxvec)
         call b%get(bs%x_ori_bin_sparse(k)%values)
         call b%get(bs%x_ori_bin_sparse(k)%len)
         call b%get(bs%x_ori_bin_sparse(k)%len_sparse)
      end do

      call b%get(bs%x_ori_poi)
      call b%get(bs%x_poi)
      call b%get(bs%x_split)
      call b%get(bs%x_pos)
      call b%get(bs%x)
      call b%get(bs%y)
      call b%get(bs%int)
      call b%get(bs%err)

    end associate
  end subroutine unpackBinarySplit_wuGenericNeo2Workunit

end module wuGenericNeo2Workunit_module
