module collop
  use rkstep_mod, only : lag,leg,asource,anumm,denmm,ailmm,weightlag,anumm_ms, denmm_ms
  use collop_compute, only : compute_collop, m_d, m_C
  use mpiprovider_module

  implicit none
  
  !**********************************************************
  ! From old module, mainly for compatibility
  !**********************************************************
  integer, parameter, private :: dp = kind(1.0d0)
  integer, parameter, private :: dummy_read = 20
  real(kind=dp), public       :: z_eff = 1.0_dp
  logical, public             :: collop_talk      =  .true. 
  logical, public             :: collop_talk_much =  .true.
  character(len=100), public  :: collop_path

  !**********************************************************
  ! Number of species
  !**********************************************************
  integer :: num_spec

  contains
    
    subroutine collop_construct()     
      use nrtype, only : pi
      integer          :: openstatus, readstatus, k, m
      character(len=1) :: dummy

      !**********************************************************
      ! Test configuration two species
      !**********************************************************
      num_spec = mpro%getNumProcs()
      !write (*,*) "Number of species ", num_spec

      if(allocated(anumm_ms)) deallocate(anumm_ms)
      allocate(anumm_ms(0:lag,0:lag,0:num_spec-1))

      if(allocated(denmm_ms)) deallocate(denmm_ms)
      allocate(denmm_ms(0:lag,0:lag,0:num_spec-1))

      if(allocated(asource)) deallocate(asource)
      allocate(asource(0:lag,3))
      
      if(allocated(ailmm)) deallocate(ailmm)
      allocate(ailmm(0:lag,0:lag,0:leg))
      
      if(allocated(weightlag)) deallocate(weightlag)
      allocate(weightlag(3,0:lag))

    end subroutine collop_construct

    subroutine collop_set_species(ispec, opt_talk)
      integer :: ispec
      logical, optional :: opt_talk
      logical :: talk

      talk = .true.
      if (present(opt_talk)) talk = opt_talk
      
      write (*,*) "Setting species to ", ispec

      !**********************************************************
      ! Switch collision operator matrices
      !**********************************************************
      anumm(0:lag, 0:lag) => anumm_ms(:,:,ispec)      
      denmm(0:lag, 0:lag) => denmm_ms(:,:,ispec)

      !**********************************************************
      ! Switch collisionality parameter
      !**********************************************************
      ! Not yet implemented
      
    end subroutine collop_set_species
    
    subroutine collop_load()
      real(kind=dp), dimension(:), allocatable :: asource_temp

      !**********************************************************
      ! Test configuration dd and dC
      !**********************************************************   
      call collop_set_species(0, .false.)
      call compute_collop('d', 'd' , m_d, m_d, 1d0, 1d0, asource, weightlag, anumm, denmm, ailmm)

      if (mpro%getNumProcs() .ge. 2) then
         call collop_set_species(1, .false.)
         call compute_collop('d', 'C' , m_d, m_C, 1d0, 1d0, asource, weightlag, anumm, denmm, ailmm)
      end if

      if (mpro%getNumProcs() .ge. 3) then
         write (*,*) "More than two species not yet supported"
         stop
      end if

      !**********************************************************
      ! Swap sources for NEO-2 convention
      !**********************************************************
      allocate(asource_temp(0:lag))
      asource_temp = asource(:, 2)
      asource(:,2) = asource(:,3)
      asource(:,3) = asource_temp

      asource_temp = weightlag(2,:)
      weightlag(2,:) = weightlag(3,:)
      weightlag(3,:) = asource_temp
      deallocate(asource_temp)

      !write (*,*) weightlag(1,:)
      !write (*,*) weightlag(2,:)
      !write (*,*) weightlag(3,:)
      !stop
      
      !**********************************************************
      ! Set to main species
      !**********************************************************
      call collop_set_species(0)
    end subroutine collop_load

    subroutine collop_unload()

    end subroutine collop_unload
  
    subroutine collop_deconstruct()

    end subroutine collop_deconstruct
    
end module collop
