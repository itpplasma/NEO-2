module collop
  use rkstep_mod!, only : lag,leg,asource,anumm,denmm,ailmm,weightlag,anumm_a, denmm_a, anumm_aa, denmm_aa
  use collop_compute, only : compute_collop, compute_collop_inf, compute_source, m_d, m_C, m_ele, write_collop
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
      num_spec = 1 !mpro%getNumProcs()
      !write (*,*) "Number of species ", num_spec

    end subroutine collop_construct

    subroutine collop_set_species(ispec, opt_talk)
      integer :: ispec
      logical, optional :: opt_talk
      logical :: talk

      talk = .true.
      if (present(opt_talk)) talk = opt_talk
      
      !write (*,*) "Setting species to ", ispec

      !**********************************************************
      ! Switch collision operator matrices
      !**********************************************************
      anumm(0:lag, 0:lag) => anumm_a(:,:,ispec)      
      denmm(0:lag, 0:lag) => denmm_a(:,:,ispec)
      ailmm(0:lag, 0:lag, 0:leg) => ailmm_aa(:,:,:,mpro%getRank(),ispec)

      !**********************************************************
      ! Switch collisionality parameter
      !**********************************************************
      ! Not yet implemented
      
    end subroutine collop_set_species
    
    subroutine collop_load()
      real(kind=dp), dimension(:), allocatable :: asource_temp
      real(kind=dp), dimension(:,:), allocatable :: anumm_inf
      
      !**********************************************************
      ! Allocation of matrices
      !**********************************************************
      if(allocated(anumm_aa)) deallocate(anumm_aa)
      allocate(anumm_aa(0:lag,0:lag,0:num_spec-1,0:num_spec-1))

      if(allocated(anumm_a)) deallocate(anumm_a)
      allocate(anumm_a(0:lag,0:lag,0:num_spec-1))
      
      if(allocated(denmm_aa)) deallocate(denmm_aa)
      allocate(denmm_aa(0:lag,0:lag,0:num_spec-1,0:num_spec-1))

      if(allocated(denmm_a)) deallocate(denmm_aa)
      allocate(denmm_a(0:lag,0:lag,0:num_spec-1))
      
      if(allocated(asource)) deallocate(asource)
      allocate(asource(0:lag,3))
      
      if(allocated(ailmm_aa)) deallocate(ailmm_aa)
      allocate(ailmm_aa(0:lag,0:lag,0:leg,0:num_spec-1,0:num_spec-1))
      
      if(allocated(weightlag)) deallocate(weightlag)
      allocate(weightlag(3,0:lag))

      if (allocated(anumm_inf)) deallocate(anumm_inf)
      allocate(anumm_inf(0:lag, 0:lag))

      !**********************************************************
      ! Compute sources
      !**********************************************************
      call compute_source(asource, weightlag)

      !**********************************************************
      ! Compute collision operator
      !**********************************************************
      call compute_collop_inf('e', 'e', m_ele, m_ele, 1d0, 1d0, anumm_aa(:,:,0,0), anumm_inf, &
           denmm_aa(:,:,0,0), ailmm_aa(:,:,:,0,0))

      !**********************************************************
      ! Sum up matrices
      !**********************************************************
      anumm_a(:,:,0) = anumm_aa(:,:,0,0) + Z_eff * anumm_inf(:,:)
      denmm_a(:,:,0) = denmm_aa(:,:,0,0)
     
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
      
      !**********************************************************
      ! Set pointers to main species
      !**********************************************************
      call collop_set_species(0)

      !**********************************************************
      ! Write to screen
      !**********************************************************
      !write (*,*) asource
      !write (*,*) anumm
      !write (*,*) denmm
      !write (*,*) ailmm
      !write (*,*) weightlag
      
    end subroutine collop_load

    subroutine collop_unload()
      deallocate(anumm_aa)
      deallocate(denmm_aa)
      deallocate(anumm_a)
      deallocate(denmm_a)
      deallocate(asource)
      deallocate(weightlag)
    end subroutine collop_unload
  
    subroutine collop_deconstruct()

    end subroutine collop_deconstruct
    
end module collop
