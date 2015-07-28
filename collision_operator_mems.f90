module collop
  use rkstep_mod, only : lag,leg,asource,anumm,denmm,ailmm,weightlag,anumm_ms, denmm_ms
  use collop_compute, only : compute_collop, m_d, m_C

  implicit none
  
  !**********************************************************
  ! From old module, mainly for compatibility
  !**********************************************************
  integer, parameter, private :: dp = kind(1.0d0)
  integer, parameter, private :: dummy_read = 20
  real(kind=dp), public :: z_eff = 1.0_dp
  logical, public :: collop_talk      =  .true. 
  logical, public :: collop_talk_much =  .true.
  character(len=100), public :: collop_path

  !**********************************************************
  ! Number of species
  !**********************************************************
  integer :: num_spec
  
  contains
    
    subroutine collop_construct()     
      num_spec = 2
      
      if(allocated(anumm_ms)) deallocate(anumm_ms)
      allocate(anumm_ms(0:lag,0:lag,0:num_spec-1))
      !allocate(anumm(0:lag,0:lag))
      if(allocated(denmm_ms)) deallocate(denmm_ms)
      allocate(denmm_ms(0:lag,0:lag,0:num_spec-1))
      !allocate(denmm(0:lag,0:lag))
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
      anumm(0:lag, 0:lag) => anumm_ms(:,:,ispec)      
      denmm(0:lag, 0:lag) => denmm_ms(:,:,ispec)
    end subroutine collop_set_species
    
    subroutine collop_load()

      !**********************************************************
      ! Test configuration dd and dC
      !**********************************************************   
      call collop_set_species(0, .false.)
      call compute_collop('d', 'd' , m_d, m_d, 1d0, 1d0, asource, anumm, denmm, ailmm)
      
      call collop_set_species(1, .false.)
      call compute_collop('d', 'C' , m_d, m_C, 1d0, 1d0, asource, anumm, denmm, ailmm)

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
