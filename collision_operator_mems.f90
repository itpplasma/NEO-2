module collop
  use rkstep_mod, only : lag,leg,asource,anumm,denmm,ailmm,weightlag,anumm_ms, denmm_ms

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

    subroutine collop_set_species(ispec)
      integer :: ispec

      write (*,*) "Setting species to ", ispec
      anumm => anumm_ms(:,:,ispec)      
      denmm => denmm_ms(:,:,ispec)
    end subroutine collop_set_species
    
    subroutine collop_load()

      
      
      call collop_set_species(0)
    end subroutine collop_load

    subroutine collop_unload()

    end subroutine collop_unload
  
    subroutine collop_deconstruct()

    end subroutine collop_deconstruct
    
end module collop
