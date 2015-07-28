module collop

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
      use rkstep_mod, only : lag,leg,asource,anumm,denmm,ailmm,weightlag
      
      if(allocated(anumm)) deallocate(anumm)
      allocate(anumm(0:lag,0:lag,0:num_spec-1))
      !allocate(anumm(0:lag,0:lag))
      if(allocated(denmm)) deallocate(denmm)
      allocate(denmm(0:lag,0:lag,0:num_spec-1))
      !allocate(denmm(0:lag,0:lag))
      if(allocated(asource)) deallocate(asource)
      allocate(asource(0:lag,3))
      if(allocated(ailmm)) deallocate(ailmm)
      allocate(ailmm(0:lag,0:lag,0:leg))
      if(allocated(weightlag)) deallocate(weightlag)
      allocate(weightlag(3,0:lag))
    end subroutine collop_construct

    subroutine collop_load()
      use rkstep_mod, only : lag,leg,asource,anumm,denmm,ailmm,weightlag
      
    end subroutine collop_load

    subroutine collop_unload()

    end subroutine collop_unload
  
    subroutine collop_deconstruct()

    end subroutine collop_deconstruct
    
end module collop
