module collop
  use rkstep_mod
  !use nrtype, only, private : pi, dp
  use hdf5_tools
  use collop_compute, only : init_collop, &
       compute_source, compute_collop, gamma_ab, M_transform, M_transform_inv, &
       m_ele, m_d, m_C, compute_collop_inf
  use mpiprovider_module

  implicit none
  
  !**********************************************************
  ! From old module, mainly for compatibility
  !**********************************************************
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter, private :: dummy_read = 20
  real(kind=dp), public       :: z_eff = 1.0_dp
  logical, public             :: collop_talk      =  .true. 
  logical, public             :: collop_talk_much =  .true.
  character(len=100), public  :: collop_path
  integer                     :: collop_base_prj  = 0 ! 0...Laguerre (Default), 1...Polynomial
  integer                     :: collop_base_exp  = 0 
  real(kind=dp)               :: scalprod_alpha = 0d0
  real(kind=dp)               :: scalprod_beta  = 0d0

  !**********************************************************
  ! Number of species
  !**********************************************************
  integer :: num_spec

  real(kind=dp), dimension(:,:), allocatable :: anumm_inf

  contains
    
    subroutine collop_construct()     

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
      if (Z_eff .ne. 0) then
         ailmm(0:lag, 0:lag, 0:leg) => ailmm_aa(:,:,:,ispec,ispec)
      else
         ailmm(0:lag, 0:lag, 0:leg) => ailmm_aa(:,:,:,mpro%getRank(),ispec)
      end if
      !**********************************************************
      ! Switch collisionality parameter
      !**********************************************************
      ! Not yet implemented
      
    end subroutine collop_set_species
    
    subroutine collop_load()
      real(kind=dp), dimension(:), allocatable :: asource_temp
      real(kind=dp) :: alpha_temp, beta_temp
      integer       :: a,b

      if (Z_eff .ne. 0 ) then
         write (*,*) "Standard mode."
         num_spec = 1
      else
         write (*,*) "Test mode for two species."
         num_spec = 2
      end if
      
      !**********************************************************
      ! Allocation of matrices
      !**********************************************************
      if(allocated(anumm_aa)) deallocate(anumm_aa)
      allocate(anumm_aa(0:lag,0:lag,0:num_spec-1,0:num_spec-1))

      if(allocated(anumm_a)) deallocate(anumm_a)
      allocate(anumm_a(0:lag,0:lag,0:num_spec-1))
      
      if(allocated(anumm_lag)) deallocate(anumm_lag)
      allocate(anumm_lag(0:lag,0:lag))
      
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

      if (Z_eff .ne. 0) then

         !**********************************************************
         ! Compute collision operator with Laguerre base for
         ! eta-level positioning of flint
         !**********************************************************
         call init_collop(0, 0, 0d0, 0d0)
         call compute_source(asource, weightlag)
         call compute_collop_inf('e', 'e', m_ele, m_ele, 1d0, 1d0, anumm_aa(:,:,0,0), anumm_inf, &
              denmm_aa(:,:,0,0), ailmm_aa(:,:,:,0,0))
         anumm_lag(:,:) = anumm_aa(:,:,0,0) + Z_eff * anumm_inf(:,:)

         !**********************************************************
         ! Now compute collision operator with desired base
         !**********************************************************
         call init_collop(collop_base_prj, collop_base_exp, scalprod_alpha, scalprod_beta)

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

      else

         write (*,*) "Multispecies test mode."
         
         !**********************************************************
         ! Now compute collision operator with desired base
         !**********************************************************
         call init_collop(collop_base_prj, collop_base_exp, scalprod_alpha, scalprod_beta)

         !**********************************************************
         ! Compute sources
         !**********************************************************
         call compute_source(asource, weightlag)

         !**********************************************************
         ! Compute collision operator
         !**********************************************************
         call compute_collop('d', 'd', m_d, m_d, 1d0, 1d0, anumm_aa(:,:,0,0), &
              denmm_aa(:,:,0,0), ailmm_aa(:,:,:,0,0))
         call compute_collop('d', 'C', m_d, m_C, 1d0, 1d0, anumm_aa(:,:,0,1), &
              denmm_aa(:,:,0,1), ailmm_aa(:,:,:,0,1))
         call compute_collop('C', 'C', m_c, m_C, 1d0, 1d0, anumm_aa(:,:,1,1), &
              denmm_aa(:,:,1,1), ailmm_aa(:,:,:,1,1))
         call compute_collop('C', 'd', m_C, m_d, 1d0, 1d0, anumm_aa(:,:,1,0), &
              denmm_aa(:,:,1,0), ailmm_aa(:,:,:,1,0))
         
         !**********************************************************
         ! Sum up matrices
         !**********************************************************
         anumm_a = 0d0

         do a = 1, num_spec
            do b = 1, num_spec
               if (a .ne. b) then
                  anumm_a(:,:,a) = anumm_a(:,:,a) + anumm_aa(:,:,a,b)
                  denmm_a(:,:,a) = denmm_a(:,:,a) + denmm_aa(:,:,a,b)
               end if
            end do
         end do
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

      if (mpro%isMaster()) call write_collop('collop.h5')
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

    subroutine write_collop(h5filename)
      character(len=*) :: h5filename
      integer(HID_T)   :: h5id_collop, h5id_meta, h5id_species

      call h5_create(h5filename, h5id_collop)
      !call h5_define_group(h5id_collop, trim(tag_a) //'-'// trim(tag_b), h5id_species)
      call h5_define_group(h5id_collop, 'meta', h5id_meta)

      call h5_add(h5id_meta, 'lag', lag)
      call h5_add(h5id_meta, 'leg', leg)
      call h5_add(h5id_meta, 'scalprod_alpha', scalprod_alpha)
      call h5_add(h5id_meta, 'scalprod_beta',  scalprod_beta)
      !call h5_add(h5id_meta, 'm_a', m_a)
      !call h5_add(h5id_meta, 'm_b', m_b)
      !call h5_add(h5id_meta, 'T_a', T_a)
      !call h5_add(h5id_meta, 'T_b', T_b)
      call h5_add(h5id_meta, 'gamma_ab', gamma_ab)
      !call h5_add(h5id_meta, 'tag_a', tag_a)
      !call h5_add(h5id_meta, 'tag_b', tag_b)
      call h5_add(h5id_meta, 'collop_base_prj', collop_base_prj)
      call h5_add(h5id_meta, 'collop_base_exp', collop_base_exp)
      call h5_close_group(h5id_meta)
      
      call h5_add(h5id_collop, 'asource', asource, lbound(asource), ubound(asource))
      call h5_add(h5id_collop, 'anumm_a', anumm_a, lbound(anumm_a), ubound(anumm_a))
      call h5_add(h5id_collop, 'denmm_a', denmm_a, lbound(denmm_a), ubound(denmm_a))
      call h5_add(h5id_collop, 'anumm_aa', anumm_aa, lbound(anumm_aa), ubound(anumm_aa))
      call h5_add(h5id_collop, 'denmm_aa', denmm_aa, lbound(denmm_aa), ubound(denmm_aa))
      call h5_add(h5id_collop, 'ailmm_aa', ailmm_aa, lbound(ailmm_aa), ubound(ailmm_aa))
      call h5_add(h5id_collop, 'anumm', anumm, lbound(anumm), ubound(anumm))
      call h5_add(h5id_collop, 'anumm_lag', anumm_lag, lbound(anumm_lag), ubound(anumm_lag))
      call h5_add(h5id_collop, 'anumm_aa0', anumm_aa(:,:,0,0), lbound(anumm_aa(:,:,0,0)), ubound(anumm_aa(:,:,0,0)))
      call h5_add(h5id_collop, 'anumm_inf', anumm_inf, lbound(anumm_inf), ubound(anumm_inf))
      call h5_add(h5id_collop, 'denmm', denmm, lbound(denmm), ubound(denmm))
      call h5_add(h5id_collop, 'ailmm', ailmm, lbound(ailmm), ubound(ailmm))
      call h5_add(h5id_collop, 'weightlag', weightlag, lbound(weightlag), ubound(weightlag))
      call h5_add(h5id_collop, 'M_transform_', M_transform, lbound(M_transform), ubound(M_transform))
      call h5_add(h5id_collop, 'M_transform_inv', M_transform_inv, lbound(M_transform_inv), ubound(M_transform_inv))

      call h5_close(h5id_collop)

    end subroutine write_collop
    
end module collop
