MODULE collop
  USE rkstep_mod
  !use nrtype, only, private : pi, dp
  USE hdf5_tools
  USE collop_compute, ONLY : init_collop, &
       compute_source, compute_collop, gamma_ab, M_transform, M_transform_inv, &
       m_ele, m_d, m_C, compute_collop_inf
  USE mpiprovider_module

  IMPLICIT NONE
  
  !**********************************************************
  ! From old module, mainly for compatibility
  !**********************************************************
  INTEGER, PARAMETER :: dp = KIND(1.0d0)
  INTEGER, PARAMETER, PRIVATE :: dummy_read = 20
  REAL(kind=dp), PUBLIC       :: z_eff = 1.0_dp
  LOGICAL, PUBLIC             :: collop_talk      =  .TRUE. 
  LOGICAL, PUBLIC             :: collop_talk_much =  .TRUE.
  CHARACTER(len=100), PUBLIC  :: collop_path
  INTEGER                     :: collop_base_prj  = 0 ! 0...Laguerre (Default), 1...Polynomial
  INTEGER                     :: collop_base_exp  = 0 
  REAL(kind=dp)               :: scalprod_alpha = 0d0
  REAL(kind=dp)               :: scalprod_beta  = 0d0

  !**********************************************************
  ! Number of species
  !**********************************************************
  INTEGER :: num_spec

  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: anumm_inf

  CONTAINS
    
    SUBROUTINE collop_construct()     

    END SUBROUTINE collop_construct

    SUBROUTINE collop_set_species(ispec, opt_talk)
      INTEGER :: ispec
      LOGICAL, OPTIONAL :: opt_talk
      LOGICAL :: talk

      talk = .TRUE.
      IF (PRESENT(opt_talk)) talk = opt_talk
      
      !write (*,*) "Setting species to ", ispec

      !**********************************************************
      ! Switch collision operator matrices
      !**********************************************************
      anumm(0:lag, 0:lag) => anumm_a(:,:,ispec)      
      denmm(0:lag, 0:lag) => denmm_a(:,:,ispec)
      IF (Z_eff .NE. 0) THEN
         ailmm(0:lag, 0:lag, 0:leg) => ailmm_aa(:,:,:,ispec,ispec)
      ELSE
         ailmm(0:lag, 0:lag, 0:leg) => ailmm_aa(:,:,:,mpro%getRank(),ispec)
      END IF
      !**********************************************************
      ! Switch collisionality parameter
      !**********************************************************
      ! Not yet implemented
      
    END SUBROUTINE collop_set_species
    
    SUBROUTINE collop_load()
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: asource_temp
      REAL(kind=dp) :: alpha_temp, beta_temp
      INTEGER       :: a,b

      IF (Z_eff .NE. 0 ) THEN
         WRITE (*,*) "Standard mode."
         num_spec = 1
      ELSE
         WRITE (*,*) "Test mode for two species."
         num_spec = 2
      END IF
      
      !**********************************************************
      ! Allocation of matrices
      !**********************************************************
      IF(ALLOCATED(anumm_aa)) DEALLOCATE(anumm_aa)
      ALLOCATE(anumm_aa(0:lag,0:lag,0:num_spec-1,0:num_spec-1))

      IF(ALLOCATED(anumm_a)) DEALLOCATE(anumm_a)
      ALLOCATE(anumm_a(0:lag,0:lag,0:num_spec-1))
      
      IF(ALLOCATED(anumm_lag)) DEALLOCATE(anumm_lag)
      ALLOCATE(anumm_lag(0:lag,0:lag))
      
      IF(ALLOCATED(denmm_aa)) DEALLOCATE(denmm_aa)
      ALLOCATE(denmm_aa(0:lag,0:lag,0:num_spec-1,0:num_spec-1))

      IF(ALLOCATED(denmm_a)) DEALLOCATE(denmm_aa)
      ALLOCATE(denmm_a(0:lag,0:lag,0:num_spec-1))
      
      IF(ALLOCATED(asource)) DEALLOCATE(asource)
      ALLOCATE(asource(0:lag,3))
      
      IF(ALLOCATED(ailmm_aa)) DEALLOCATE(ailmm_aa)
      ALLOCATE(ailmm_aa(0:lag,0:lag,0:leg,0:num_spec-1,0:num_spec-1))
      
      IF(ALLOCATED(weightlag)) DEALLOCATE(weightlag)
      ALLOCATE(weightlag(3,0:lag))

      IF (ALLOCATED(anumm_inf)) DEALLOCATE(anumm_inf)
      ALLOCATE(anumm_inf(0:lag, 0:lag))

      IF (Z_eff .NE. 0) THEN

         !**********************************************************
         ! Compute collision operator with Laguerre base for
         ! eta-level positioning of flint
         !**********************************************************
         CALL init_collop(0, 0, 0d0, 0d0)
         CALL compute_source(asource, weightlag)
         CALL compute_collop_inf('e', 'e', m_ele, m_ele, 1d0, 1d0, anumm_aa(:,:,0,0), anumm_inf, &
              denmm_aa(:,:,0,0), ailmm_aa(:,:,:,0,0))
         anumm_lag(:,:) = anumm_aa(:,:,0,0) + Z_eff * anumm_inf(:,:)

         !**********************************************************
         ! Now compute collision operator with desired base
         !**********************************************************
         CALL init_collop(collop_base_prj, collop_base_exp, scalprod_alpha, scalprod_beta)

         !**********************************************************
         ! Compute sources
         !**********************************************************
         CALL compute_source(asource, weightlag)

         !**********************************************************
         ! Compute collision operator
         !**********************************************************
         CALL compute_collop_inf('e', 'e', m_ele, m_ele, 1d0, 1d0, anumm_aa(:,:,0,0), anumm_inf, &
              denmm_aa(:,:,0,0), ailmm_aa(:,:,:,0,0))

         !**********************************************************
         ! Sum up matrices
         !**********************************************************
         anumm_a(:,:,0) = anumm_aa(:,:,0,0) + Z_eff * anumm_inf(:,:)
         denmm_a(:,:,0) = denmm_aa(:,:,0,0)

      ELSE

         WRITE (*,*) "Multispecies test mode."
         
         !**********************************************************
         ! Now compute collision operator with desired base
         !**********************************************************
         CALL init_collop(collop_base_prj, collop_base_exp, scalprod_alpha, scalprod_beta)

         !**********************************************************
         ! Compute sources
         !**********************************************************
         CALL compute_source(asource, weightlag)

         !**********************************************************
         ! Compute collision operator
         !**********************************************************
         CALL compute_collop('d', 'd', m_d, m_d, 1d0, 1d0, anumm_aa(:,:,0,0), &
              denmm_aa(:,:,0,0), ailmm_aa(:,:,:,0,0))
         CALL compute_collop('d', 'C', m_d, m_C, 1d0, 1d0, anumm_aa(:,:,0,1), &
              denmm_aa(:,:,0,1), ailmm_aa(:,:,:,0,1))
         CALL compute_collop('C', 'C', m_c, m_C, 1d0, 1d0, anumm_aa(:,:,1,1), &
              denmm_aa(:,:,1,1), ailmm_aa(:,:,:,1,1))
         CALL compute_collop('C', 'd', m_C, m_d, 1d0, 1d0, anumm_aa(:,:,1,0), &
              denmm_aa(:,:,1,0), ailmm_aa(:,:,:,1,0))
         
         !**********************************************************
         ! Sum up matrices
         !**********************************************************
         anumm_a = 0d0
         denmm_a = 0d0
         DO a = 0, num_spec-1
            DO b = 0, num_spec-1
               !IF (a .NE. b) THEN
                  anumm_a(:,:,a) = anumm_a(:,:,a) + anumm_aa(:,:,a,b)
                  denmm_a(:,:,a) = denmm_a(:,:,a) + denmm_aa(:,:,a,b)
               !END IF
            END DO
         END DO
      END IF
     
      !**********************************************************
      ! Swap sources for NEO-2 convention
      !**********************************************************
      ALLOCATE(asource_temp(0:lag))
      asource_temp = asource(:, 2)
      asource(:,2) = asource(:,3)
      asource(:,3) = asource_temp

      asource_temp = weightlag(2,:)
      weightlag(2,:) = weightlag(3,:)
      weightlag(3,:) = asource_temp
      DEALLOCATE(asource_temp)
      
      !**********************************************************
      ! Set pointers to main species
      !**********************************************************
      CALL collop_set_species(0)

      !**********************************************************
      ! Write to screen
      !**********************************************************
      !write (*,*) asource
      !write (*,*) anumm
      !write (*,*) denmm
      !write (*,*) ailmm
      !write (*,*) weightlag

      !IF (mpro%isMaster()) CALL write_collop('collop.h5')
    END SUBROUTINE collop_load

    SUBROUTINE collop_unload()
      DEALLOCATE(anumm_aa)
      DEALLOCATE(denmm_aa)
      DEALLOCATE(anumm_a)
      DEALLOCATE(denmm_a)
      DEALLOCATE(asource)
      DEALLOCATE(weightlag)
    END SUBROUTINE collop_unload
  
    SUBROUTINE collop_deconstruct()
      
    END SUBROUTINE collop_deconstruct

    SUBROUTINE write_collop(h5filename)
      CHARACTER(len=*) :: h5filename
      INTEGER(HID_T)   :: h5id_collop, h5id_meta, h5id_species

      WRITE (*,*) "Rank: ", mpro%getRank()

      CALL h5_create(h5filename, h5id_collop)
      !call h5_define_group(h5id_collop, trim(tag_a) //'-'// trim(tag_b), h5id_species)
      CALL h5_define_group(h5id_collop, 'meta', h5id_meta)

      CALL h5_add(h5id_meta, 'lag', lag)
      CALL h5_add(h5id_meta, 'leg', leg)
      CALL h5_add(h5id_meta, 'scalprod_alpha', scalprod_alpha)
      CALL h5_add(h5id_meta, 'scalprod_beta',  scalprod_beta)
      !call h5_add(h5id_meta, 'm_a', m_a)
      !call h5_add(h5id_meta, 'm_b', m_b)
      !call h5_add(h5id_meta, 'T_a', T_a)
      !call h5_add(h5id_meta, 'T_b', T_b)
      CALL h5_add(h5id_meta, 'gamma_ab', gamma_ab)
      !call h5_add(h5id_meta, 'tag_a', tag_a)
      !call h5_add(h5id_meta, 'tag_b', tag_b)
      CALL h5_add(h5id_meta, 'collop_base_prj', collop_base_prj)
      CALL h5_add(h5id_meta, 'collop_base_exp', collop_base_exp)
      CALL h5_close_group(h5id_meta)
      
      CALL h5_add(h5id_collop, 'asource', asource, LBOUND(asource), UBOUND(asource))
      CALL h5_add(h5id_collop, 'anumm_a', anumm_a, LBOUND(anumm_a), UBOUND(anumm_a))
      CALL h5_add(h5id_collop, 'denmm_a', denmm_a, LBOUND(denmm_a), UBOUND(denmm_a))
      CALL h5_add(h5id_collop, 'anumm_aa', anumm_aa, LBOUND(anumm_aa), UBOUND(anumm_aa))
      CALL h5_add(h5id_collop, 'denmm_aa', denmm_aa, LBOUND(denmm_aa), UBOUND(denmm_aa))
      CALL h5_add(h5id_collop, 'ailmm_aa', ailmm_aa, LBOUND(ailmm_aa), UBOUND(ailmm_aa))
      CALL h5_add(h5id_collop, 'anumm', anumm, LBOUND(anumm), UBOUND(anumm))
      CALL h5_add(h5id_collop, 'anumm_lag', anumm_lag, LBOUND(anumm_lag), UBOUND(anumm_lag))
      CALL h5_add(h5id_collop, 'anumm_aa0', anumm_aa(:,:,0,0), LBOUND(anumm_aa(:,:,0,0)), UBOUND(anumm_aa(:,:,0,0)))
      CALL h5_add(h5id_collop, 'anumm_inf', anumm_inf, LBOUND(anumm_inf), UBOUND(anumm_inf))
      CALL h5_add(h5id_collop, 'denmm', denmm, LBOUND(denmm), UBOUND(denmm))
      CALL h5_add(h5id_collop, 'ailmm', ailmm, LBOUND(ailmm), UBOUND(ailmm))
      CALL h5_add(h5id_collop, 'weightlag', weightlag, LBOUND(weightlag), UBOUND(weightlag))
      CALL h5_add(h5id_collop, 'M_transform_', M_transform, LBOUND(M_transform), UBOUND(M_transform))
      CALL h5_add(h5id_collop, 'M_transform_inv', M_transform_inv, LBOUND(M_transform_inv), UBOUND(M_transform_inv))

      CALL h5_close(h5id_collop)

    END SUBROUTINE write_collop
    
END MODULE collop
