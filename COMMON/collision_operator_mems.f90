module collop
  use rkstep_mod
  !use nrtype, only, private : pi, dp
  use hdf5_tools
  use collop_compute, only : init_collop, &
       compute_source, compute_collop, gamma_ab, M_transform, M_transform_inv, &
       m_ele, m_d, m_C, m_alp, m_W, compute_collop_inf, C_m, compute_xmmp, &
       compute_collop_lorentz, nu_D_hat, phi_exp, d_phi_exp, dd_phi_exp, &
       compute_collop_rel
  use mpiprovider_module
  ! WINNY
  use collisionality_mod, only : collpar,collpar_min,collpar_max, &
       v_max_resolution, v_min_resolution, phi_x_max, isw_lorentz, conl_over_mfp, &
       isw_relativistic, T_e, &
       lsw_multispecies, isw_coul_log, num_spec, conl_over_mfp_spec, collpar_spec, &
       z_spec, m_spec, T_spec, n_spec
  use device_mod, only : device

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
  logical                     :: collop_only_precompute = .false.

  real(kind=dp), dimension(:,:), allocatable :: anumm_inf
  real(kind=dp), dimension(:,:), allocatable :: x1mm, x2mm

  integer, dimension(2) :: spec_ind_det, spec_ind_in

  contains
    
    subroutine collop_construct()     

    end subroutine collop_construct

    subroutine collop_set_species(ispec, opt_talk)
      integer :: ispec
      logical, optional :: opt_talk
      ! parameter
      real(kind=dp), parameter :: pi=3.14159265358979d0
      logical :: talk

      talk = .true.
      if (present(opt_talk)) talk = opt_talk
      
      !write (*,*) "Setting species to ", ispec

      !**********************************************************
      ! Switch collision operator matrices
      !**********************************************************
      anumm(0:lag, 0:lag) => anumm_a(:,:,ispec)      
      denmm(0:lag, 0:lag) => denmm_a(:,:,ispec)
      if (.not. lsw_multispecies) then 
         ailmm(0:lag, 0:lag, 0:leg) => ailmm_aa(:,:,:,ispec,ispec)
      else
         ailmm(0:lag, 0:lag, 0:leg) => ailmm_aa(:,:,:,mpro%getRank(),ispec)
      end if

      !! Modification by Andreas F. Martitsch (21.02.2017)
      ! ToDo: Delete this block - switch moved to flint_prepare for consistency!
!!$      !**********************************************************
!!$      ! Switch collisionality parameter
!!$      !**********************************************************
!!$      ! negative input for conl_over_mfp should provide collpar directly
!!$      if (.not. collop_only_precompute) then
!!$         conl_over_mfp=conl_over_mfp_spec(ispec)
!!$         if (conl_over_mfp .gt. 0.0d0) then
!!$            collpar=4.d0/(2.d0*pi*device%r0)*conl_over_mfp
!!$         else
!!$            collpar=-conl_over_mfp
!!$         end if
!!$      end if
      !! End Modification by Andreas F. Martitsch (21.02.2017)
      
    end subroutine collop_set_species
    
    subroutine collop_load()
      real(kind=dp), dimension(:), allocatable :: asource_temp

      real(kind=dp) :: alpha_temp, beta_temp
      integer       :: a,b,ispec

      real(kind=dp) :: taa_ov_tab_temp
      real(kind=dp) :: coll_a_temp, coll_b_temp
      real(kind=dp) :: za_temp, zb_temp
      
      if (.not. lsw_multispecies) then
         write (*,*) "Single species mode."
      else
         write (*,*) "Multispecies mode."
      end if

      ! species index
      ispec = mpro%getRank()
      
      !**********************************************************
      ! Allocation of matrices
      !**********************************************************
      if(allocated(Amm)) deallocate(Amm)
      allocate(Amm(0:lag,0:lag))

      if(allocated(anumm_aa)) deallocate(anumm_aa)
      allocate(anumm_aa(0:lag,0:lag,0:num_spec-1,0:num_spec-1))

      if(allocated(anumm_a)) deallocate(anumm_a)
      allocate(anumm_a(0:lag,0:lag,0:num_spec-1))
      
      if(allocated(anumm_lag)) deallocate(anumm_lag)
      allocate(anumm_lag(0:lag,0:lag))
      
      if(allocated(denmm_aa)) deallocate(denmm_aa)
      allocate(denmm_aa(0:lag,0:lag,0:num_spec-1,0:num_spec-1))

      if(allocated(denmm_a)) deallocate(denmm_a)
      allocate(denmm_a(0:lag,0:lag,0:num_spec-1))
      
      if(allocated(asource)) deallocate(asource)
      allocate(asource(0:lag,3))

      if(allocated(x1mm)) deallocate(x1mm)
      allocate(x1mm(0:lag,0:lag))

      if(allocated(x2mm)) deallocate(x2mm)
      allocate(x2mm(0:lag,0:lag))
      
      if(allocated(ailmm_aa)) deallocate(ailmm_aa)
      allocate(ailmm_aa(0:lag,0:lag,0:leg,0:num_spec-1,0:num_spec-1))
      
      if(allocated(weightlag)) deallocate(weightlag)
      allocate(weightlag(3,0:lag))

      if(allocated(weightden)) deallocate(weightden)
      allocate(weightden(0:lag))

      if(allocated(weightparflow)) deallocate(weightparflow)
      allocate(weightparflow(0:lag))

      if(allocated(weightenerg)) deallocate(weightenerg)
      allocate(weightenerg(0:lag))

      if (allocated(anumm_inf)) deallocate(anumm_inf)
      allocate(anumm_inf(0:lag, 0:lag))

      if (.not. lsw_multispecies) then
         !**********************************************************
         ! Now compute collision operator with desired base
         !**********************************************************
         call init_collop(collop_base_prj, collop_base_exp, scalprod_alpha, scalprod_beta)
         
         ! WINNY - for flint
         ! without any formula at the moment
         ! This 2.4 was anumm(1,1) 
         ! collpar_max = collpar * sqrt(2.4d0)
         ! I just used that collpar is proportional to 1/x with x = v/v_th
         ! v_max_resolution can be set in neo2.in in section collisions
         ! collpar_min = collpar_max / v_max_resolution**3

         ! New version with deflection frequency
         if (isw_lorentz .eq. 1) then
            collpar_min = collpar
            collpar_max = collpar
         else
            collpar_max = collpar * nu_D_hat(v_min_resolution)
            collpar_min = collpar * nu_D_hat(v_max_resolution)
         end if

         ! Non-relativistic Limit according to Trubnikov
         if (isw_relativistic .eq. 0) then

            !**********************************************************
            ! Compute sources
            !**********************************************************
            call compute_source(asource, weightlag, weightden, weightparflow, &
              weightenerg, Amm)
            write (*,*) "Weightden: ", weightden

            !**********************************************************
            ! Compute x1mm and x2mm
            !**********************************************************
            call compute_xmmp(x1mm, x2mm)

            !**********************************************************
            ! Compute collision operator
            !**********************************************************
            call compute_collop_inf('e', 'e', m_ele, m_ele, 1d0, 1d0, anumm_aa(:,:,0,0), anumm_inf, &
                 denmm_aa(:,:,0,0), ailmm_aa(:,:,:,0,0))

         ! Relativistic collision operator accordning to Braams and Karney
         elseif (isw_relativistic .ge. 1) then
            call compute_collop_rel(isw_relativistic, T_e, asource, weightlag, weightden, weightparflow, &
                 weightenerg, Amm, anumm_aa(:,:,0,0), anumm_inf, denmm_aa(:,:,0,0), ailmm_aa(:,:,:,0,0))
            !stop
         else
            write (*,*) "Relativistic switch ", isw_relativistic, " not defined."
         end if

         !**********************************************************
         ! Sum up matrices
         !**********************************************************
         anumm_a(:,:,0) = anumm_aa(:,:,0,0) + z_eff * anumm_inf(:,:)
         denmm_a(:,:,0) = denmm_aa(:,:,0,0)

         !**********************************************************
         ! Set pointers to main species
         !**********************************************************
         call collop_set_species(0)
         
      else

         write (*,*) "Multispecies test mode."
         
         !**********************************************************
         ! Now compute collision operator with desired base
         !**********************************************************
         call init_collop(collop_base_prj, collop_base_exp, scalprod_alpha, scalprod_beta)

         ! New version with deflection frequency 
         ! At the momement only for self-collisions !!!!!!!
         if (isw_lorentz .eq. 1) then
            print *,"collision_operator_mems.f90:&
                 &Lorentz collision model for multi-species mode not available!"
            stop
         else
            ! (ToDo: species-dependent eta-grid refinement - re-discretization)
            !collpar_max = collpar * nu_D_hat(v_min_resolution)
            !collpar_min = collpar * nu_D_hat(v_max_resolution)

            ! For testing make eta-grid species independent
            ! (collpar_max and collpar_min set for main species)
            write (*,*) "For testing make eta-grid species independent."
            collpar_max = collpar_spec(0) * nu_D_hat(v_min_resolution)
            collpar_min = collpar_spec(0) * nu_D_hat(v_max_resolution)
         end if
         
         !**********************************************************
         ! Compute sources
         !**********************************************************
         call compute_source(asource, weightlag, weightden, weightparflow, &
              weightenerg, Amm)

         !**********************************************************
         ! Compute x1mm and x2mm
         !**********************************************************
         call compute_xmmp(x1mm, x2mm)
         
         !**********************************************************
         ! Compute collision operator
         !**********************************************************
         ! old: without MPI parallelization
!!$         do a = 0, num_spec-1
!!$            do b = 0, num_spec-1
!!$               
!!$               spec_ind_in = (/ a, b /)
!!$               call collop_exists(m_spec(a), m_spec(b), T_spec(a), T_spec(b), &
!!$                    spec_ind_in, spec_ind_det)
!!$
!!$               if ( all(spec_ind_det .eq. -1) ) then
!!$                  !print *,'a, b, spec_ind_det: ',a,b,spec_ind_det
!!$                  !print *,'Do the computation.'
!!$                  call compute_collop('a', 'b', m_spec(a), m_spec(b), T_spec(a), T_spec(b), &
!!$                       anumm_aa(:,:,a,b), denmm_aa(:,:,a,b), ailmm_aa(:,:,:,a,b))
!!$               else
!!$                  !print *,'a, b, spec_ind_det: ',a,b,spec_ind_det
!!$                  !print *,'Load matrices.'
!!$                  anumm_aa(:,:,a,b) = anumm_aa(:,:,spec_ind_det(1),spec_ind_det(2))
!!$                  denmm_aa(:,:,a,b) = denmm_aa(:,:,spec_ind_det(1),spec_ind_det(2))
!!$                  ailmm_aa(:,:,:,a,b) = ailmm_aa(:,:,:,spec_ind_det(1),spec_ind_det(2))
!!$               end if
!!$               
!!$            end do
!!$         end do
         ! new: with MPI parallelization
         b = mpro%getRank()
         do a = 0, num_spec-1
            spec_ind_in = (/ a, b /)
            call collop_exists_mpi(m_spec(a), T_spec(a), spec_ind_in, spec_ind_det)
            if ( all(spec_ind_det .eq. -1) ) then
               !print *,'a, b, spec_ind_det: ',a,b,spec_ind_det
               !print *,'Do the computation.'
               call compute_collop('a', 'b', m_spec(a), m_spec(b), T_spec(a), T_spec(b), &
                    anumm_aa(:,:,a,b), denmm_aa(:,:,a,b), ailmm_aa(:,:,:,a,b))
            else
               !print *,'a, b, spec_ind_det: ',a,b,spec_ind_det
               !print *,'Load matrices.'
               anumm_aa(:,:,a,b) = anumm_aa(:,:,spec_ind_det(1),spec_ind_det(2))
               denmm_aa(:,:,a,b) = denmm_aa(:,:,spec_ind_det(1),spec_ind_det(2))
               ailmm_aa(:,:,:,a,b) = ailmm_aa(:,:,:,spec_ind_det(1),spec_ind_det(2))
            end if
         end do
         call mpro%allgather(anumm_aa(:,:,:,b),anumm_aa)
         call mpro%allgather(denmm_aa(:,:,:,b),denmm_aa)
         call mpro%allgather(ailmm_aa(:,:,:,:,b),ailmm_aa)

         !call compute_collop('d', 'd', m_d, m_d, 1d0, 1d0, anumm_aa(:,:,0,0), &
         !     denmm_aa(:,:,0,0), ailmm_aa(:,:,:,0,0))
         
         !call compute_collop('d', 'C', m_d, m_C, 1d0, 1d0, anumm_aa(:,:,0,1), &
         !     denmm_aa(:,:,0,1), ailmm_aa(:,:,:,0,1))
         !call compute_collop('C', 'C', m_C, m_C, 1d0, 1d0, anumm_aa(:,:,1,1), &
         !     denmm_aa(:,:,1,1), ailmm_aa(:,:,:,1,1))
         !call compute_collop('C', 'd', m_C, m_d, 1d0, 1d0, anumm_aa(:,:,1,0), &
         !     denmm_aa(:,:,1,0), ailmm_aa(:,:,:,1,0))
         
         !call compute_collop('d', 'alp', m_d, m_alp, 1d0, 1d0, anumm_aa(:,:,0,1), &
         !     denmm_aa(:,:,0,1), ailmm_aa(:,:,:,0,1))
         !call compute_collop('alp', 'alp', m_alp, m_alp, 1d0, 1d0, anumm_aa(:,:,1,1), &
         !     denmm_aa(:,:,1,1), ailmm_aa(:,:,:,1,1))
         !call compute_collop('alp', 'd', m_alp, m_d, 1d0, 1d0, anumm_aa(:,:,1,0), &
         !     denmm_aa(:,:,1,0), ailmm_aa(:,:,:,1,0))
         
         !call compute_collop('d', 'W', m_d, m_W, 1d0, 1d0, anumm_aa(:,:,0,1), &
         !     denmm_aa(:,:,0,1), ailmm_aa(:,:,:,0,1))
         !call compute_collop('W', 'W', m_W, m_W, 1d0, 1d0, anumm_aa(:,:,1,1), &
         !     denmm_aa(:,:,1,1), ailmm_aa(:,:,:,1,1))
         !call compute_collop('W', 'd', m_W, m_d, 1d0, 1d0, anumm_aa(:,:,1,0), &
         !     denmm_aa(:,:,1,0), ailmm_aa(:,:,:,1,0))

         !!m_ele = m_ele*0.5d0
         !call compute_collop('d', 'e', m_d, m_ele, 1d0, 1d0, anumm_aa(:,:,0,1), &
         !     denmm_aa(:,:,0,1), ailmm_aa(:,:,:,0,1))
         !call compute_collop('e', 'e', m_ele, m_ele, 1d0, 1d0, anumm_aa(:,:,1,1), &
         !     denmm_aa(:,:,1,1), ailmm_aa(:,:,:,1,1))
         !call compute_collop('e', 'd', m_ele, m_d, 1d0, 1d0, anumm_aa(:,:,1,0), &
         !     denmm_aa(:,:,1,0), ailmm_aa(:,:,:,1,0))

         !call compute_collop('W', 'e', m_W, m_ele, 1d0, 1d0, anumm_aa(:,:,0,1), &
         !     denmm_aa(:,:,0,1), ailmm_aa(:,:,:,0,1))
         !call compute_collop('e', 'e', m_ele, m_ele, 1d0, 1d0, anumm_aa(:,:,1,1), &
         !     denmm_aa(:,:,1,1), ailmm_aa(:,:,:,1,1))
         !call compute_collop('e', 'W', m_ele, m_W, 1d0, 1d0, anumm_aa(:,:,1,0), &
         !     denmm_aa(:,:,1,0), ailmm_aa(:,:,:,1,0))
         !stop
         
         !**********************************************************
         ! Sum up matrices
         !**********************************************************
         anumm_a = 0d0
         denmm_a = 0d0
         !print *,num_spec
         !stop
         do a = 0, num_spec-1
            coll_a_temp = collpar_spec(a)
            za_temp = z_spec(a)
            do b = 0, num_spec-1
               coll_b_temp = collpar_spec(b)
               zb_temp = z_spec(b)
               ! isw_coul_log = 0: Coulomb logarithm set as species independent
               !                   (overrides values for n_spec)
               ! isw_coul_log = 1: Coulomb logarithm computed for each species
               !                   using n_spec, T_spec
               !                   (overrides values for collisionality parameters)
               if (isw_coul_log .eq. 0) then
                  taa_ov_tab_temp = (coll_b_temp/coll_a_temp) * &
                       (((T_spec(b)*za_temp)/(T_spec(a)*zb_temp))**2)
               else
                  print *,"collision_operator_mems.f90: &
                       &species-dependent Coulomb logarithm not yet implemented!"
                  print *,"Please use switch isw_coul_log = 0."
                  stop
               end if
               !print *,'taa_ov_tab: ',a,b,taa_ov_tab_temp
               !print *,coll_a_temp,coll_b_temp
               !print *,za_temp,zb_temp
               anumm_a(:,:,a) = &
                    anumm_a(:,:,a) + anumm_aa(:,:,a,b) * taa_ov_tab_temp
               denmm_a(:,:,a) = &
                    denmm_a(:,:,a) + denmm_aa(:,:,a,b) * taa_ov_tab_temp
               ailmm_aa(:,:,:,a,b) = &
                    ailmm_aa(:,:,:,a,b) * taa_ov_tab_temp
            end do
         end do

         !**********************************************************
         ! Set pointers to main species
         !**********************************************************
         call collop_set_species(ispec)
         
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
      ! Write to screen
      !**********************************************************
      !write (*,*) asource
      !write (*,*) anumm
      !write (*,*) denmm
      !write (*,*) ailmm
      !write (*,*) weightlag

      if (mpro%isMaster()) call write_collop()

    contains
      !
      subroutine collop_exists(ma,mb,Ta,Tb,ind_in,ind_out)
        ! input / output
        real(kind=dp), intent(in) :: ma, mb, Ta, Tb
        integer, dimension(2), intent(in)  :: ind_in
        integer, dimension(2), intent(out) :: ind_out
        ! internal
        integer :: a, b
        !
        ! ind_out = -1: matrix elements not available for
        !               species (ma,Ta;mb,Tb) - do the computation 
        ind_out = -1
        do a=0,ind_in(1)
           if ( (ma.ne.m_spec(a)) .or. (Ta.ne.T_spec(a)) ) cycle
           do b=0,ind_in(2)
              if ( (mb.ne.m_spec(b)) .or. (Tb.ne.T_spec(b)) ) cycle
              ind_out=(/a,b/)
              !print *,'ind_in: ',ind_in
              !print *,'ind_out: ',ind_out
              if ( all(ind_out.eq.ind_in) ) then
                 ind_out = -1
              end if
              return
           end do
        end do
        !
      end subroutine collop_exists
      !
      subroutine collop_exists_mpi(ma,Ta,ind_in,ind_out)
        ! input / output
        real(kind=dp), intent(in) :: ma, Ta
        integer, dimension(2), intent(in)  :: ind_in
        integer, dimension(2), intent(out) :: ind_out
        ! internal
        integer :: a, b
        !
        ! ind_out = -1: matrix elements not available for
        !               species (ma,Ta;mb,Tb) - do the computation 
        ind_out = -1
        do a=0,ind_in(1)
           if ( (ma.ne.m_spec(a)) .or. (Ta.ne.T_spec(a)) ) cycle
           ind_out=(/a,ind_in(2)/)
           !print *,'ind_in: ',ind_in
           !print *,'ind_out: ',ind_out
           if ( all(ind_out.eq.ind_in) ) then
              ind_out = -1
           end if
           return
        end do
        !
      end subroutine collop_exists_mpi
      !
    end subroutine collop_load

    subroutine collop_unload()
      deallocate(anumm_aa)
      deallocate(denmm_aa)
      deallocate(anumm_a)
      deallocate(denmm_a)
      deallocate(asource)
      deallocate(weightlag)
      deallocate(x1mm)
      deallocate(x2mm)
      deallocate(conl_over_mfp_spec)
    end subroutine collop_unload
  
    subroutine collop_deconstruct()
      
    end subroutine collop_deconstruct

    subroutine write_collop()
      integer(HID_T)   :: h5id_collop, h5id_meta
      integer          :: m, mp, l, xi, n_x
      integer          :: f = 4234
      real(kind=dp), dimension(:), allocatable :: x
      real(kind=dp), dimension(:,:), allocatable :: phi_x, dphi_x, ddphi_x

      call h5_create('collop.h5', h5id_collop)
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
      !call h5_add(h5id_collop, 'anumm_lag', anumm_lag, lbound(anumm_lag), ubound(anumm_lag))
      call h5_add(h5id_collop, 'anumm_aa0', anumm_aa(:,:,0,0), lbound(anumm_aa(:,:,0,0)), ubound(anumm_aa(:,:,0,0)))
      if (z_eff .ne. 0) call h5_add(h5id_collop, 'anumm_inf', anumm_inf, lbound(anumm_inf), ubound(anumm_inf))
      call h5_add(h5id_collop, 'denmm', denmm, lbound(denmm), ubound(denmm))
      call h5_add(h5id_collop, 'ailmm', ailmm, lbound(ailmm), ubound(ailmm))
      call h5_add(h5id_collop, 'weightlag', weightlag, lbound(weightlag), ubound(weightlag))
      call h5_add(h5id_collop, 'M_transform_', M_transform, lbound(M_transform), ubound(M_transform))
      call h5_add(h5id_collop, 'M_transform_inv', M_transform_inv, lbound(M_transform_inv), ubound(M_transform_inv))
      call h5_add(h5id_collop, 'C_m', C_m, lbound(C_m), ubound(C_m))

      !**********************************************************
      ! Write test functions
      !**********************************************************
      n_x = 999
      allocate(x(n_x))
      allocate(phi_x(0:lag, 1:n_x))
      allocate(dphi_x(0:lag, 1:n_x))
      allocate(ddphi_x(0:lag, 1:n_x))
      
      do m = 0, lag
         do xi = 1, n_x
            x(xi) = 10d0/(n_x-1) * (xi-1) 
            phi_x(m,xi)   = phi_exp(m,x(xi))
            dphi_x(m,xi)  = d_phi_exp(m,x(xi))
            ddphi_x(m,xi) = dd_phi_exp(m,x(xi))
            !write (*,*) x(xi), m, phi_exp(m,x(xi)), d_phi_exp(m,x(xi)), dd_phi_exp(m,x(xi))
         end do
      end do

      call h5_add(h5id_collop, 'x', x, lbound(x), ubound(x))
      call h5_add(h5id_collop, 'phi_x', phi_x, lbound(phi_x), ubound(phi_x))
      call h5_add(h5id_collop, 'dphi_x', dphi_x, lbound(dphi_x), ubound(dphi_x))
      call h5_add(h5id_collop, 'ddphi_x', ddphi_x, lbound(ddphi_x), ubound(ddphi_x))
      
      call h5_close(h5id_collop)

      !**********************************************************
      ! ASCII
      !**********************************************************

      open(f, file='SourceAa123m_Cm.dat', status='replace')
      do m=0,19
         write (f,'(A)') '!' 
      end do
      write (f,'(I0)') lag
      do m=0,lag
         write (f,'(4(es23.15E02))') asource(m, 1), asource(m, 3), asource(m, 2), C_m(m)
      end do
      close(f)

      open(f, file='MatrixNu_mmp-gee.dat', status='replace')
      do m=0,19
         write (f,'(A)') '!' 
      end do
      write (f,'(I0)') lag
      write (f,'(I0)') lag
      do m=0,lag
         do mp=0,lag
            write (f,'(1(es23.15E02))', advance='NO') anumm_aa(m, mp, 0, 0)
         end do
         write (f, '(A)', advance='NO') NEW_line('A')
      end do
      close(f)    

      if (z_eff .ne. 0) then
         open(f, file='MatrixNu_mmp-ginf.dat', status='replace')
         do m=0,19
            write (f,'(A)') '!' 
         end do
         write (f,'(I0)') lag
         write (f,'(I0)') lag
         do m=0,lag
            do mp=0,lag
               write (f,'(1(es23.15E02))', advance='NO') anumm_inf(m, mp)
            end do
            write (f, '(A)', advance='NO') NEW_line('A')
         end do
         close(f)    
      end if
   
      open(f, file='MatrixD_mmp-gee.dat', status='replace')
      do m=0,19
         write (f,'(A)') '!' 
      end do
      write (f,'(I0)') lag
      write (f,'(I0)') lag
      do m=0,lag
         do mp=0,lag
            write (f,'(1(es23.15E02))', advance='NO') denmm(m, mp)
         end do
         write (f, '(A)', advance='NO') NEW_line('A')
      end do
      close(f)

      open(f, file='MatrixI_mmp-gee.dat', status='replace')
      do m=0,19
         write (f,'(A)') '!' 
      end do
      write (f,'(I0)') lag
      write (f,'(I0)') lag
      write (f,'(I0)') leg
      do l = 0,leg
         write (f,'(I0)') l
         do m=0,lag
            do mp=0,lag
               write (f,'(1(es23.15E02))', advance='NO') ailmm_aa(m, mp, l, 0, 0)
            end do
            write (f, '(A)', advance='NO') NEW_line('A')
         end do
      end do
      close(f)    

      !stop
    end subroutine write_collop
    
end module collop
