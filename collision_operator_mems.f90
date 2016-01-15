module collop
  use rkstep_mod
  !use nrtype, only, private : pi, dp
  use hdf5_tools
  use collop_compute, only : init_collop, &
       compute_source, compute_collop, gamma_ab, M_transform, M_transform_inv, &
       m_ele, m_d, m_C, compute_collop_inf, C_m, compute_xmmp, nu_D_hat, phi_exp, d_phi_exp, dd_phi_exp
  use mpiprovider_module
  ! WINNY
  use collisionality_mod, only : collpar,collpar_min,collpar_max, &
       v_max_resolution, v_min_resolution, phi_x_max, isw_lorentz

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
  real(kind=dp), dimension(:,:), allocatable :: x1mm, x2mm

  logical :: use_oldcollop = .false.
  
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
      
      if (allocated(anumm_inf)) deallocate(anumm_inf)
      allocate(anumm_inf(0:lag, 0:lag))

      if (Z_eff .ne. 0) then
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
         
         !**********************************************************
         ! Compute sources
         !**********************************************************
         call compute_source(asource, weightlag, weightden)
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
         call compute_source(asource, weightlag, weightden)

         !**********************************************************
         ! Compute x1mm and x2mm
         !**********************************************************
         call compute_xmmp(x1mm, x2mm)
         
         !**********************************************************
         ! Compute collision operator
         !**********************************************************
         call compute_collop('d', 'd', m_d, m_d, 1d0, 1d0, anumm_aa(:,:,0,0), &
              denmm_aa(:,:,0,0), ailmm_aa(:,:,:,0,0))
         call compute_collop('d', 'C', m_d, m_C, 1d0, 1d0, anumm_aa(:,:,0,1), &
              denmm_aa(:,:,0,1), ailmm_aa(:,:,:,0,1))
         call compute_collop('C', 'C', m_C, m_C, 1d0, 1d0, anumm_aa(:,:,1,1), &
              denmm_aa(:,:,1,1), ailmm_aa(:,:,:,1,1))
         call compute_collop('C', 'd', m_C, m_d, 1d0, 1d0, anumm_aa(:,:,1,0), &
              denmm_aa(:,:,1,0), ailmm_aa(:,:,:,1,0))
         
         !**********************************************************
         ! Sum up matrices
         !**********************************************************
         anumm_a = 0d0
         denmm_a = 0d0
         do a = 1, num_spec
            do b = 1, num_spec
               !if (a .ne. b) then
                  anumm_a(:,:,a) = anumm_a(:,:,a) + anumm_aa(:,:,a,b)
                  denmm_a(:,:,a) = denmm_a(:,:,a) + denmm_aa(:,:,a,b)
               !end if
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

      if (mpro%isMaster()) call write_collop()
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
    end subroutine collop_unload
  
    subroutine collop_deconstruct()
      
    end subroutine collop_deconstruct

    subroutine write_collop()
      integer(HID_T)   :: h5id_collop, h5id_meta, h5id_species
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
      call h5_add(h5id_collop, 'anumm_lag', anumm_lag, lbound(anumm_lag), ubound(anumm_lag))
      call h5_add(h5id_collop, 'anumm_aa0', anumm_aa(:,:,0,0), lbound(anumm_aa(:,:,0,0)), ubound(anumm_aa(:,:,0,0)))
      call h5_add(h5id_collop, 'anumm_inf', anumm_inf, lbound(anumm_inf), ubound(anumm_inf))
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
         write (f,*) NEW_line('A')
      end do
      close(f)    
      
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
