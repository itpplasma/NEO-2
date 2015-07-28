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

  integer, private :: collop_unit = 300
  character(len=100), public         :: collop_file_Cm  = 'VectorCm.dat'
  integer, public                    :: collop_C_m_min = 0    
  integer, public                    :: collop_C_m_max = 0    
  real(kind=dp), allocatable, public :: collop_C(:)
  

  !**********************************************************
  ! Number of species
  !**********************************************************
  integer :: num_spec

  ! private only
  private full_file
  private full_file_1
  interface full_file
     module procedure full_file_1
  end interface full_file
  
  private check_unit
  private check_unit_1
  interface check_unit
     module procedure check_unit_1
  end interface check_unit
  
  contains

    ! ---------------------------------------------------------------------------
    ! private routines
    ! ---------------------------------------------------------------------------
    function full_file_1(name_in) result(name_out)
      character(len=100), intent(in)  :: name_in
      character(len=200) :: name_out
      
      write(name_out,'(200A)') trim(adjustl(collop_path)),trim(adjustl(name_in))
      
    end function full_file_1
    
    ! ---------------------------------------------------------------------------
    subroutine check_unit_1
      logical :: opened
      do
         inquire(unit=collop_unit,opened=opened)
         if(.not. opened) exit
         collop_unit = collop_unit + 1
      end do
    end subroutine check_unit_1
    
    subroutine collop_construct()     
      use nrtype, only : pi
      integer          :: openstatus, readstatus, k, m
      character(len=1) :: dummy

      num_spec = 2

      call check_unit      
      
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

      ! open
      open(unit=collop_unit,file=full_file(collop_file_Cm), & 
           status='old',action='read',iostat=openstatus)
      if (openstatus .ne. 0) then
         print *, 'Can not open file: ',full_file(collop_file_Cm)
         stop
      end if
      
      ! read initial lines
      do k = 1, dummy_read
         read(collop_unit,*,iostat=readstatus) dummy
         if (readstatus .ne. 0) then
            print *, 'Can not read file: ',full_file(collop_file_Cm),' (initial lines)'
            stop
         end if
      end do
      
      ! read upper size
      read(collop_unit,*,iostat=readstatus) collop_C_m_max
      if (readstatus .ne. 0) then
         print *, 'Can not read file: ',full_file(collop_file_Cm),' (size)'
         stop
      end if

      if (allocated(collop_C))  deallocate(collop_C)
      allocate(collop_C(collop_C_m_min:collop_C_m_max))

      ! read data
      do m = collop_C_m_min,collop_C_m_max
         read(collop_unit,*,iostat=readstatus) collop_C(m)
         if (readstatus .ne. 0) then
            print *, 'Can not read file: ',full_file(collop_file_Cm),' at m=',m
            stop
         end if
      end do

      weightlag=0.d0
      !    weightlag(1,0)=SQRT(6.d0)*pi/4.d0               !***change: 08.08.07 =>
      weightlag(1,0)=SQRT(6.d0*pi)/4.d0
      !    weightlag(2,0:lag)=collop_C(0:lag)              !***change: 19.09.07
      weightlag(2,0:lag)=collop_C(0:lag)/sqrt(pi)
      !    weightlag(3,0)=5.d0*SQRT(6.d0*pi)/16.d0         !***change: 19.09.07
      weightlag(3,0)=5.d0*SQRT(6.d0*pi)/8.d0
      !    IF(lag.GT.0) weightlag(3,1)=SQRT(15.d0*pi)/8.d0 !***change: 08.08.07 =>
      !    IF(lag.GT.0) weightlag(3,1)=-SQRT(15.d0*pi)/8.d0!***change: 19.09.07
      IF(lag.GT.0) weightlag(3,1)=-SQRT(15.d0*pi)/4.d0
      
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
      real(kind=dp), dimension(:), allocatable :: temp

      !**********************************************************
      ! Test configuration dd and dC
      !**********************************************************   
      call collop_set_species(0, .false.)
      call compute_collop('d', 'd' , m_d, m_d, 1d0, 1d0, asource, anumm, denmm, ailmm)
      
      call collop_set_species(1, .false.)
      call compute_collop('d', 'C' , m_d, m_C, 1d0, 1d0, asource, anumm, denmm, ailmm)

      allocate(temp(0:lag))
      temp = asource(:, 2)
      asource(:,2) = asource(:,3)
      asource(:,3) = temp
      deallocate(temp)
      
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
