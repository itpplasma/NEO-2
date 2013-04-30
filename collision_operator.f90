MODULE collop
  ! Changelog:
  !
  ! New files from Georg are treated now (Winny, 30.03.2006)
  !
  ! additional:   collop_C
  ! changed name: collop_nu is now collop_nu_gee
  ! additional:   collop_nu_ginf
  ! changed name: collop_D is now collop_D_gee
  ! changed name: collop_I is now collop_I_gee

  ! Module to compute various matrix elements for the Coulomb collision operator
  !
  ! Public routines:
  !  collop_construct   (reading, allocating)
  !  collop_deconstruct (deallocating)
  !
  ! Public variables:
  !
  ! A(m,k), k = 1,2,3
  !  collop_A,collop_A_m_min,collop_A_m_max
  ! A1(m), A2(m), A3(m) Alternative
  !  collop_A1,collop_A2,collop_A3,collop_A_m_min,collop_A_m_max
  ! C(m)
  !  collop_C,collop_C_m_min,collop_C_m_max
  ! nu_gee(m,mp), nu_ginf(m,mp), 
  !  collop_nu_gee,collop_nu_inf
  !  collop_nu_m_min,collop_nu_m_max,collop_nu_mp_min,collop_nu_mp_max
  ! D_gee(m,mp)
  !  collop_D_gee
  !  collop_D_m_min,collop_D_m_max,collop_D_mp_min,collop_D_mp_max
  ! I_gee(m,mp,lp)
  !  collop_I_gee
  !  collop_I_m_min,collop_I_m_max,collop_I_mp_min,collop_I_mp_max,collop_I_lp_min,collop_I_lp_max
  !
  ! Pathname (Default values are set)
  !  collop_path
  ! Filenames (Default values are set)
  !  collop_file_A
  !  collop_file_nu_gee
  !  collop_file_nu_ginf
  !  collop_file_D
  !  collop_file_I

  IMPLICIT NONE
  
  ! ---------------------------------------------------------------------------
  ! private parameters
  INTEGER, PARAMETER, PRIVATE :: dp = KIND(1.0d0)
  INTEGER, PARAMETER, PRIVATE :: dummy_read = 20
  
  ! ---------------------------------------------------------------------------
  ! public variables
  !
  REAL(kind=dp), PUBLIC :: z_eff = 1.0_dp

  ! talk or be silent (not used at the moment)
  LOGICAL, PUBLIC :: collop_talk      =  .TRUE. 
  LOGICAL, PUBLIC :: collop_talk_much =  .TRUE. 

  ! path to data
#if defined(VSC2)
  CHARACTER(len=100), PUBLIC :: collop_path = '/home/lv70337/gernot_k/Neo2/data-MatrixElements/'
#else
#if defined(ZID)
  CHARACTER(len=100), PUBLIC :: collop_path = /home/gernot_k/Neo2/data-MatrixElements/'
#else
  ! This is ITP-Cluster
  CHARACTER(len=100), PUBLIC :: collop_path = '/afs/itp.tugraz.at/proj/plasma/DOCUMENTS/Neo2/data-MatrixElements/'
#endif
#endif

! file names
  CHARACTER(len=100), PUBLIC :: collop_file_A        = 'SourceAa123m_Cm.dat'
  CHARACTER(len=100), PUBLIC :: collop_file_nu_gee   = 'MatrixNu_mmp-gee.dat'
  CHARACTER(len=100), PUBLIC :: collop_file_nu_ginf  = 'MatrixNu_mmp-ginf.dat'
  CHARACTER(len=100), PUBLIC :: collop_file_D        = 'MatrixD_mmp-gee.dat'
  CHARACTER(len=100), PUBLIC :: collop_file_I        = 'MatrixI_mmp-gee.dat'


  ! data arrays
  ! source terms
  ! either 
  !        collop_A1(:),collop_A2(:),collop_A3(:)
  ! or
  !        collop_A(:,1),collop_A(:,2),collop_A(:,3)
  !
  INTEGER, PUBLIC                    :: collop_A_m_min = 0    
  INTEGER, PUBLIC                    :: collop_A_m_max = 0    
  REAL(kind=dp), ALLOCATABLE, PUBLIC :: collop_A(:,:)
  REAL(kind=dp), ALLOCATABLE, PUBLIC :: collop_A1(:)
  REAL(kind=dp), ALLOCATABLE, PUBLIC :: collop_A2(:)
  REAL(kind=dp), ALLOCATABLE, PUBLIC :: collop_A3(:)
  INTEGER, PUBLIC                    :: collop_C_m_min = 0    
  INTEGER, PUBLIC                    :: collop_C_m_max = 0    
  REAL(kind=dp), ALLOCATABLE, PUBLIC :: collop_C(:)

  ! nu(m,mp)
  INTEGER, PUBLIC                    :: collop_nu_m_min  = 0      
  INTEGER, PUBLIC                    :: collop_nu_m_max  = 0      
  INTEGER, PUBLIC                    :: collop_nu_mp_min = 0      
  INTEGER, PUBLIC                    :: collop_nu_mp_max = 0      
  REAL(kind=dp), ALLOCATABLE, PUBLIC :: collop_nu_gee(:,:)
  REAL(kind=dp), ALLOCATABLE, PUBLIC :: collop_nu_ginf(:,:)

  ! D(m,mp)
  INTEGER, PUBLIC                    :: collop_D_m_min  = 0      
  INTEGER, PUBLIC                    :: collop_D_m_max  = 0      
  INTEGER, PUBLIC                    :: collop_D_mp_min = 0      
  INTEGER, PUBLIC                    :: collop_D_mp_max = 0      
  REAL(kind=dp), ALLOCATABLE, PUBLIC :: collop_D_gee(:,:)

  ! I(m,mp,lp)
  INTEGER, PUBLIC                    :: collop_I_m_min  = 0      
  INTEGER, PUBLIC                    :: collop_I_m_max  = 0      
  INTEGER, PUBLIC                    :: collop_I_mp_min = 0      
  INTEGER, PUBLIC                    :: collop_I_mp_max = 0      
  INTEGER, PUBLIC                    :: collop_I_lp_min = 0      
  INTEGER, PUBLIC                    :: collop_I_lp_max = 0      
  REAL(kind=dp), ALLOCATABLE, PUBLIC :: collop_I_gee(:,:,:)
  ! ---------------------------------------------------------------------------
  ! private variables 
  !
  ! unit
  INTEGER, PRIVATE :: collop_unit = 300

  ! ---------------------------------------------------------------------------
  ! interfaces

  ! public
  PUBLIC collop_construct
  PRIVATE collop_construct_all
  INTERFACE collop_construct
     MODULE PROCEDURE collop_construct_all
  END INTERFACE

  PUBLIC collop_deconstruct
  PRIVATE collop_deconstruct_all
  INTERFACE collop_deconstruct
     MODULE PROCEDURE collop_deconstruct_all
  END INTERFACE

  PUBLIC collop_load
  PRIVATE collop_load_all
  INTERFACE collop_load
     MODULE PROCEDURE collop_load_all
  END INTERFACE

  PUBLIC collop_unload
  PRIVATE collop_unload_all
  INTERFACE collop_unload
     MODULE PROCEDURE collop_unload_all
  END INTERFACE

  ! private only
  PRIVATE full_file
  PRIVATE full_file_1
  INTERFACE full_file
     MODULE PROCEDURE full_file_1
  END INTERFACE

  PRIVATE check_unit
  PRIVATE check_unit_1
  INTERFACE check_unit
     MODULE PROCEDURE check_unit_1
  END INTERFACE


CONTAINS

  ! ---------------------------------------------------------------------------
  ! public routines
  ! ---------------------------------------------------------------------------
  SUBROUTINE collop_construct_all
    INTEGER :: k,m,ll,lp
    INTEGER :: openstatus,readstatus
    CHARACTER(len=1) :: dummy

    CALL check_unit

    ! -------------------------------------------------------------------------
    ! A(m,k) k = 1,2,3; C(m)
    !
    ! deallocate
    IF (ALLOCATED(collop_A))  DEALLOCATE(collop_A)
    IF (ALLOCATED(collop_A1)) DEALLOCATE(collop_A1)
    IF (ALLOCATED(collop_A2)) DEALLOCATE(collop_A2)
    IF (ALLOCATED(collop_A3)) DEALLOCATE(collop_A3)
    IF (ALLOCATED(collop_C))  DEALLOCATE(collop_C)
    ! open
    OPEN(unit=collop_unit,file=full_file(collop_file_A), & 
         status='old',action='read',iostat=openstatus)
    IF (openstatus .NE. 0) THEN
       PRINT *, 'Can not open file: ',full_file(collop_file_A)
       STOP
    END IF
    ! read initial lines
    DO k = 1, dummy_read
       READ(collop_unit,*,iostat=readstatus) dummy
       IF (readstatus .NE. 0) THEN
          PRINT *, 'Can not read file: ',full_file(collop_file_A),' (initial lines)'
          STOP
       END IF
    END DO
    ! read upper size
    READ(collop_unit,*,iostat=readstatus) collop_A_m_max
    IF (readstatus .NE. 0) THEN
       PRINT *, 'Can not read file: ',full_file(collop_file_A),' (size)'
       STOP
    END IF
    collop_C_m_max = collop_A_m_max
    ! allocate
    ALLOCATE(collop_A(collop_A_m_min:collop_A_m_max,3))
    ALLOCATE(collop_A1(collop_A_m_min:collop_A_m_max))
    ALLOCATE(collop_A2(collop_A_m_min:collop_A_m_max))
    ALLOCATE(collop_A3(collop_A_m_min:collop_A_m_max))
    ALLOCATE(collop_C(collop_C_m_min:collop_C_m_max))
    ! read data
    DO m = collop_A_m_min,collop_A_m_max
       READ(collop_unit,*,iostat=readstatus) collop_A1(m),collop_A2(m),collop_A3(m),collop_C(m)
       IF (readstatus .NE. 0) THEN
          PRINT *, 'Can not read file: ',full_file(collop_file_A),' at m=',m
          STOP
       END IF
       collop_A(m,:) = (/ collop_A1(m),collop_A2(m),collop_A3(m) /)
    END DO
    ! close
    CLOSE(unit=collop_unit)
    ! -------------------------------------------------------------------------
    ! end A(m,k) and C(m)
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! nu_gee(m,mp)
    !
    ! deallocate
    IF (ALLOCATED(collop_nu_gee))  DEALLOCATE(collop_nu_gee)
    ! open
    OPEN(unit=collop_unit,file=full_file(collop_file_nu_gee), & 
         status='old',action='read',iostat=openstatus)
    IF (openstatus .NE. 0) THEN
       PRINT *, 'Can not open file: ',full_file(collop_file_nu_gee)
       STOP
    END IF
    ! read initial lines
    DO k = 1, dummy_read
       READ(collop_unit,*,iostat=readstatus) dummy
       IF (readstatus .NE. 0) THEN
          PRINT *, 'Can not read file: ',full_file(collop_file_nu_gee),' (initial lines)'
          STOP
       END IF
    END DO
    ! read upper sizes
    READ(collop_unit,*,iostat=readstatus) collop_nu_m_max
    IF (readstatus .NE. 0) THEN
       PRINT *, 'Can not read file: ',full_file(collop_file_nu_gee),' (size)'
       STOP
    END IF
    READ(collop_unit,*,iostat=readstatus) collop_nu_mp_max
    IF (readstatus .NE. 0) THEN
       PRINT *, 'Can not read file: ',full_file(collop_file_nu_gee),' (size)'
       STOP
    END IF
    ! allocate
    ALLOCATE(collop_nu_gee(collop_nu_m_min:collop_nu_m_max,collop_nu_mp_min:collop_nu_mp_max))
    ! read data
    DO m = collop_nu_m_min,collop_nu_m_max
       READ(collop_unit,*,iostat=readstatus) collop_nu_gee(m,:)
       IF (readstatus .NE. 0) THEN
          PRINT *, 'Can not read file: ',full_file(collop_file_nu_gee),' at m=',m
          STOP
       END IF
    END DO
    ! close
    CLOSE(unit=collop_unit)
    ! -------------------------------------------------------------------------
    ! end nu_gee
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! nu_ginf(m,mp)
    !
    ! deallocate
    IF (ALLOCATED(collop_nu_ginf))  DEALLOCATE(collop_nu_ginf)
    ! open
    OPEN(unit=collop_unit,file=full_file(collop_file_nu_ginf), & 
         status='old',action='read',iostat=openstatus)
    IF (openstatus .NE. 0) THEN
       PRINT *, 'Can not open file: ',full_file(collop_file_nu_ginf)
       STOP
    END IF
    ! read initial lines
    DO k = 1, dummy_read
       READ(collop_unit,*,iostat=readstatus) dummy
       IF (readstatus .NE. 0) THEN
          PRINT *, 'Can not read file: ',full_file(collop_file_nu_ginf),' (initial lines)'
          STOP
       END IF
    END DO
    ! read upper sizes
    READ(collop_unit,*,iostat=readstatus) collop_nu_m_max
    IF (readstatus .NE. 0) THEN
       PRINT *, 'Can not read file: ',full_file(collop_file_nu_ginf),' (size)'
       STOP
    END IF
    READ(collop_unit,*,iostat=readstatus) collop_nu_mp_max
    IF (readstatus .NE. 0) THEN
       PRINT *, 'Can not read file: ',full_file(collop_file_nu_ginf),' (size)'
       STOP
    END IF
    ! allocate
    ALLOCATE(collop_nu_ginf(collop_nu_m_min:collop_nu_m_max,collop_nu_mp_min:collop_nu_mp_max))
    ! read data
    DO m = collop_nu_m_min,collop_nu_m_max
       READ(collop_unit,*,iostat=readstatus) collop_nu_ginf(m,:)
       IF (readstatus .NE. 0) THEN
          PRINT *, 'Can not read file: ',full_file(collop_file_nu_ginf),' at m=',m
          STOP
       END IF
    END DO
    ! close
    CLOSE(unit=collop_unit)
    ! -------------------------------------------------------------------------
    ! end nu_ginf
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! D_gee(m,mp)
    !
    ! deallocate
    IF (ALLOCATED(collop_D_gee))  DEALLOCATE(collop_D_gee)
    ! open
    OPEN(unit=collop_unit,file=full_file(collop_file_D), & 
         status='old',action='read',iostat=openstatus)
    IF (openstatus .NE. 0) THEN
       PRINT *, 'Can not open file: ',full_file(collop_file_D)
       STOP
    END IF
    ! read initial lines
    DO k = 1, dummy_read
       READ(collop_unit,*,iostat=readstatus) dummy
       IF (readstatus .NE. 0) THEN
          PRINT *, 'Can not read file: ',full_file(collop_file_D),' (initial lines)'
          STOP
       END IF
    END DO
    ! read upper sizes
    READ(collop_unit,*,iostat=readstatus) collop_D_m_max
    IF (readstatus .NE. 0) THEN
       PRINT *, 'Can not read file: ',full_file(collop_file_D),' (size)'
       STOP
    END IF
    READ(collop_unit,*,iostat=readstatus) collop_D_mp_max
    IF (readstatus .NE. 0) THEN
       PRINT *, 'Can not read file: ',full_file(collop_file_D),' (size)'
       STOP
    END IF
    ! allocate
    ALLOCATE(collop_D_gee(collop_D_m_min:collop_D_m_max,collop_D_mp_min:collop_D_mp_max))
    ! read data
    DO m = collop_D_m_min,collop_D_m_max
       READ(collop_unit,*,iostat=readstatus) collop_D_gee(m,:)
       IF (readstatus .NE. 0) THEN
          PRINT *, 'Can not read file: ',full_file(collop_file_D),' at m=',m
          STOP
       END IF
    END DO
    ! close
    CLOSE(unit=collop_unit)
    ! -------------------------------------------------------------------------
    ! end D
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! I(m,mp,lp)
    !
    ! deallocate
    IF (ALLOCATED(collop_I_gee))  DEALLOCATE(collop_I_gee)
    ! open
    OPEN(unit=collop_unit,file=full_file(collop_file_I), & 
         status='old',action='read',iostat=openstatus)
    IF (openstatus .NE. 0) THEN
       PRINT *, 'Can not open file: ',full_file(collop_file_I)
       STOP
    END IF
    ! read initial lines
    DO k = 1, dummy_read
       READ(collop_unit,*,iostat=readstatus) dummy
       IF (readstatus .NE. 0) THEN
          PRINT *, 'Can not read file: ',full_file(collop_file_I),' (initial lines)'
          STOP
       END IF
    END DO
    ! read upper sizes
    READ(collop_unit,*,iostat=readstatus) collop_I_m_max
    IF (readstatus .NE. 0) THEN
       PRINT *, 'Can not read file: ',full_file(collop_file_I),' (size)'
       STOP
    END IF
    READ(collop_unit,*,iostat=readstatus) collop_I_mp_max
    IF (readstatus .NE. 0) THEN
       PRINT *, 'Can not read file: ',full_file(collop_file_I),' (size)'
       STOP
    END IF
    READ(collop_unit,*,iostat=readstatus) collop_I_lp_max
    IF (readstatus .NE. 0) THEN
       PRINT *, 'Can not read file: ',full_file(collop_file_I),' (size)'
       STOP
    END IF
    ! allocate
    ALLOCATE(collop_I_gee(collop_I_m_min:collop_I_m_max, &
         collop_I_mp_min:collop_I_mp_max, &
         collop_I_lp_min:collop_I_lp_max ))
    ! read data
    DO ll = collop_I_lp_min,collop_I_lp_max
       READ(collop_unit,*,iostat=readstatus) lp
       IF (readstatus .NE. 0) THEN
          PRINT *, 'Can not read file: ',full_file(collop_file_I),' at lp=',lp
          STOP
       END IF
       DO m = collop_I_m_min,collop_I_m_max
          READ(collop_unit,*,iostat=readstatus) collop_I_gee(m,:,lp)
          IF (readstatus .NE. 0) THEN
             PRINT *, 'Can not read file: ',full_file(collop_file_I),' at lp=',lp,' m=',m
             STOP
          END IF
       END DO
    END DO
    ! close
    CLOSE(unit=collop_unit)
    ! -------------------------------------------------------------------------
    ! end I
    ! -------------------------------------------------------------------------

  END SUBROUTINE collop_construct_all

  ! ---------------------------------------------------------------------------
  SUBROUTINE collop_deconstruct_all    
    ! -------------------------------------------------------------------------
    ! A(m,k) k = 1,2,3
    IF (ALLOCATED(collop_A))  DEALLOCATE(collop_A)
    IF (ALLOCATED(collop_A1)) DEALLOCATE(collop_A1)
    IF (ALLOCATED(collop_A2)) DEALLOCATE(collop_A2)
    IF (ALLOCATED(collop_A3)) DEALLOCATE(collop_A3)
    IF (ALLOCATED(collop_C))  DEALLOCATE(collop_C)
    collop_A_m_max = 0
    collop_C_m_max = 0
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! nu(m,mp)
    IF (ALLOCATED(collop_nu_gee))  DEALLOCATE(collop_nu_gee)
    IF (ALLOCATED(collop_nu_ginf))  DEALLOCATE(collop_nu_ginf)
    collop_nu_m_max  = 0
    collop_nu_mp_max = 0
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! D(m,mp)
    IF (ALLOCATED(collop_D_gee))  DEALLOCATE(collop_D_gee)
    collop_D_m_max  = 0
    collop_D_mp_max = 0
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! I(m,mp,lp)
    IF (ALLOCATED(collop_I_gee))  DEALLOCATE(collop_I_gee)
    collop_I_m_max  = 0
    collop_I_mp_max = 0
    collop_I_lp_max = 0
    ! -------------------------------------------------------------------------

  END SUBROUTINE collop_deconstruct_all

  ! ---------------------------------------------------------------------------
  ! private routines
  ! ---------------------------------------------------------------------------
  FUNCTION full_file_1(name_in) RESULT(name_out)
    CHARACTER(len=100), INTENT(in)  :: name_in
    CHARACTER(len=200) :: name_out

    WRITE(name_out,'(200A)') TRIM(ADJUSTL(collop_path)),TRIM(ADJUSTL(name_in))

  END FUNCTION full_file_1

  ! ---------------------------------------------------------------------------
  SUBROUTINE check_unit_1
    LOGICAL :: opened
    DO
       INQUIRE(unit=collop_unit,opened=opened)
       IF(.NOT. opened) EXIT
       collop_unit = collop_unit + 1
    END DO
  END SUBROUTINE check_unit_1
!
  SUBROUTINE collop_load_all
    !
    USE rkstep_mod, ONLY : lag,leg,asource,anumm,denmm,ailmm,weightlag
    !
integer :: lag1,lag2
    REAL(kind=dp), PARAMETER :: pi=3.14159265358979d0
    ! REAL(kind=dp)            :: z_eff ! moved to module
!
    IF(collop_A_m_min.NE.0    .OR.                                          &
         collop_nu_m_min.NE.0   .OR.                                          &
         collop_nu_mp_min.NE.0  .OR.                                          &
         collop_D_m_min.NE.0    .OR.                                          &
         collop_D_mp_min.NE.0   .OR.                                          &
         collop_I_m_min.NE.0    .OR.                                          &
         collop_I_mp_min.NE.0   .OR.                                          &
         collop_I_lp_min.NE.0) THEN
       PRINT *,'starting collop index is not zero:',                           &
            collop_A_m_min,collop_nu_m_min,collop_nu_mp_min,collop_D_m_min,       &
            collop_D_mp_min,collop_I_m_min,collop_I_mp_min,collop_I_lp_min
       STOP
    ENDIF
!
    IF(collop_A_m_max.LT.lag  .OR.                                          &
         collop_nu_m_max.LT.lag .OR.                                          &
         collop_nu_mp_max.LT.lag.OR.                                          &
         collop_D_m_max.LT.lag  .OR.                                          &
         collop_D_mp_max.LT.lag .OR.                                          &
         collop_I_m_max.LT.lag  .OR.                                          &
         collop_I_mp_max.LT.lag .OR.                                          &
         collop_I_lp_max.LT.leg) THEN
       PRINT *,'collop index is out of range',lag,'>',collop_A_m_max,        &
            collop_nu_m_max,collop_nu_mp_max,                            &
            collop_D_m_max,collop_D_mp_max,                              &
            collop_I_m_max,collop_I_mp_max,                              &
            ' or ',leg,' > ',collop_I_lp_max
       STOP
    ENDIF
!
    IF(ALLOCATED(anumm)) DEALLOCATE(anumm)
    ALLOCATE(anumm(0:lag,0:lag))
    IF(ALLOCATED(denmm)) DEALLOCATE(denmm)
    ALLOCATE(denmm(0:lag,0:lag))
    IF(ALLOCATED(asource)) DEALLOCATE(asource)
    ALLOCATE(asource(0:lag,3))
    IF(ALLOCATED(ailmm)) DEALLOCATE(ailmm)
    ALLOCATE(ailmm(0:lag,0:lag,0:leg))
    IF(ALLOCATED(weightlag)) DEALLOCATE(weightlag)
    ALLOCATE(weightlag(3,0:lag))
!
! Comment from 19.09.07:
! Here sources, asource, and flux convolution factors, weightlag, have 
! opposite signs compared to the short writeup of Georg. This does not change
! signs of diffusion matrix
!
    asource(0:lag,1)=collop_A(0:lag,1)
    asource(0:lag,2)=collop_A(0:lag,3)
    asource(0:lag,3)=collop_A(0:lag,2)
    anumm(0:lag,0:lag)=collop_nu_gee(0:lag,0:lag)                           &
         +z_eff*collop_nu_ginf(0:lag,0:lag)
    denmm(0:lag,0:lag)=collop_D_gee(0:lag,0:lag)
    ailmm(0:lag,0:lag,0:leg)=collop_I_gee(0:lag,0:lag,0:leg)
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
!Truncation:
!do lag1=0,lag
!  do lag2=0,lag
!     if(lag1+lag2.gt.lag) then
!       anumm(lag1,lag2)=0.d0
!       denmm(lag1,lag2)=0.d0
!       ailmm(lag1,lag2,0:leg)=0.d0
!     endif
!  enddo
!enddo
!End truncation
    !
  END SUBROUTINE collop_load_all
!
  SUBROUTINE collop_unload_all
    !
    USE rkstep_mod, ONLY : asource,anumm,denmm,ailmm,weightlag
    !
    IF(ALLOCATED(anumm)) DEALLOCATE(anumm)
    IF(ALLOCATED(denmm)) DEALLOCATE(denmm)
    IF(ALLOCATED(asource)) DEALLOCATE(asource)
    IF(ALLOCATED(ailmm)) DEALLOCATE(ailmm)
    IF(ALLOCATED(weightlag)) DEALLOCATE(weightlag)
  END SUBROUTINE collop_unload_all

END MODULE collop
