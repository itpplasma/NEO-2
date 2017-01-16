MODULE mag_sub
  !! Modifications by Andreas F. Martitsch (09.03.2014)
  ! Collection of subroutines converted to a module.
  ! This allows one to make use of generic interfaces 
  ! for module procedures (e.g., for a different number of
  ! arguments).
  ! Note: This requires changes in "flint_prepare" and 
  ! "write_volume_data" (both flint.f90), and "rhs_kin"
  ! and "magdata_for_particles". 
  !! End Modifications by Andreas F. Martitsch (09.03.2014)

  IMPLICIT NONE

  PUBLIC mag
  PRIVATE mag_a, mag_a1
  INTERFACE mag
     MODULE PROCEDURE mag_a, mag_a1
  END INTERFACE mag

  PRIVATE mag_real
  PRIVATE mag_real_a
  INTERFACE mag_real
     MODULE PROCEDURE mag_real_a
  END INTERFACE mag_real

  PRIVATE mag_homogeneous
  PRIVATE mag_homogeneous_a
  INTERFACE mag_homogeneous
     MODULE PROCEDURE mag_homogeneous_a
  END INTERFACE mag_homogeneous

  PRIVATE mag_legendre
  PRIVATE mag_legendre_a
  INTERFACE mag_legendre
     MODULE PROCEDURE mag_legendre_a
  END INTERFACE mag_legendre

  ! COMMON-block is depricated (not necessary within a module)
  ! (It seems that these variables do not interact with 
  ! external routines --> Set to private)
  ! (Furthermore asprat and qsaf0 were never defined 
  ! explicitly -> No IMPLICIT NONE in the code)
  private rbig0, asprat, qsaf0
  double precision rbig0, asprat, qsaf0

CONTAINS
  ! ------------------------------------------------------------------------
  SUBROUTINE mag_a(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
    USE mag_interface_mod, ONLY : mag_magfield,mag_coordinates
    USE magfie_mod, ONLY : magfie
    USE compute_aiota_mod, ONLY: mag_efit
    DOUBLE PRECISION x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder
    DOUBLE PRECISION hcurl
    DIMENSION x(3),bder(3),hcovar(3),hctrvr(3),hcoder(3,3),hctder(3,3)
    DIMENSION hcurl(3)

    IF (mag_coordinates .EQ. 0) THEN
       ! cylindrical
       IF (mag_magfield .EQ. 0) THEN ! homogeneous
          CALL mag_homogeneous(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
       ELSEIF (mag_magfield .EQ. 1) THEN ! W7-AS, Tokamak
          CALL mag_real(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
       ELSEIF (mag_magfield .EQ. 2) THEN ! Legendre
          CALL mag_legendre(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
       ELSEIF (mag_magfield .EQ. 3) THEN ! EFIT
          CALL mag_efit(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
       ELSE
          PRINT *, 'Not implemented: mag_magfield = ',mag_magfield
          STOP
       END IF
!!$     IF (mag_magfield .NE. 0) THEN
!!$        CALL mag_real(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
!!$     ELSE
!!$        CALL mag_homogeneous(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
!!$     END IF
    ELSE
       ! Boozer
       CALL magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
       hcoder = 0.0d0
       hctder = 0.0d0
    END IF
    !PRINT *, x,bmod
  END SUBROUTINE mag_a
  ! ------------------------------------------------------------------------
  !! Modifications by Andreas F. Martitsch (09.03.2014)
  ! Optional output (necessary for modeling the magnetic rotation)
  SUBROUTINE mag_a1(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder,bcovar_s_hat_der)
    USE mag_interface_mod, ONLY : mag_magfield,mag_coordinates
    USE magfie_mod, ONLY : magfie
    USE compute_aiota_mod, ONLY: mag_efit
    DOUBLE PRECISION x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder
    DOUBLE PRECISION hcurl,bcovar_s_hat_der
    DIMENSION x(3),bder(3),hcovar(3),hctrvr(3),hcoder(3,3),hctder(3,3)
    DIMENSION hcurl(3),bcovar_s_hat_der(3)

    IF (mag_coordinates .EQ. 0) THEN
       ! cylindrical
       IF (mag_magfield .EQ. 0) THEN ! homogeneous
          CALL mag_homogeneous(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
       ELSEIF (mag_magfield .EQ. 1) THEN ! W7-AS, Tokamak
          CALL mag_real(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
       ELSEIF (mag_magfield .EQ. 2) THEN ! Legendre
          CALL mag_legendre(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
       ELSEIF (mag_magfield .EQ. 3) THEN ! EFIT
          CALL mag_efit(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
       ELSE
          PRINT *, 'Not implemented: mag_magfield = ',mag_magfield
          STOP
       END IF
       bcovar_s_hat_der = 0.0d0 ! Not needed in this branch (hcoder-array is available)
!!$     IF (mag_magfield .NE. 0) THEN
!!$        CALL mag_real(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
!!$     ELSE
!!$        CALL mag_homogeneous(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
!!$     END IF
    ELSE
       ! Boozer
       CALL magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl,bcovar_s_hat_der)
       hcoder = 0.0d0
       hctder = 0.0d0
    END IF
    !PRINT *, x,bmod
  END SUBROUTINE mag_a1
  !! End Modifications by Andreas F. Martitsch (09.03.2014)
  ! ------------------------------------------------------------------------
  SUBROUTINE mag_real_a(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
    !
    ! Computes magnetic field MODULE normalized to axial value  - bmod,
    ! square root of determinant of the metric tensor           - sqrtg,
    ! derivatives of the logarythm of the magnetic field MODULE
    ! over coordinates                                          - bder,
    ! covariant componets of the unit vector of the magnetic
    ! field direction                                           - hcovar,
    ! contravariant components of this vector                   - hctrvr,
    ! derivatives of the covariant components of this vector    - hcoder(i,j)
    ! here hcoder(i,j) is the derivative of hcovar(j) over x(i)
    ! for given set of coordinates x(i).
    ! Oder of coordinates is the following: x(1)=R (big radius), 
    ! x(2)=phi (toroidal angle), x(3)=z (altitude).
    !
    !  Input parameters:
    !            formal:  x                -    array of coordinates
    !            COMMON:  rbig0            -    big radius of mag. axis
    !  Output parameters:
    !            formal:  bmod
    !                     sqrtg
    !                     bder
    !                     hcovar
    !                     hctrvr
    !                     hcoder
    !
    !  Called routines:  GBhs,GBRZd 
    !
    USE mag_interface_mod, ONLY : magnetic_device
    !
    DOUBLE PRECISION x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder
    DOUBLE PRECISION rbig0,hr,hf,hz
    !
    DOUBLE PRECISION ri,fii,zi,br,bf,bz
    DOUBLE PRECISION BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ
    DOUBLE PRECISION BRK,BZK,BRRK,BRZK,BZRK,BZZK
    !
    double precision rbig
    !
    DIMENSION x(3),bder(3),hcovar(3),hctrvr(3),hcoder(3,3),hctder(3,3)
    !
    !! Modifications by Andreas F. Martitsch (09.03.2014)
    ! COMMON-block is depricated (not necessary within a module)
    ! COMMON/magpar/rbig0,asprat,qsaf0
    !! End Modifications by Andreas F. Martitsch (09.03.2014)
    !
    rbig=MAX(x(1),1d-12)
    !
    ! computation of gb in cylindrical co-ordinates cccccccc
    ri=rbig
    fii=x(2)
    zi=x(3)

    IF (magnetic_device .EQ. 0) THEN ! Tokamak
       CALL GBas_0(ri,fii,zi,br,bf,bz,  &
            BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ)
    ELSEIF (magnetic_device .EQ. 1) THEN ! W7-AS
       CALL GBas_1(ri,fii,zi,br,bf,bz,  &
            BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ)
    ELSEIF (magnetic_device .EQ. 2) THEN ! Legendre        !<-in
       CALL GBhs_l(ri,fii,zi,br,bf,bz,  &                  !<-in
            BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ)           !<-in
       CALL GBRZd_l(ri,zi,BRK,BZK,BRRK,BRZK,BZRK,BZZK)     !<-in
       BR=BR+BRK                                           !<-in
       BZ=BZ+BZK                                           !<-in
       BRR=BRR+BRRK                                        !<-in
       BRZ=BRZ+BRZK                                        !<-in
       BZR=BZR+BZRK                                        !<-in
       BZZ=BZZ+BZZK                                        !<-in
    ELSE
       PRINT *, 'Magnetic Device not implemented'
       STOP
    END IF

    ! END of gb computation cccccccccc
    bmod=dsqrt(br**2+bf**2+bz**2)
    sqrtg=rbig
    hr=br/bmod
    hf=bf/bmod
    hz=bz/bmod
    !
    bder(1)=(brr*hr+bfr*hf+bzr*hz)/bmod
    bder(2)=(brf*hr+bff*hf+bzf*hz)/bmod
    bder(3)=(brz*hr+bfz*hf+bzz*hz)/bmod
    !
    hcovar(1)=hr
    hcovar(2)=hf*rbig
    hcovar(3)=hz
    !
    hctrvr(1)=hr
    hctrvr(2)=hf/rbig
    hctrvr(3)=hz
    !
    hcoder(1,1)=brr/bmod-hcovar(1)*bder(1)
    hcoder(2,1)=brf/bmod-hcovar(1)*bder(2)
    hcoder(3,1)=brz/bmod-hcovar(1)*bder(3)
    hcoder(1,2)=(rbig*bfr+bf)/bmod-hcovar(2)*bder(1)
    hcoder(2,2)=rbig*bff/bmod-hcovar(2)*bder(2)
    hcoder(3,2)=rbig*bfz/bmod-hcovar(2)*bder(3)
    hcoder(1,3)=bzr/bmod-hcovar(3)*bder(1)
    hcoder(2,3)=bzf/bmod-hcovar(3)*bder(2)
    hcoder(3,3)=bzz/bmod-hcovar(3)*bder(3)
    !
    hctder(1,1)=brr/bmod-hctrvr(1)*bder(1)
    hctder(2,1)=brf/bmod-hctrvr(1)*bder(2)
    hctder(3,1)=brz/bmod-hctrvr(1)*bder(3)
    hctder(1,2)=(bfr-bf/rbig)/(rbig*bmod)-hctrvr(2)*bder(1)
    hctder(2,2)=bff/(rbig*bmod)-hctrvr(2)*bder(2)
    hctder(3,2)=bfz/(rbig*bmod)-hctrvr(2)*bder(3)
    hctder(1,3)=bzr/bmod-hctrvr(3)*bder(1)
    hctder(2,3)=bzf/bmod-hctrvr(3)*bder(2)
    hctder(3,3)=bzz/bmod-hctrvr(3)*bder(3)
    !
    RETURN
  END SUBROUTINE mag_real_a
  !
  ! ------------------------------------------------------------------------
!!$SUBROUTINE mag_efit(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
!!$  !
!!$  ! This is the EFIT-Version of mag_real
!!$  !
!!$  ! Computes magnetic field MODULE normalized to axial value  - bmod,
!!$  ! square root of determinant of the metric tensor           - sqrtg,
!!$  ! derivatives of the logarythm of the magnetic field MODULE
!!$  ! over coordinates                                          - bder,
!!$  ! covariant componets of the unit vector of the magnetic
!!$  ! field direction                                           - hcovar,
!!$  ! contravariant components of this vector                   - hctrvr,
!!$  ! derivatives of the covariant components of this vector    - hcoder(i,j)
!!$  ! here hcoder(i,j) is the derivative of hcovar(j) over x(i)
!!$  ! for given set of coordinates x(i).
!!$  ! Oder of coordinates is the following: x(1)=R (big radius), 
!!$  ! x(2)=phi (toroidal angle), x(3)=z (altitude).
!!$  !
!!$  !  Input parameters:
!!$  !            formal:  x                -    array of coordinates
!!$  !  Output parameters:
!!$  !            formal:  bmod
!!$  !                     sqrtg
!!$  !                     bder
!!$  !                     hcovar
!!$  !                     hctrvr
!!$  !                     hcoder
!!$  !
!!$  !  Called routines:  field 
!!$  !
!!$  !USE mag_interface_mod, ONLY : magnetic_device
!!$  !
!!$  !implicit none
!!$
!!$  DOUBLE PRECISION x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder
!!$  DOUBLE PRECISION rbig0
!!$  DOUBLE PRECISION hr,hf,hz
!!$  !
!!$  DOUBLE PRECISION ri,fii,zi,br,bf,bz
!!$  DOUBLE PRECISION BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ
!!$  DOUBLE PRECISION BRK,BZK,BRRK,BRZK,BZRK,BZZK
!!$  !
!!$  DIMENSION x(3),bder(3),hcovar(3),hctrvr(3),hcoder(3,3),hctder(3,3)
!!$  !
!!$  !COMMON/magpar/rbig0,asprat,qsaf0
!!$  !
!!$  rbig=MAX(x(1),1d-12)
!!$  !
!!$  ! computation of gb in cylindrical co-ordinates cccccccc
!!$  ri=rbig
!!$  fii=x(2)
!!$  zi=x(3)
!!$
!!$  !call field(R,phi,Z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,   &
!!$  !     dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!!$  call field(ri,fii,zi,br,bf,bz, &
!!$       BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ)
!!$
!!$  ! END of gb computation cccccccccc
!!$  bmod=dsqrt(br**2+bf**2+bz**2)
!!$  sqrtg=rbig
!!$  hr=br/bmod
!!$  hf=bf/bmod
!!$  hz=bz/bmod
!!$  !
!!$  bder(1)=(brr*hr+bfr*hf+bzr*hz)/bmod
!!$  bder(2)=(brf*hr+bff*hf+bzf*hz)/bmod
!!$  bder(3)=(brz*hr+bfz*hf+bzz*hz)/bmod
!!$  !
!!$  hcovar(1)=hr
!!$  hcovar(2)=hf*rbig
!!$  hcovar(3)=hz
!!$  !
!!$  hctrvr(1)=hr
!!$  hctrvr(2)=hf/rbig
!!$  hctrvr(3)=hz
!!$  !
!!$  hcoder(1,1)=brr/bmod-hcovar(1)*bder(1)
!!$  hcoder(2,1)=brf/bmod-hcovar(1)*bder(2)
!!$  hcoder(3,1)=brz/bmod-hcovar(1)*bder(3)
!!$  hcoder(1,2)=(rbig*bfr+bf)/bmod-hcovar(2)*bder(1)
!!$  hcoder(2,2)=rbig*bff/bmod-hcovar(2)*bder(2)
!!$  hcoder(3,2)=rbig*bfz/bmod-hcovar(2)*bder(3)
!!$  hcoder(1,3)=bzr/bmod-hcovar(3)*bder(1)
!!$  hcoder(2,3)=bzf/bmod-hcovar(3)*bder(2)
!!$  hcoder(3,3)=bzz/bmod-hcovar(3)*bder(3)
!!$  !
!!$  hctder(1,1)=brr/bmod-hctrvr(1)*bder(1)
!!$  hctder(2,1)=brf/bmod-hctrvr(1)*bder(2)
!!$  hctder(3,1)=brz/bmod-hctrvr(1)*bder(3)
!!$  hctder(1,2)=(bfr-bf/rbig)/(rbig*bmod)-hctrvr(2)*bder(1)
!!$  hctder(2,2)=bff/(rbig*bmod)-hctrvr(2)*bder(2)
!!$  hctder(3,2)=bfz/(rbig*bmod)-hctrvr(2)*bder(3)
!!$  hctder(1,3)=bzr/bmod-hctrvr(3)*bder(1)
!!$  hctder(2,3)=bzf/bmod-hctrvr(3)*bder(2)
!!$  hctder(3,3)=bzz/bmod-hctrvr(3)*bder(3)
!!$  !
!!$  RETURN
!!$END SUBROUTINE mag_efit
!!$!
!!$! ------------------------------------------------------------------------
  !
  SUBROUTINE mag_homogeneous_a(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
    ! same as above for the homogeneous case
    DOUBLE PRECISION x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder
    DOUBLE PRECISION rbig0,hr,hf,hz
    !
    DOUBLE PRECISION ri,fii,zi,br,bf,bz
    DOUBLE PRECISION BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ
    DOUBLE PRECISION BRK,BZK,BRRK,BRZK,BZRK,BZZK
    !
    DIMENSION x(3),bder(3),hcovar(3),hctrvr(3),hcoder(3,3),hctder(3,3)
    !
    !! Modifications by Andreas F. Martitsch (09.03.2014)
    ! COMMON-block is depricated (not necessary within a module)
    ! COMMON/magpar/rbig0,asprat,qsaf0
    !! End Modifications by Andreas F. Martitsch (09.03.2014)
    !
    rbig0=200.d0
    !
    ! computation of gb in cylindrical co-ordinates cccccccc
    ri=rbig0
    fii=x(2)
    zi=x(3)

    ! end of gb computation
    bmod=1.d0
    sqrtg=rbig0
    hr=0.d0
    hf=1.d0
    hz=0.d0
    !
    bder(1)=0.d0
    bder(2)=0.d0
    bder(3)=0.d0
    !
    hcovar(1)=0.d0
    hcovar(2)=hf*rbig0
    hcovar(3)=0.d0
    !
    hctrvr(1)=0.d0
    hctrvr(2)=hf/rbig0
    hctrvr(3)=0.d0
    !
    hcoder=0.d0
    !
    hctder=0.d0
    !
    RETURN
  END SUBROUTINE mag_homogeneous_a
  ! ------------------------------------------------------------------------
  SUBROUTINE mag_legendre_a(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
    !
    ! Computes magnetic field module in units of the magnetic code  - bmod,
    !square root of determinant of the metric tensor               - sqrtg,
    !derivatives of the logarythm of the magnetic field MODULE
    ! over coordinates                                              - bder,
    ! covariant componets of the unit vector of the magnetic
    ! field direction                                               - hcovar,
    ! contravariant components of this vector                       - hctrvr,
    !contravariant component of the curl of this vector            - hcurl
    ! Order of coordinates is the following: x(1)=R (big radius), 
    ! x(2)=phi (toroidal angle), x(3)=Z (altitude).
    !
    ! Input parameters:
    ! formal:  x                -    array of coordinates
    !  Output parameters:
    !            formal:  bmod
    !                     sqrtg
    !                     bder
    !                     hcovar
    !                     hctrvr
    !                     hcurl
    !
    !  Called routines:  GBhs_l,GBRZd_l 
    !
    DOUBLE PRECISION, PARAMETER :: phi_shift=3.14159265358979d0/6.d0
    DOUBLE PRECISION x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
    DOUBLE PRECISION hr,hf,hz
    !
    DOUBLE PRECISION ri,fii,zi,br,bf,bz
    DOUBLE PRECISION BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ
    DOUBLE PRECISION BRK,BZK,BRRK,BRZK,BZRK,BZZK
    !
    double precision rbig
    !
    DIMENSION x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
    !
    rbig=MAX(x(1),1d-12)
    !
    ! computation of gb in cylindrical co-ordinates cccccccc
    ri=rbig
    fii=x(2)+phi_shift
    zi=x(3)

    CALL GBhs_l(ri,fii,zi,br,bf,bz, &
         BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ)

    CALL GBRZd_l(ri,zi,BRK,BZK,BRRK,BRZK,BZRK,BZZK)
    br=br+BRK
    bz=bz+BZK
    BRR=BRR+BRRK
    BRZ=BRZ+BRZK
    BZR=BZR+BZRK
    BZZ=BZZ+BZZK
    ! end of gb computation cccccccccc
    bmod=dsqrt(br**2+bf**2+bz**2)
    sqrtg=rbig
    hr=br/bmod
    hf=bf/bmod
    hz=bz/bmod
    !
    bder(1)=(brr*hr+bfr*hf+bzr*hz)/bmod
    bder(2)=(brf*hr+bff*hf+bzf*hz)/bmod
    bder(3)=(brz*hr+bfz*hf+bzz*hz)/bmod
    !
    hcovar(1)=hr
    hcovar(2)=hf*rbig
    hcovar(3)=hz
    !
    hctrvr(1)=hr
    hctrvr(2)=hf/rbig
    hctrvr(3)=hz
    !
    hcurl(1)=((bzf-rbig*bfz)/bmod+ &
         hcovar(2)*bder(3)-hcovar(3)*bder(2))/sqrtg
    hcurl(2)=((brz-bzr)/bmod+ &
         hcovar(3)*bder(1)-hcovar(1)*bder(3))/sqrtg
    hcurl(3)=((bf+rbig*bfr-brf)/bmod+ &
         hcovar(1)*bder(2)-hcovar(2)*bder(1))/sqrtg
    !
    !PRINT *, x,bmod
    RETURN

  END SUBROUTINE mag_legendre_a
  ! ------------------------------------------------------------------------
END MODULE mag_sub
