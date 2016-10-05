!
  MODULE compute_aiota_mod
!
    INTEGER :: nmax=0
    DOUBLE PRECISION :: dz_dphi
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DYDX,YT,DYT,DYM
!
  CONTAINS
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      SUBROUTINE RK4D(Y,N,X,H,DERIVS)
!
!
      DOUBLE PRECISION Y,X,H,HH,H6,XH
      EXTERNAL DERIVS
      DIMENSION Y(N)
!
      IF(n.NE.nmax) THEN
        IF(ALLOCATED(dydx)) DEALLOCATE(dydx)
        IF(ALLOCATED(yt)) DEALLOCATE(yt)
        IF(ALLOCATED(dyt)) DEALLOCATE(dyt)
        IF(ALLOCATED(dym)) DEALLOCATE(dym)
        nmax=n
        ALLOCATE(DYDX(NMAX),YT(NMAX),DYT(NMAX),DYM(NMAX))
      ENDIF
!
      HH=H*0.5D0
      H6=H/6.D0
      XH=X+HH
      CALL DERIVS(N,X,Y,DYDX)
      DO 11 I=1,N
        YT(I)=Y(I)+HH*DYDX(I)
11    CONTINUE
      CALL DERIVS(N,XH,YT,DYT)
      DO 12 I=1,N
        YT(I)=Y(I)+HH*DYT(I)
12    CONTINUE
      CALL DERIVS(N,XH,YT,DYM)
      DO 13 I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
13    CONTINUE
      CALL DERIVS(N,X+H,YT,DYT)
      DO 14 I=1,N
        Y(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.D0*DYM(I))
14    CONTINUE
        X=X+H
      RETURN
      END SUBROUTINE RK4D  
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  SUBROUTINE rhs1(ndim,phi,y,dery)
!
  IMPLICIT NONE
!
  INTEGER :: ndim ! = 5
!
  DOUBLE PRECISION :: phi,y,dery
  DOUBLE PRECISION x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder
  DIMENSION y(ndim),dery(ndim)
  DIMENSION x(3),bder(3),hcovar(3),hctrvr(3),hcoder(3,3),hctder(3,3)
!
  x(1)=y(1)
  x(2)=phi
  x(3)=y(2)
!
  CALL mag_efit(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
!
  dery(1)=hctrvr(1)/hctrvr(2)
  dery(2)=hctrvr(3)/hctrvr(2)
  dery(3)=y(1)*hctrvr(3)/hctrvr(2)
  IF(ndim.EQ.5) THEN
    dery(4)=y(1)
    dery(5)=y(2)
  ELSEIF(ndim.EQ.4) THEN
    dery(4)=bmod*y(1)*y(2)*hctrvr(1)
  ENDIF
!
  dz_dphi=dery(2)
!
  RETURN
  END SUBROUTINE rhs1
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE mag_efit(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
!
! Computes magnetic field module normalized to axial value  - bmod,
! square root of determinant of the metric tensor           - sqrtg,
! derivatives of the logarythm of the magnetic field module
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
!  Output parameters:
!            formal:  bmod
!                     sqrtg
!                     bder
!                     hcovar
!                     hctrvr
!                     hcoder
!
!  Called routines:  field_eq
!
      DOUBLE PRECISION x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder
      DOUBLE PRECISION hr,hf,hz
!
      DOUBLE PRECISION ri,fii,zi,br,bf,bz, &
      BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ, &
      BRK,BZK,BRRK,BRZK,BZRK,BZZK
!
      !! Modification by Andreas F. Martitsch (16.07.2015)
      ! fixed Warning: Possible change of value in conversion from
      ! REAL(8) to REAL(4) at (1)
      DOUBLE PRECISION rbig
      !! End Modification by Andreas F. Martitsch (16.07.2015)
!
      DIMENSION x(3),bder(3),hcovar(3),hctrvr(3),hcoder(3,3),hctder(3,3)
!
      rbig=MAX(x(1),1d-12)
!
!cccccc computation of gb in cylindrical co-ordinates cccccccc
      ri=rbig
      fii=x(2)
      zi=x(3)

      CALL field_eq(ri,fii,zi,br,bf,bz, &
      BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ)

!ccccc end of gb computation cccccccccc
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
      END SUBROUTINE mag_efit
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  SUBROUTINE compute_aiota(R,raxis,zaxis,aiota,ierr)
!
  USE field_eq_mod, ONLY : icall_eq,rtf,btf,nrad,nzet,rad,zet             &
                         , psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
  USE field_c_mod,  ONLY : icall_c
!
  IMPLICIT NONE
!
  DOUBLE PRECISION, PARAMETER :: pi=3.14159265358979d0
!
  INTEGER :: ndim,ndimc,nstep,nmap,ntotstep,niter,iter,i,j,ind,m,n,i1
  INTEGER :: nsurf,isurf,nsurfmax,ndim_fc,ntheta,nsqpsi,nsubstep
  INTEGER :: k1,k2,k3,k4,numbig,npoisep,nlabel,k,ierr
  DOUBLE PRECISION :: R,aiota
  DOUBLE PRECISION :: h,phi,phibeg,rbeg,zbeg,raxis,zaxis,hr,sig,eps_four
  DOUBLE PRECISION :: hh,hh0,sig0,rmn,rmx,zmn,zmx,psi_axis,dphidpsi
  DOUBLE PRECISION :: rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ    &
                      ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: y
!
! Initialization of the field:
!
  rrr=1.d0
  ppp=0.d0
  zzz=0.d0
!
  CALL field(rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
            ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
! End of initialization
!
! Computation box:
  rmn=rad(1)
  rmx=rad(nrad)
  zmn=zet(1)
  zmx=zet(nzet)
!
  rbeg=0.5d0*(rmn+rmx)
  zbeg=0.5d0*(zmn+zmx)
!
  rrr=rbeg
  ppp=0.d0
  zzz=zbeg
!
  nstep=360
  nmap=10
  niter=10
  ntotstep=nstep*nmap
  h=2.d0*pi/nstep
!
! Search for the magnetic axis
!
  ndim=5
  ALLOCATE(y(ndim))
  phi=0.d0
  y=0.d0
  y(1)=rbeg
  y(2)=zbeg
  DO iter=1,niter
    phibeg=phi
    y(4:5)=0.d0
    DO i=1,ntotstep
      CALL RK4D(y,ndim,phi,h,rhs1)
    ENDDO
    y(1:2)=y(4:5)/(phi-phibeg)
  ENDDO
  raxis=y(1)
  zaxis=y(2)
!
! Computation of iota, flux label, surface area
!
  ierr=0
  phi=0.d0
  y(1)=R
  y(2)=zaxis
  y(3)=0.d0
  y(4)=0.d0
  CALL RK4D(y,ndim,phi,h,rhs1)
  sig=y(2)-zaxis
  DO WHILE(sig*(y(2)-zaxis).GT.0.d0)
    CALL RK4D(y,ndim,phi,h,rhs1)
    IF( y(1).LT.rmn .OR. y(1).GT.rmx .OR.            &
        y(2).LT.zmn .OR. y(2).GT.zmx ) THEN
      ierr=1
      RETURN
    ENDIF
  ENDDO
  sig=y(2)-zaxis
  DO WHILE(sig*(y(2)-zaxis).GT.0.d0)
    CALL RK4D(y,ndim,phi,h,rhs1)
  ENDDO
! Newton method
  DO iter=1,niter
    hh=(zaxis-y(2))/dz_dphi
    sig0=sig
    CALL RK4D(y,ndim,phi,hh,rhs1)
    IF( y(1).LT.rmn .OR. y(1).GT.rmx .OR.            &
        y(2).LT.zmn .OR. y(2).GT.zmx ) THEN
      ierr=1
      RETURN
    ENDIF
  ENDDO
!
  aiota=2.d0*pi/phi
!
  END SUBROUTINE compute_aiota
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  END MODULE compute_aiota_mod
