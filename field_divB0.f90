MODULE period_mod 
  DOUBLE PRECISION :: per_phi, per_tht
END MODULE period_mod

MODULE input_files
  CHARACTER*1024 :: eqfile, cfile, gfile,pfile,convexfile,fluxdatapath
  INTEGER :: iunit=1738
  INTEGER :: ieqfile=1
!
  DATA eqfile  /'ASDEX/d3d-087506.03687.equ'/
  DATA cfile   /'DATA/ccoil.dat'/
!  data gfile   /'gfiles/shot115452/g115452.03525'/
!  data pfile   /'Conly/probe_g129_bfield.out'/  
END MODULE input_files
!
MODULE field_c_mod
  INTEGER :: icall_c=0
  INTEGER :: ntor=16
  INTEGER :: nr,np,nz
  INTEGER :: icftype
END MODULE field_c_mod
!
MODULE field_eq_mod
  INTEGER :: icall_eq=0
  INTEGER :: nrad,nzet,icp,nwindow_r,nwindow_z
  REAL(kind=8), PARAMETER                      :: pi=3.14159265358979d0
  REAL(kind=8) :: psib,btf,rtf,hrad,hzet
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE    :: psi, psi0
  REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE  :: splpsi
  REAL(kind=8), DIMENSION(:), ALLOCATABLE      :: rad, zet, xi,f
  INTEGER, DIMENSION(:), ALLOCATABLE           :: imi,ima,jmi,jma
  INTEGER, DIMENSION(:,:), ALLOCATABLE         :: ipoint
  DOUBLE PRECISION :: psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
END MODULE field_eq_mod
!
MODULE field_mod
  INTEGER          :: icall=0
  INTEGER          :: ipert,iequil
  DOUBLE PRECISION :: ampl
END MODULE field_mod
!
MODULE inthecore_mod
  LOGICAL :: prop=.TRUE.
  INTEGER :: npoi,ijumpb,ibeg,iend
  DOUBLE PRECISION, PARAMETER :: epssep=1.d-6
  DOUBLE PRECISION :: rc,zc,twopi,sig,psi_sep,psi_cut,sigpsi,cutoff
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rho2i,theti
  INTEGER          :: incore
  DOUBLE PRECISION :: vacf,dvacdr,dvacdz,d2vacdr2,d2vacdrdz,d2vacdz2
  DOUBLE PRECISION :: plaf,dpladr,dpladz,d2pladr2,d2pladrdz,d2pladz2
END MODULE inthecore_mod
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! -----------------------------------------------------------------
SUBROUTINE field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
                ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  USE period_mod
  USE input_files
  USE field_c_mod,   ONLY : ntor,icftype
  USE field_mod
  USE inthecore_mod, ONLY : incore,cutoff
  USE field_eq_mod, ONLY : nwindow_r,nwindow_z
!
  IMPLICIT NONE
!
  DOUBLE PRECISION :: r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
                     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  DOUBLE PRECISION :: rm,zm,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc   &
                     ,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc
!
  IF(icall .EQ. 0) THEN
     icall = 1
     OPEN(iunit, file='field_divB0.inp')
     READ(iunit,*) ipert        ! 0=eq only, 1=vac, 2=vac+plas no derivatives, 
                                ! 3=plas+vac with derivatives
     READ(iunit,*) iequil       ! 0=perturbation alone, 1=with equilibrium
     READ(iunit,*) ampl         ! amplitude of perturbation, a.u.
     READ(iunit,*) ntor         ! number of toroidal harmonics
     READ(iunit,*) cutoff       ! inner cutoff in psi/psi_a units
     READ(iunit,*) icftype      ! type of coil file
     READ(iunit,*) gfile        ! equilibrium file
     READ(iunit,*) pfile        ! coil        file
     READ(iunit,*) convexfile   ! convex file for stretchcoords
     READ(iunit,*) fluxdatapath ! directory with data in flux coord.
     READ(iunit,*) nwindow_r    ! widow size for filtering of psi array over R
     READ(iunit,*) nwindow_z    ! widow size for filtering of psi array over Z
     READ(iunit,*,err=1) ieqfile! equilibrium file type (0 - old, 1 - EFIT)
1    CLOSE(iunit)
     PRINT *, 'Perturbation field',ipert,'Ampl',ampl
  ENDIF

  CALL stretch_coords(r,z,rm,zm)

  CALL field_eq(rm,p,zm,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
               ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ) 
!
  IF(iequil.EQ.0) THEN
    Br=0.d0
    Bp=0.d0
    Bz=0.d0
    dBrdR=0.d0
    dBrdp=0.d0
    dBrdZ=0.d0
    dBpdR=0.d0
    dBpdp=0.d0
    dBpdZ=0.d0
    dBzdR=0.d0
    dBzdp=0.d0
    dBzdZ=0.d0
  ENDIF
!
  IF(ipert.GT.0) THEN
!
    IF(ipert.GT.1) THEN
      CALL inthecore(rm,zm)
    ELSE
      incore=-1
    ENDIF
!
! vacuum perturbation coil field:
!
    CALL field_c(rm,p,zm,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc   &
                ,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc) 
!
    Br = Br + Brc*ampl
    Bp = Bp + Bpc*ampl
    Bz = Bz + Bzc*ampl
    dBrdR = dBrdR + dBrdRc*ampl
    dBrdp = dBrdp + dBrdpc*ampl
    dBrdZ = dBrdZ + dBrdZc*ampl
    dBpdR = dBpdR + dBpdRc*ampl
    dBpdp = dBpdp + dBpdpc*ampl
    dBpdZ = dBpdZ + dBpdZc*ampl
    dBzdR = dBzdR + dBzdRc*ampl
    dBzdp = dBzdp + dBzdpc*ampl
    dBzdZ = dBzdZ + dBzdZc*ampl
!
    IF(incore.GT.-1) THEN
! perturbation coil field with plasma shielding:
!
      IF(ipert.EQ.2) THEN
!
        CALL field_fourier(rm,p,zm,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc           &
                          ,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc)
!
        Br = Br + Brc*ampl
        Bp = Bp + Bpc*ampl
        Bz = Bz + Bzc*ampl
!
      ELSEIF(ipert.EQ.3) THEN
!
        CALL field_fourier_derivs(rm,p,zm,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc    &
                                 ,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc)
!
        Br = Br + Brc*ampl
        Bp = Bp + Bpc*ampl
        Bz = Bz + Bzc*ampl
        dBrdR = dBrdR + dBrdRc*ampl
        dBrdp = dBrdp + dBrdpc*ampl
        dBrdZ = dBrdZ + dBrdZc*ampl
        dBpdR = dBpdR + dBpdRc*ampl
        dBpdp = dBpdp + dBpdpc*ampl
        dBpdZ = dBpdZ + dBpdZc*ampl
        dBzdR = dBzdR + dBzdRc*ampl
        dBzdp = dBzdp + dBzdpc*ampl
        dBzdZ = dBzdZ + dBzdZc*ampl
!
      ENDIF
!
    ENDIF 
!
   END IF
!
   RETURN
 END SUBROUTINE field
! ========================================================================
SUBROUTINE field_eq(rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &             
                   ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  USE input_files
  USE field_eq_mod
!
  IMPLICIT NONE
!
  INTEGER :: ierr,i,j
!
  DOUBLE PRECISION :: rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
                     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
!
!-------first call: read data from disk-------------------------------
  IF(icall_eq .LT. 1) THEN
!
  IF(ieqfile.EQ.0) THEN
    CALL read_dimeq0(nrad,nzet)
  ELSE
    CALL read_dimeq1(nrad,nzet)
  ENDIF
!
    ALLOCATE(rad(nrad),zet(nzet))
    ALLOCATE(psi0(nrad,nzet),psi(nrad,nzet))

  IF(ieqfile.EQ.0) THEN
    CALL read_eqfile0(nrad, nzet, psib, btf, rtf, rad, zet, psi)
  ELSE
    CALL read_eqfile1(nrad, nzet, psib, btf, rtf, rad, zet, psi)
  ENDIF
!
! Filtering:
!open(191,file='psi_orig.dat')
!do i=1,nrad
!write (191,*) psi(i,:)
!enddo
!close(191)
    DO i=1,nzet
      CALL window_filter(nrad,nwindow_r,psi(:,i),psi0(:,i))
    ENDDO
!
    DO i=1,nrad
      CALL window_filter(nzet,nwindow_z,psi0(i,:),psi(i,:))
    ENDDO
!open(191,file='psi_filt.dat')
!do i=1,nrad
!write (191,*) psi(i,:)
!enddo
!close(191)
!stop
! End filtering
!     allocate(xi(nzet),f(nzet))
!     npoint = nzet
!     xi = zet
!     do i=1,nrad
!        f = psi(i,:)
!        call leastsq(npoint,xi,f)
!        psi0(i,:) = f
!     enddo
!     deallocate(xi,f)

!     allocate(xi(nrad),f(nrad))
!     npoint = nrad
!     xi = rad
!     do i=1,nzet
!        f = psi0(:,i)
!        call leastsq(npoint,xi,f)
!        psi(:,i) = f
!     enddo
!
    rad = rad*100. ! cm
    zet = zet*100. ! cm
    rtf = rtf*100. ! cm
    psi = psi*1.e8
    psib= psib*1.e8
    btf = btf*1.e4
!
    psi=psi+psib
!
    hrad = rad(2) - rad(1)
    hzet = zet(2) - zet(1)
!
! rectangular domain:
    ALLOCATE( imi(nzet),ima(nzet),jmi(nrad),jma(nrad) )    
    imi = 1
    ima = nrad
    jmi = 1
    jma = nzet
!
!  Computation of the number of data in splpsi
    icp = 0
    DO i=1,nzet
      IF ( imi(i) .GT. 0 .AND. ima(i) .GT. 0 ) THEN
         icp = icp + ima(i) - imi(i) + 1
      ENDIF
    ENDDO
    WRITE(6,*) 'number of points in the table:  ',icp
!
    ALLOCATE( splpsi(6,6,icp), ipoint(nrad,nzet) )
!
    CALL s2dcut(nrad,nzet,hrad,hzet,psi,imi,ima,jmi,jma,icp,splpsi,ipoint)
!
    IF(icall_eq.EQ.-1) THEN
! Quit after initialization with zero field
      Brad=0.d0
      Bphi=0.d0
      Bzet=0.d0
      dBrdR=0.d0
      dBrdp=0.d0
      dBrdZ=0.d0
      dBpdR=0.d0
      dBpdp=0.d0
      dBpdZ=0.d0
      dBzdR=0.d0
      dBzdp=0.d0
      dBzdZ=0.d0
      icall_eq = 1 
      RETURN
    ENDIF
    icall_eq = 1
  ENDIF
!
! ------- end first call ----------------------------------------------
!
    CALL spline(nrad,nzet,rad,zet,hrad,hzet,icp,splpsi,ipoint,rrr,zzz, &
                psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)
!
    Brad = -dpsidz/rrr
    Bzet =  dpsidr/rrr
    Bphi = btf*rtf/rrr
!
! axisymmetric case:
    dBrdp = 0.
    dBpdp = 0.
    dBzdp = 0.
!
    dBpdR = -btf*rtf/rrr**2
    dBpdZ = 0.
!
    dBrdR = -d2psidrdz/rrr+dpsidz/rrr**2
    dBzdZ =  d2psidrdz/rrr
    dBrdZ = -d2psidz2/rrr
    dBzdR =  d2psidr2/rrr-dpsidr/rrr**2
!
  RETURN
END SUBROUTINE field_eq

! ----------- Runov's Original Method --------------------------------
SUBROUTINE read_dimeq0(nrad,nzet)
  USE input_files
  INTEGER :: nrad, nzet

     OPEN(11,file=eqfile)
     READ(11,*)   
     READ(11,*)   
     READ(11,*)
     READ(11,*)   
     READ(11,*)
     READ(11,*)   
     READ(11,*)   
     READ(11,*)
     READ(11,*)   
     READ(11,*) 

     READ(11,111) nrad
     READ(11,111) nzet
111  FORMAT(12x,i3)
 
     CLOSE(11)
  RETURN
END SUBROUTINE read_dimeq0

SUBROUTINE read_eqfile0(nrad, nzet, psib, btf, rtf, rad, zet, psi)
  USE input_files
  INTEGER :: nrad, nzet, dum
  REAL(kind=8) :: psib, btf, rtf
  REAL(kind=8) :: rad(nrad), zet(nzet)
  REAL(kind=8) :: psi(nrad,nzet)

     OPEN(11,file=eqfile)
     READ(11,*)   
     READ(11,*)   
     READ(11,*)
     READ(11,*)   
     READ(11,*)
     READ(11,*)   
     READ(11,*)   
     READ(11,*)
     READ(11,*)   
     READ(11,*) 

     READ(11,111) dum !nrad
     READ(11,111) dum !nzet
     READ(11,112) psib
     READ(11,112) btf
     READ(11,112) rtf
     READ(11,*) 
     READ(11,*) 

     PRINT *, nrad, nzet, psib, btf, rtf
 
!     nrad = nrad - 3
!     nzet = nzet - 4

!!$     read(11,113)dummy,dummy,(rad(i),i=1,nrad)
!!$     read(11,*) 
!!$     read(11,*) 
!!$     read(11,113)dummy,dummy,(zet(i),i=1,nzet)
!!$     read(11,*) 
!!$     read(11,*) 
!!$     read(11,113)(dummy,dummy,(psi(j,k),j=1,nrad),dummy,k=1,2)
!!$     read(11,113)(dummy,dummy,(psi(j,k),j=1,nrad),dummy,k=1,nzet)

     READ(11,113)(rad(i),i=1,nrad)
     READ(11,*) 
     READ(11,*) 
     READ(11,113)(zet(i),i=1,nzet)
     READ(11,*) 
     READ(11,*) 
     READ(11,113)((psi(j,k),j=1,nrad),k=1,nzet)

!!$     do k=1,nzet
!!$        write(41,*)(psi(j,k),j=1,nrad)
!!$     enddo
    CLOSE(11)
    RETURN

111  FORMAT(12x,i3)
112  FORMAT(12x,f21.2)
113  FORMAT(5(e17.4))
END SUBROUTINE read_eqfile0


! ----------- Read gfile directly --------------------------------
SUBROUTINE read_dimeq1(nwEQD,nhEQD)
  USE input_files
  IMPLICIT NONE
  INTEGER :: nwEQD, nhEQD,i
  INTEGER :: idum
  CHARACTER*10 CASE(6)
!
     OPEN(unit=iunit,file=TRIM(gfile),status='old',action='read')
     READ(iunit,2000)(CASE(i),i=1,6),idum,nwEQD,nhEQD
     CLOSE(iunit)
     OPEN(unit=iunit,file='out.06')
     WRITE(iunit,*) 'READ_DIMEQ1: ',nwEQD,nhEQD
     CLOSE(iunit)
  RETURN

2000  FORMAT(6a8,3i4)
55    PRINT *, 'READ_EQDIM1: Early EOF in',TRIM(gfile); STOP
250   PRINT *, 'READ_EQDIM1: Error reading ',TRIM(gfile); STOP
END SUBROUTINE read_dimeq1


SUBROUTINE read_eqfile1(nwEQD,nhEQD,psiSep, bt0, rzero, rad, zet, psiRZ)
  USE input_files
  IMPLICIT NONE
  INTEGER :: nwEQD, nhEQD
  INTEGER :: gunit, idum
  CHARACTER*10 CASE(6)
  INTEGER :: i,j
  REAL (kind=8) :: xdim,zdim,r1,zmid,rmaxis,zmaxis,xdum
  REAL (kind=8) :: bt0, rzero, plas_cur, psiAxis, psiSep
  REAL (kind=8), DIMENSION(nwEQD) :: fpol,pres,ffprim,pprime,qpsi
  REAL (kind=8), DIMENSION(nwEQD,nhEQD) :: psiRZ
  REAL (kind=8) :: rad(nwEQD), zet(nhEQD)

  INTEGER :: n_bndyxy,nlimEQD
  REAL (kind=8), DIMENSION(:), ALLOCATABLE :: LCFS, limEQD

      gunit=iunit

      OPEN(unit=gunit,file=TRIM(gfile),status='old',action='read')

! Equilibrium Parameters
      READ(gunit,2000)(CASE(i),i=1,6),idum,nwEQD,nhEQD
      WRITE(*,*) 'READ_EQFILE1: ',TRIM(gfile),nwEQD,nhEQD
      READ(gunit,2010,END=55,err=250)xdim,zdim,rzero,r1,zmid
      WRITE(*,*) xdim, zdim, rzero, r1, zmid
      READ(gunit,2010,END=55,err=250)rmaxis,zmaxis,psiAxis,psiSep,bt0
      WRITE(*,*) rmaxis,zmaxis,psiAxis,psiSep,bt0
      READ(gunit,2010,END=55,err=250)plas_cur,psiAxis,xdum,rmaxis,xdum
      WRITE(*,*) plas_cur,psiAxis,xdum,rmaxis,xdum
      READ(gunit,2010,END=55,err=250)zmaxis,xdum,psiSep,xdum,xdum
      WRITE(*,*) zmaxis,xdum,psiSep,xdum,xdum
      READ(gunit,2010,END=55,err=250)(fpol(i),i=1,nwEQD)
      READ(gunit,2010,END=55,err=250)(pres(i),i=1,nwEQD)
      READ(gunit,2010,END=55,err=250)(ffprim(i),i=1,nwEQD)
      READ(gunit,2010,END=55,err=250)(pprime(i),i=1,nwEQD)
      READ(gunit,2010,END=55,err=250)((psiRZ(i,j),i=1,nwEQD),j=1,nhEQD)
      READ(gunit,2010,END=55,err=250)(qpsi(i),i=1,nwEQD)
      PRINT *, 'Equilibrium Done.', TRIM(gfile)
! Boundary Data
      READ(gunit,*,END=55,err=250)n_bndyxy,nlimEQD    
      ALLOCATE(LCFS(2*n_bndyxy))
      ALLOCATE(limEQD(2*nlimEQD))               
      READ(gunit,2010,END=55,err=250)(LCFS(i),i=1,2*n_bndyxy)
      READ(gunit,2010,END=55,err=250)(limEQD(i),i=1,2*nlimEQD)
!      print *, 'Boundary Done.'
      CLOSE(gunit)

      CALL set_eqcoords(nwEQD,nhEQD,xdim,zdim,r1,zmid,rad,zet)           
  RETURN

2000  FORMAT(6a8,3i4)
2010  FORMAT(5(e16.9))
55    PRINT *, 'READ_EQFILE1: Early EOF in',TRIM(gfile); STOP
250   PRINT *, 'READ_EQFILE1: Error reading ',TRIM(gfile); STOP

END SUBROUTINE read_eqfile1


SUBROUTINE set_eqcoords(nwEQD,nhEQD,xdim,zdim,r1,zmid,rad,zet)
  IMPLICIT NONE
  INTEGER :: j,k,nwEQD,nhEQD
  REAL (kind=8) :: xdim,zdim,r1,zmid,z1
  REAL (kind=8) :: rad(nwEQD), zet(nhEQD)

  DO j=1,nwEQD
    rad(j) = r1 + (j-1)*(xdim/(nwEQD-1))
  END DO

  z1 = zmid - zdim/2.0 ! check this definition wrt zmid
  DO k=1,nhEQD ! runov chooses lower, probe chooses upper
    zet(k) = z1 + (k-1)*(zdim/(nhEQD-1))
  END DO

!      print *, 'set_coords done.'
!      print *, 'rad'
!      print 2010, (rad(j),j=1,nwEQD)
!      print *, ''
!      print *, 'zet'
!      print 2010, (zet(k),k=1,nhEQD)

  RETURN
2010  FORMAT(5(e16.9))
END SUBROUTINE set_eqcoords

! ===========================================================================
!
SUBROUTINE field_c(rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &             
                  ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  USE input_files
  USE field_c_mod
!
  IMPLICIT NONE
!
  DOUBLE PRECISION, PARAMETER :: pi=3.14159265358979d0
!
  DOUBLE PRECISION :: rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
                     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  DOUBLE PRECISION :: rmin,pmin,zmin,rmax,pmax,zmax,hrm1,hpm1,hzm1
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Br,Bp,Bz
!
!-------first call: read data from disk-------------------------------
  IF(icall_c .LT. 1) THEN
    PRINT *,'coils: file type = ',icftype
    IF(icftype.EQ.1) THEN
      nr=129  !64
      np=37   !37
      nz=129  !64
    ELSEIF(icftype.EQ.2) THEN
      nr=129
      np=33
      nz=129
    ELSEIF(icftype.EQ.3) THEN
      nr=129
      np=37
      nz=131
      icftype=1
    ELSEIF(icftype.EQ.4) THEN
      CALL read_sizes(nr,np,nz)
    ELSE
      PRINT *,'field_c: unknown coil file type'
      STOP
    ENDIF
    ALLOCATE(Br(nr,np,nz),Bp(nr,np,nz),Bz(nr,np,nz))
!
!    call read_field0(rad,phi,zet,rmin,pmin,zmin,hrm1,hpm1,hzm1,Br,Bp,Bz)
!    call read_field1(icftype,nr,np,nz,rmin,rmax,pmin,pmax,zmin,zmax,Br,Bp,Bz)
!
!
    IF(icftype.LT.4) THEN
      CALL read_field2(icftype,nr,np,nz,rmin,rmax,pmin,pmax,zmin,zmax,Br,Bp,Bz)
    ELSE
      CALL read_field4(nr,np,nz,rmin,rmax,pmin,pmax,zmin,zmax,Br,Bp,Bz)
    ENDIF
!
    PRINT *,'coils: nr,np,nz = ',nr,np,nz
    PRINT *,'coils: rmin,rmax = ',rmin,rmax
    PRINT *,'coils: zmin,zmax = ',zmin,zmax
    PRINT *,'coils: pmin,pmax = ',pmin,pmax
!
    CALL vector_potentials(nr,np,nz,ntor,rmin,rmax,pmin,pmax,zmin,zmax,br,bp,bz)
!
    DEALLOCATE(Br,Bp,Bz)
!
    IF(icall_c.EQ.-1) THEN
! Quit after initialization with zero field
      Brad=0.d0
      Bphi=0.d0
      Bzet=0.d0
      dBrdR=0.d0
      dBrdp=0.d0
      dBrdZ=0.d0
      dBpdR=0.d0
      dBpdp=0.d0
      dBpdZ=0.d0
      dBzdR=0.d0
      dBzdp=0.d0
      dBzdZ=0.d0
      icall_c = 1
      RETURN
    ENDIF
    icall_c = 1
  ENDIF
  !------- end first call ----------------------------------------------
!
  CALL  field_divfree(rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ    &
                     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  RETURN
END SUBROUTINE field_c



! ===========================================================================
SUBROUTINE read_field0(rad,phi,zet,rmin,pmin,zmin,hrm1,hpm1,hzm1,Br,Bp,Bz)
!
  USE input_files
  PARAMETER(nr=64,np=37,nz=64)
  !! Modification by Andreas F. Martitsch (16.07.2015)
  ! previous (Warning: Change of value in conversion from REAL(8) to REAL(4))
  !REAL, PARAMETER :: pi=3.14159265358979d0
  ! corrected:
  REAL, PARAMETER :: pi=3.14159265358979e0
  !! End Modification by Andreas F. Martitsch (16.07.2015)
  PARAMETER (mp=4) ! power of Lagrange's polynomial =3
  DIMENSION Bz(nr,np,nz)
  DIMENSION Br(nr,np,nz),Bp(nr,np,nz)
  DIMENSION rad(nr), phi(np), zet(nz)
  DIMENSION xp(mp),yp(mp),zp(mp),fp(mp,mp,mp)
  INTEGER indx(mp), indy(mp), indz(mp)
  DATA icall/0/
  SAVE
!
!-------first call: read data from disk-------------------------------
     OPEN(1,file=cfile,status='old',action='read')  
     READ(1,*)   
     READ(1,*)
     READ(1,*)   
     READ(1,*)
     READ(1,*)

!---Input B      -->T = V*s/m/m
     DO j=1,np-1 !only npmax-1 points are given
        DO k=nz,1,-1  !reverse order of probe data
           DO i=1,nr
              READ(1,*) Br(i,j,k), Bp(i,j,k), Bz(i,j,k)
              
              !! Modification by Andreas F. Martitsch (16.07.2015)
              ! previous (Warning: Change of value in conversion from REAL(8) to REAL(4))
              !Br(i,j,k) = Br(i,j,k)*1.d4
              !Bp(i,j,k) = Bp(i,j,k)*1.d4
              !Bz(i,j,k) = Bz(i,j,k)*1.d4
              ! corrected:
              Br(i,j,k) = Br(i,j,k)*1.e4
              Bp(i,j,k) = Bp(i,j,k)*1.e4
              Bz(i,j,k) = Bz(i,j,k)*1.e4
              !! End Modification by Andreas F. Martitsch (16.07.2015)

           ENDDO
           READ(1,*)
        ENDDO
        READ(1,*)
     ENDDO
     CLOSE(1)
     !
     rmin = 84.
     rmax = 254.
     zmin = -160.
     zmax = 160.
     pmin = 0.
     pmax = 2.*pi

 
     hrad = (rmax - rmin)/(nr-1)  
     hphi = (pmax - pmin)/(np-1)
     hzet = (zmax - zmin)/(nz-1)

     DO i=1,nr
        rad(i) = rmin + hrad*(i-1)
     ENDDO
     DO i=1,np
        phi(i) = pmin + hphi*(i-1)
     ENDDO
     DO i=1,nz
        zet(i) = zmin + hzet*(i-1)
     ENDDO

     DO i=1,nr
        DO k=1,nz
           Br(i,np,k) = Br(i,1,k)
           Bp(i,np,k) = Bp(i,1,k)
           Bz(i,np,k) = Bz(i,1,k)
        ENDDO
     ENDDO
END SUBROUTINE read_field0
!
SUBROUTINE read_field1(icftype,nr,np,nz,rmin,rmax,pmin,pmax,zmin,zmax,Br,Bp,Bz)
!
  USE input_files
!
  IMPLICIT NONE
!
  DOUBLE PRECISION, PARAMETER :: pi=3.14159265358979d0
!
  INTEGER :: nr,np,nz,i,j,k,icftype
  DOUBLE PRECISION :: rmin,pmin,zmin,rmax,pmax,zmax,xdim,zdim,zmid,dum
  DOUBLE PRECISION, DIMENSION(nr,np,nz) :: Br,Bp,Bz
!
  OPEN(iunit,file=TRIM(pfile),status='old',action='read')
  READ(iunit,*)   
  READ(iunit,*)   
  READ(iunit,*)   
  READ(iunit,*)  !PROBE 
  READ(iunit,*)   
  READ(iunit,*)   
  IF(icftype.EQ.2) THEN
    READ(iunit,*)    !New Format
    READ(iunit,*)    !New Format
  ENDIF

!---Input B      -->T = V*s/m/m
  DO j=1,np-1   !only npmax-1 points are given
     DO k=nz,1,-1  !reverse order of probe data
        DO i=1,nr
           IF(icftype.EQ.1) THEN
!					Old Format
             READ(iunit,*) Br(i,j,k),Bp(i,j,k),Bz(i,j,k)
           ELSEIF(icftype.EQ.2) THEN
!					New Format
             READ(iunit,*) dum,dum,dum,Br(i,j,k),Bp(i,j,k),Bz(i,j,k),dum,dum
           ENDIF
!
				  !Convert to CGS
           Br(i,j,k) = Br(i,j,k)*1.d4
           Bp(i,j,k) = Bp(i,j,k)*1.d4
           Bz(i,j,k) = Bz(i,j,k)*1.d4
        ENDDO
        READ(iunit,*)
     ENDDO
     READ(iunit,*)
  ENDDO
  CLOSE(iunit)
!
  xdim=170.d0
  rmin=84.d0
  rmax=rmin+xdim
!
  pmin = 0.
  pmax = 2.*pi
!
  zdim=320.d0
  zmid=0.d0
  zmin=zmid - zdim/2.d0
  zmax=zmid + zdim/2.d0
!     
  DO i=1,nr
     DO k=1,nz
        Br(i,np,k) = Br(i,1,k)
        Bp(i,np,k) = Bp(i,1,k)
        Bz(i,np,k) = Bz(i,1,k)
     ENDDO
  ENDDO
END SUBROUTINE read_field1
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
SUBROUTINE stretch_coords(r,z,rm,zm)
  USE input_files, ONLY : iunit,convexfile
  IMPLICIT NONE
  INTEGER icall, i, j, nrz ! number of points "convex wall" in input file
  INTEGER, PARAMETER :: nrzmx=100 ! possible max. of nrz
  INTEGER, PARAMETER :: nrhotht=360 
  INTEGER :: iflag
  REAL(kind=8), PARAMETER :: pi = 3.14159265358979d0
  REAL(kind=8) R0,Rw, Zw, htht, Rl, Zl, a, b, r, z, rm, zm, rho, tht, rho_c, delta, dummy
  REAL(kind=8), DIMENSION(0:1000):: rad_w, zet_w ! points "convex wall" 
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: rho_w, tht_w 
  REAL(kind=8), DIMENSION(nrhotht) :: rho_wall, tht_wall ! polar coords of CW 
  DATA icall /0/, delta/1./
  SAVE
!----------- 1st call --------------------------------------------------------
  IF(icall .EQ. 0) THEN
     icall = 1
     nrz = 0
     rad_w = 0.
     zet_w = 0.
     OPEN(iunit,file=TRIM(convexfile))
     DO i=1,nrzmx
        READ(iunit,*,END=10)rad_w(i),zet_w(i)
        nrz = nrz + 1
     ENDDO
10   CONTINUE
     CLOSE(iunit)

     ALLOCATE(rho_w(0:nrz+1), tht_w(0:nrz+1))
     R0 = (MAXVAL(rad_w(1:nrz)) +  MINVAL(rad_w(1:nrz)))*0.5d0
     DO i=1,nrz
        rho_w(i) = SQRT( (rad_w(i)-R0)**2 + zet_w(i)**2 )
        tht_w(i) = ATAN2(zet_w(i),(rad_w(i)-R0))
        IF(tht_w(i) .LT. 0.d0) tht_w(i) = tht_w(i) + 2.d0*pi
     ENDDO
!
     DO
       iflag=0
       DO i=1,nrz-1
         IF(tht_w(i).GT.tht_w(i+1)) THEN
           iflag=1
           dummy=rad_w(i+1)
           rad_w(i+1)=rad_w(i)
           rad_w(i)=dummy
           dummy=zet_w(i+1)
           zet_w(i+1)=zet_w(i)
           zet_w(i)=dummy
           dummy=rho_w(i+1)
           rho_w(i+1)=rho_w(i)
           rho_w(i)=dummy
           dummy=tht_w(i+1)
           tht_w(i+1)=tht_w(i)
           tht_w(i)=dummy
         ENDIF
       ENDDO
       IF(iflag.EQ.0) EXIT
     ENDDO
     rad_w(0)=rad_w(nrz)
     zet_w(0)=zet_w(nrz)
     tht_w(0)=tht_w(nrz)-2.d0*pi
     rho_w(0)=rho_w(nrz)
     rad_w(nrz+1)=rad_w(1)
     zet_w(nrz+1)=zet_w(1)
     tht_w(nrz+1)=tht_w(1)+2.d0*pi
     rho_w(nrz+1)=rho_w(1)
!
     htht = 2.*pi/(nrhotht-1)
     DO i=2,nrhotht
        tht_wall(i) = htht*(i-1)
        DO j=0,nrz
           IF(tht_wall(i).GE.tht_w(j) .AND. tht_wall(i).LE.tht_w(j+1)) THEN
              IF( ABS((rad_w(j+1) - rad_w(j))/rad_w(j)) .GT. 1.e-3) THEN
                 a = (zet_w(j+1) - zet_w(j))/(rad_w(j+1) - rad_w(j))
                 b = zet_w(j) - a*(rad_w(j) - R0)
                 Rw = b/(TAN(tht_wall(i)) - a) + R0
                 Zw = a*(Rw - R0) + b
              ELSE
                 a = (rad_w(j+1) - rad_w(j))/(zet_w(j+1) - zet_w(j))
                 b = rad_w(j)-R0 - a*zet_w(j)
                 Zw = b/(1./TAN(tht_wall(i)) - a)
                 Rw = a*Zw + b + R0
              ENDIF
           ENDIF
        ENDDO
        rho_wall(i) = SQRT((Rw-R0)**2 + Zw**2)
     ENDDO
     tht_wall(1) = 0.
     rho_wall(1) = rho_wall(nrhotht)
!  do i=1,nrhotht
!     write(19,*)tht_wall(i), rho_wall(i)   
!  enddo
  ENDIF
!----------- end of the 1st call --------------------------------------------
  rm = r
  zm = z
  rho = SQRT((r-R0)**2 + z**2)
  IF(z.EQ.0.d0.AND.r-r0.EQ.0.d0) THEN
    tht=0.d0
  ELSE
    tht = ATAN2(z,(r-R0))
  ENDIF
  IF(tht .LT. 0.) tht = tht + 2.*pi
  i = INT(tht/htht) + 1
  rho_c = (rho_wall(i+1) - rho_wall(i))/(tht_wall(i+1) - tht_wall(i))   &
       *(tht - tht_wall(i)) + rho_wall(i)
!print *,rho,rho_c,i,tht
  IF(rho .GE. rho_c) THEN
     rho = rho_c + delta*ATAN2((rho-rho_c), delta)
     rm = rho*COS(tht) + R0
     zm = rho*SIN(tht)
  ENDIF

  RETURN
END SUBROUTINE stretch_coords
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  SUBROUTINE inthecore(R,Z)
!
  USE inthecore_mod
  USE input_files,  ONLY : iunit,fluxdatapath
  USE field_eq_mod, ONLY : psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
!
  IMPLICIT NONE
!
  INTEGER :: i
  DOUBLE PRECISION :: R,Z,rho2,thet,scalp,xx,yy
  DOUBLE PRECISION :: weight,dweight,ddweight
  DOUBLE PRECISION, DIMENSION(4) :: x,y
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ri,zi
!
  IF(prop) THEN
    prop=.FALSE.
    OPEN(iunit,file=TRIM(fluxdatapath)//'/separ.dat')
    READ(iunit,*) x(1:2)
    psi_sep=x(2)+(x(1)-x(2))*(1.d0-epssep)
    psi_cut=x(2)+(x(1)-x(2))*cutoff
    sigpsi=SIGN(1.d0,psi_sep-psi_cut)
    npoi=0
    DO WHILE(npoi.GE.0)
      npoi=npoi+1
      READ(iunit,*,END=1) rc
    ENDDO
 1  ALLOCATE(ri(0:npoi),zi(0:npoi),rho2i(0:npoi),theti(0:npoi))
    ri=0.d0
    zi=0.d0
    CLOSE(iunit)
    OPEN(iunit,file=TRIM(fluxdatapath)//'/separ.dat')
    READ(iunit,*)
    DO i=1,npoi-1
      READ(iunit,*) ri(i),zi(i)
    ENDDO
    CLOSE(iunit)
    rc=SUM(ri(1:npoi-1))/(npoi-1)
    zc=SUM(zi(1:npoi-1))/(npoi-1)
    rho2i=(ri-rc)**2+(zi-zc)**2
    theti=ATAN2(zi-zc,ri-rc)
    sig=theti(2)-theti(1)
    DO i=2,npoi-2
      IF((theti(i+1)-theti(i))*sig.LT.0.d0) THEN
        ijumpb=i
        EXIT
      ENDIF
    ENDDO
    twopi=8.d0*ATAN2(1.d0,1.d0)
    ri(1:npoi-1-ijumpb)=rho2i(ijumpb+1:npoi-1)
    ri(npoi-ijumpb:npoi-1)=rho2i(1:ijumpb)
    rho2i=ri
    ri(1:npoi-1-ijumpb)=theti(ijumpb+1:npoi-1)
    ri(npoi-ijumpb:npoi-1)=theti(1:ijumpb)
    theti=ri
    DEALLOCATE(ri,zi)
    sig=theti(2)-theti(1)
    rho2i(npoi)=rho2i(1)
    theti(npoi)=theti(1)+SIGN(twopi,sig)
    rho2i(0)=rho2i(npoi-1)
    theti(0)=theti(npoi-1)-SIGN(twopi,sig)
  ENDIF
!
  rho2=(r-rc)**2+(z-zc)**2
  thet=ATAN2(z-zc,r-rc)
!
  ibeg=0
  iend=npoi
  DO WHILE(ibeg+1.LT.iend)
    i=(ibeg+iend)/2
    IF((thet-theti(i))*sig.GT.0.d0) THEN
      ibeg=i
    ELSE
      iend=i
    ENDIF
  ENDDO
  x=theti(ibeg-1:iend+1)
  y=rho2i(ibeg-1:iend+1)
!
  xx=thet
  yy=y(1)*(xx-x(2))/(x(1)-x(2))*(xx-x(3))/(x(1)-x(3))*(xx-x(4))/(x(1)-x(4)) &
    +y(2)*(xx-x(3))/(x(2)-x(3))*(xx-x(4))/(x(2)-x(4))*(xx-x(1))/(x(2)-x(1)) &
    +y(3)*(xx-x(4))/(x(3)-x(4))*(xx-x(1))/(x(3)-x(1))*(xx-x(2))/(x(3)-x(2)) &
    +y(4)*(xx-x(1))/(x(4)-x(1))*(xx-x(2))/(x(4)-x(2))*(xx-x(3))/(x(4)-x(3))
!
  IF(rho2.GT.yy) THEN
    incore=-1
    RETURN
  ELSEIF((psif-psi_cut)*sigpsi.LT.0.d0) THEN
    incore=1
    RETURN
  ENDIF
!
  incore=0
!
  CALL localizer(psi_cut,psi_sep,psif,weight,dweight,ddweight)
!
  plaf=weight
  dpladr=dweight*dpsidr
  dpladz=dweight*dpsidz
  d2pladr2=ddweight*dpsidr**2+dweight*d2psidr2
  d2pladrdz=ddweight*dpsidr*dpsidz+dweight*d2psidrdz
  d2pladz2=ddweight*dpsidz**2+dweight*d2psidz2
!
  vacf=1.d0-plaf
  dvacdr=-dpladr
  dvacdz=-dpladz
  d2vacdr2=-d2pladr2
  d2vacdrdz=-d2pladrdz
  d2vacdz2=-d2pladz2
!
  RETURN
  END SUBROUTINE inthecore
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  SUBROUTINE localizer(x1,x2,x,weight,dweight,ddweight)
!
  IMPLICIT NONE
!
  DOUBLE PRECISION, PARAMETER :: c1=-6.283185307179586d0,c2=-1.414213562373095d0
!
  DOUBLE PRECISION :: x1,x2,x,t,weight,dweight,ddweight,exin
!
  t=(x-x1)/(x2-x1)
!
  IF(t.LE.0.d0) THEN
    weight=1.d0
    dweight=0.d0
    ddweight=0.d0
  ELSEIF(t.GE.1.d0) THEN
    weight=0.d0
    dweight=0.d0
    ddweight=0.d0
  ELSE
    exin=EXP(c2/t)
    weight=EXP(c1/(1.d0-t)*exin)
    dweight=weight*c1*(1.d0/(1.d0-t)-c2/t**2)*exin/(1.d0-t)
    ddweight=dweight*c1*(1.d0/(1.d0-t)-c2/t**2)*exin/(1.d0-t) &
            +weight*c1*(1.d0/(1.d0-t)**2+2.d0*c2/t**3)*exin/(1.d0-t) &
            +weight*c1*(1.d0/(1.d0-t)-c2/t**2)**2*exin/(1.d0-t)
  ENDIF
!
  dweight=dweight/(x2-x1)
  ddweight=ddweight/(x2-x1)**2
!
  RETURN
  END SUBROUTINE localizer
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  SUBROUTINE window_filter(n,nw,arr_in,arr_out)
!
  IMPLICIT NONE
!
  INTEGER :: n,nw,nwa,i
  DOUBLE PRECISION, DIMENSION(n) :: arr_in,arr_out
!
  DO i=1,n
    nwa=MIN(nw,i-1,n-i)
    arr_out(i)=SUM(arr_in(i-nwa:i+nwa))/(2*nwa+1)
  ENDDO
!
  RETURN
  END SUBROUTINE window_filter
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
SUBROUTINE read_field2(icftype,nr,np,nz,rmin,rmax,pmin,pmax,zmin,zmax,Br,Bp,Bz)
!
  USE input_files
!
  IMPLICIT NONE
!
  DOUBLE PRECISION, PARAMETER :: pi=3.14159265358979d0
!
  INTEGER :: nr,np,nz,i,j,k,icftype
  DOUBLE PRECISION :: rmin,pmin,zmin,rmax,pmax,zmax,xdim,zdim,zmid,dum
  DOUBLE PRECISION, DIMENSION(nr,np,nz) :: Br,Bp,Bz
!
  OPEN(iunit,file=TRIM(pfile),status='old',action='read')

!---Input B      -->T = V*s/m/m
  DO j=1,np-1   !only npmax-1 points are given
     DO k=nz,1,-1  !reverse order of probe data
        DO i=1,nr
           IF(icftype.EQ.1) THEN
!                                       Old Format
             READ(iunit,*) Br(i,j,k),Bp(i,j,k),Bz(i,j,k)
           ELSEIF(icftype.EQ.2) THEN
!                                       New Format
             READ(iunit,*) dum,dum,dum,Br(i,j,k),Bp(i,j,k),Bz(i,j,k),dum,dum
           ENDIF
!
                                  !Convert to CGS
           Br(i,j,k) = Br(i,j,k)*1.d4
           Bp(i,j,k) = Bp(i,j,k)*1.d4
           Bz(i,j,k) = Bz(i,j,k)*1.d4
        ENDDO
     ENDDO
  ENDDO
  CLOSE(iunit)
!
  xdim=300.d0
  rmin=100.d0
  rmax=rmin+xdim
!
  pmin = 0.
  pmax = 2.*pi
!
  zdim=400.d0
  zmid=0.d0
  zmin=zmid - zdim/2.d0
  zmax=zmid + zdim/2.d0
!
  DO i=1,nr
     DO k=1,nz
        Br(i,np,k) = Br(i,1,k)
        Bp(i,np,k) = Bp(i,1,k)
        Bz(i,np,k) = Bz(i,1,k)
     ENDDO
  ENDDO
END SUBROUTINE read_field2
!
SUBROUTINE read_sizes(nr,np,nz)
!
  USE input_files, ONLY : iunit,pfile
!
  IMPLICIT NONE
  INTEGER :: nr,np,nz
!
  OPEN(iunit,file=pfile)
  READ(iunit,*) nr,np,nz
  CLOSE(iunit)
!
END SUBROUTINE read_sizes
!
SUBROUTINE read_field4(nr,np,nz,rmin,rmax,pmin,pmax,zmin,zmax,Br,Bp,Bz)
!
  USE input_files, ONLY : iunit,pfile
!
  IMPLICIT NONE
!
  INTEGER :: nr,np,nz,i,j,k
  DOUBLE PRECISION :: rmin,rmax,pmin,pmax,zmin,zmax
  DOUBLE PRECISION, DIMENSION(nr,np,nz)       :: Br,Bp,Bz
!
  OPEN(iunit,file=pfile) 
  READ(iunit,*) nr,np,nz
  READ(iunit,*) rmin,rmax
  READ(iunit,*) pmin,pmax
  READ(iunit,*) zmin,zmax
  DO i=1,nr
     DO j=1,np
        DO k=1,nz
           READ(iunit,*) Br(i,j,k),Bp(i,j,k),Bz(i,j,k)
        ENDDO
     ENDDO
  ENDDO
  CLOSE(iunit)
!
END SUBROUTINE read_field4
