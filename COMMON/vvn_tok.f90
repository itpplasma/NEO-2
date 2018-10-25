
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
SUBROUTINE stevvo_0(RT0,R0i,L1i,cbfi,BY0i,bf0)

  use circtokfield, only : rbig, btor

  IMPLICIT NONE

  INTEGER :: L1i
  DOUBLE PRECISION :: RT0,R0i,cbfi,BY0i,bf0

  RT0=rbig
  L1i=1
  bf0=btor

END SUBROUTINE stevvo_0
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
SUBROUTINE GBas_0(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,                      &
                   dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!  subroutine field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,                      &
!                   dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  USE circtokfield
  USE coordinates_of_old_model
!
  IMPLICIT NONE
!
  DOUBLE PRECISION :: r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,                   &
                      dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ  
!
! rqt,thqt - small radius and theta of quasitoroidal coordinates
  DOUBLE PRECISION :: rqt,thqt
  DOUBLE PRECISION :: x,x2,z2,rho,rqt2,rqt3,rqt4,rho2,btor_rbig
! tfl,theta - toroidal flux and theta of flux coordinates
  DOUBLE PRECISION :: tfl,theta
  DOUBLE PRECISION :: drqt_dr,drqt_dz
  DOUBLE PRECISION :: ddrqt_dr2,ddrqt_dr_dz,ddrqt_dz2
  DOUBLE PRECISION :: dthqt_dr,dthqt_dz
  DOUBLE PRECISION :: ddthqt_dr2,ddthqt_dr_dz,ddthqt_dz2
!
  DOUBLE PRECISION :: dtfl,ddtfl
  DOUBLE PRECISION :: dtheta_drqt,dtheta_dthqt
  DOUBLE PRECISION :: ddtheta_drqt2,ddtheta_drqt_dthqt,ddtheta_dthqt2
!
  DOUBLE PRECISION :: dtfl_dr,dtfl_dz
  DOUBLE PRECISION :: ddtfl_dr2,ddtfl_dr_dz,ddtfl_dz2
  DOUBLE PRECISION :: dtheta_dr,dtheta_dz
  DOUBLE PRECISION :: ddtheta_dr2,ddtheta_dr_dz,ddtheta_dz2
!
  DOUBLE PRECISION :: psi,dpsi_dtfl,dpsi_dtheta
  DOUBLE PRECISION :: ddpsi_dtfl2,ddpsi_dtfl_dtheta,ddpsi_dtheta2
  DOUBLE PRECISION :: ddpsi_dtfl_dp,ddpsi_dtheta_dp
!
  DOUBLE PRECISION :: tfl_e
!
!
! Relations between the cylindrical and flux coordinates in
! the tokamak with circular concentric magnetic surfaces (zero beta)
!
  x=r-rbig
  IF(z.GT.rbig.OR.x.GT.rbig) THEN
    PRINT *,'r > R_0'
    RETURN
  ENDIF
!
  x2=x**2
  z2=z**2
  rqt2=x2+z2
  rqt=SQRT(rqt2)
  rqt3=rqt2*rqt
  rqt4=rqt2**2
  thqt=ATAN2(z,x)
!
  rho2=rbig**2-rqt2
  rho=SQRT(rho2)
!
  drqt_dr=x/rqt
  drqt_dz=z/rqt
  ddrqt_dr2=z2/rqt3
  ddrqt_dr_dz=-x*z/rqt3
  ddrqt_dz2=x2/rqt3
!
  dthqt_dr=-z/rqt2
  dthqt_dz=x/rqt2
  ddthqt_dr2=2.d0*x*z/rqt4
  ddthqt_dr_dz=(z2-x2)/rqt4
  ddthqt_dz2=-ddthqt_dr2
!
  btor_rbig=btor*rbig
  tfl=btor_rbig*(rbig-rho)
  dtfl=btor_rbig*rqt/rho
  ddtfl=btor*(rbig/rho)**3
!
  theta=2.d0*ATAN(SQRT((rbig-rqt)/(rbig+rqt))*TAN(0.5d0*thqt))
  dtheta_drqt=-rbig*z/(R*rho*rqt)
  dtheta_dthqt=rho/R
  ddtheta_drqt2=(x/(R*rqt2)-1.d0/rho2)*rbig*z/(R*rho)
  ddtheta_drqt_dthqt=-(rqt/rho+x*rho/(R*rqt))/R
  ddtheta_dthqt2=rho*z/R**2
!
  dtfl_dr=dtfl*drqt_dr
  dtfl_dz=dtfl*drqt_dz
  ddtfl_dr2=ddtfl*drqt_dr**2+dtfl*ddrqt_dr2
  ddtfl_dr_dz=ddtfl*drqt_dr*drqt_dz+dtfl*ddrqt_dr_dz
  ddtfl_dz2=ddtfl*drqt_dz**2+dtfl*ddrqt_dz2
  dtheta_dr=dtheta_drqt*drqt_dr+dtheta_dthqt*dthqt_dr
  dtheta_dz=dtheta_drqt*drqt_dz+dtheta_dthqt*dthqt_dz
  ddtheta_dr2=ddtheta_drqt2*drqt_dr**2                                      &
             +2.d0*ddtheta_drqt_dthqt*drqt_dr*dthqt_dr                      &
             +ddtheta_dthqt2*dthqt_dr**2                                    &
             +dtheta_drqt*ddrqt_dr2+dtheta_dthqt*ddthqt_dr2
  ddtheta_dr_dz=ddtheta_drqt2*drqt_dr*drqt_dz                               &
             +ddtheta_drqt_dthqt*(drqt_dr*dthqt_dz+drqt_dz*dthqt_dr)        &
             +ddtheta_dthqt2*dthqt_dr*dthqt_dz                              &
             +dtheta_drqt*ddrqt_dr_dz+dtheta_dthqt*ddthqt_dr_dz
  ddtheta_dz2=ddtheta_drqt2*drqt_dz**2                                      &
             +2.d0*ddtheta_drqt_dthqt*drqt_dz*dthqt_dz                      &
             +ddtheta_dthqt2*dthqt_dz**2                                    &
             +dtheta_drqt*ddrqt_dz2+dtheta_dthqt*ddthqt_dz2
!
! end of relations for the tokamak with circular concentric surfaces
!
! Toroidal field
!
    Bp=btor_rbig/R
    dBpdR=-Bp/R
    dBpdZ=0.d0
    dBpdp=0.d0
!
! Poloidal field
!
! Hamiltonian form (H=psi)
!
! tfl_e - toroidal flux at the edge
  tfl_e=btor_rbig*(rbig-SQRT(rbig**2-rsmall**2))
  CALL hamilt_field(tfl_e,tfl,theta,p,dpsi_dtfl,dpsi_dtheta,              &
                    ddpsi_dtfl2,ddpsi_dtfl_dtheta,ddpsi_dtheta2,          &
                    ddpsi_dtfl_dp,ddpsi_dtheta_dp)
!
  Br=-(dpsi_dtfl*dtfl_dz+dpsi_dtheta*dtheta_dz)/R
  Bz=(dpsi_dtfl*dtfl_dr+dpsi_dtheta*dtheta_dr)/R
!
  dBrdR=-(ddpsi_dtfl2*dtfl_dr*dtfl_dz                                     &
        + ddpsi_dtfl_dtheta*(dtfl_dr*dtheta_dz+dtfl_dz*dtheta_dr)         &
        + ddpsi_dtheta2*dtheta_dr*dtheta_dz                               &
        + dpsi_dtfl*ddtfl_dr_dz+dpsi_dtheta*ddtheta_dr_dz)/R-Br/R
  dBrdZ=-(ddpsi_dtfl2*dtfl_dz**2+2.d0*ddpsi_dtfl_dtheta*dtfl_dz*dtheta_dz &
        + ddpsi_dtheta2*dtheta_dz**2                                      &
        + dpsi_dtfl*ddtfl_dz2+dpsi_dtheta*ddtheta_dz2)/R
  dBrdp=-(ddpsi_dtfl_dp*dtfl_dz+ddpsi_dtheta_dp*dtheta_dz)/R
!
  dBzdR= (ddpsi_dtfl2*dtfl_dr**2+2.d0*ddpsi_dtfl_dtheta*dtfl_dr*dtheta_dr &
       +  ddpsi_dtheta2*dtheta_dr**2                                      &
       +  dpsi_dtfl*ddtfl_dr2+dpsi_dtheta*ddtheta_dr2)/R-Bz/R
  dBzdZ= (ddpsi_dtfl2*dtfl_dr*dtfl_dz                                     &
       +  ddpsi_dtfl_dtheta*(dtfl_dr*dtheta_dz+dtfl_dz*dtheta_dr)         &
       +  ddpsi_dtheta2*dtheta_dr*dtheta_dz                               &
       +  dpsi_dtfl*ddtfl_dr_dz+dpsi_dtheta*ddtheta_dr_dz)/R
  dBzdp= (ddpsi_dtfl_dp*dtfl_dr+ddpsi_dtheta_dp*dtheta_dr)/R
!
! coordinates of the old model (QT with 0<r<1)
  rad_old=SQRT(tfl/tfl_e)
  theta_old=MODULO(theta,3.14159265358979d0*2.d0)
!
!print *,rad_old,p,theta
!print *,(Br*cos(thqt)+Bz*sin(thqt))*r/0.1d0/Bp,(Bz*cos(thqt)-Br*sin(thqt))*rho/rqt/Bp
!pause
!write (123,*) R,Z
  RETURN
END SUBROUTINE GBas_0
!  end subroutine field
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
SUBROUTINE hamilt_field(tfl_e,tfl,theta,p,dpsi_dtfl,dpsi_dtheta,        &
                          ddpsi_dtfl2,ddpsi_dtfl_dtheta,ddpsi_dtheta2,    &
                          ddpsi_dtfl_dp,ddpsi_dtheta_dp)
!
  USE vvn_period
  USE field_param
! Winny
  USE mag_interface_mod, only : aiota_tokamak
! Winny end
!
  IMPLICIT NONE
!
  DOUBLE PRECISION :: tfl_e,tfl,theta,p,dpsi_dtfl,dpsi_dtheta,        &
                      ddpsi_dtfl2,ddpsi_dtfl_dtheta,ddpsi_dtheta2,    &
                      ddpsi_dtfl_dp,ddpsi_dtheta_dp
!
  INTEGER :: m,n,iseed,nmin,nmax
  DOUBLE PRECISION :: s,aiota,daiota,sqs
  DOUBLE PRECISION :: alpha,sinalp,cosalp
  DOUBLE PRECISION :: xi,psimn,aiota_min,aiota_max
  DOUBLE PRECISION :: sin_mt_np,cos_mt_np,sin_mt_np_0,cos_mt_np_0
!
  IF(prop) THEN
    prop=.FALSE.
    OPEN(71,file='param.inp')
    READ (71,*) mmx
    READ (71,*) nmx
    READ (71,*) br0
    READ (71,*) delta_m
    READ (71,*) delta_n
    READ (71,*) r0
    READ (71,*) aiota0
    READ (71,*) aiota_pr
    READ (71,*) iseed
    READ (71,*) aiota_min
    READ (71,*) aiota_max
    READ (71,*) m_isl
    READ (71,*) n_isl
    READ (71,*) eps_isl
    CLOSE(71)
    OPEN(71,file='random.dat')
    ALLOCATE(sinmt(-mmx:mmx),cosmt(-mmx:mmx))
    ALLOCATE(sinnp(-nmx:nmx),cosnp(-nmx:nmx))
    ALLOCATE(expmn(-mmx:mmx,-nmx:nmx))
    ALLOCATE(sipsmn(-mmx:mmx,-nmx:nmx))
    ALLOCATE(copsmn(-mmx:mmx,-nmx:nmx))
!    call g05cbf(iseed)
    DO m=-mmx,mmx
      IF(m.GT.0) THEN
        nmin=MAX(-nmx,INT(aiota_min*m+1.d0))
        nmax=MIN(nmx,INT(aiota_max*m))
      ELSEIF(m.LT.0) THEN
        nmin=MAX(-nmx,INT(aiota_max*m))
        nmax=MIN(nmx,INT(aiota_min*m-1.d0))
      ELSE
        nmin=0
        nmax=0
      ENDIF
!      do n=-nmx,nmx
      DO n=nmin,nmax
!        expmn(m,n)=br0*exp(-(m/delta_m)**2-(n/delta_n)**2)/delta_m
        expmn(m,n)=br0/delta_m
!        xi=g05caf(xi)
        READ(71,*)xi
        psimn=2.d0*3.14159265358979d0*xi
        sipsmn(m,n)=SIN(psimn)
        copsmn(m,n)=COS(psimn)
      ENDDO
    ENDDO
    CLOSE(71)
  ENDIF
!
! s - dimensionless flux label (toroidal flux normalizad to 1 at the edge)
  s=tfl/tfl_e
  sqs=SQRT(s)
!    Unperturbed field:
! aiota - rotational transform angle in 2 Pi units
!  aiota=2.d0*s**2
!  daiota=4.d0*s
!!  aiota=sqs
!!  daiota=0.5d0/sqs
!Shearless field with q=3
! Winny - aiota 
  !aiota=1.d0/3.d0
  !aiota=0.35145
  aiota = aiota_tokamak
! Winny - end
  daiota=0.d0

!  dpsi_dtfl=-aiota
  dpsi_dtfl=aiota
  dpsi_dtheta=0.d0
!  ddpsi_dtfl2=-daiota/tfl_e
  ddpsi_dtfl2=daiota/tfl_e
  ddpsi_dtheta2=0.d0
  ddpsi_dtfl_dtheta=0.d0
  ddpsi_dtfl_dp=0.d0
  ddpsi_dtheta_dp=0.d0

!  dpsi_dtfl=-aiota
  dpsi_dtfl=aiota
  dpsi_dtheta=0.d0
!  ddpsi_dtfl2=-daiota/tfl_e
  ddpsi_dtfl2=daiota/tfl_e
  ddpsi_dtheta2=0.d0
  ddpsi_dtfl_dtheta=0.d0
  ddpsi_dtfl_dp=0.d0
  ddpsi_dtheta_dp=0.d0

!We don't need any perturbation here:
  RETURN


! $\psi=2\Phi_e(s^{3/2}/3-b_{r0} s^{1/2}\sum_{m,n} \sin(m\theta - n \varphi))$
!    Perturbation field:
  DO m=-mmx,mmx
    sinmt(m)=SIN(m*theta)
    cosmt(m)=COS(m*theta)
  ENDDO
  DO n=-nmx,nmx
    sinnp(n)=SIN(n*p)
    cosnp(n)=COS(n*p)
  ENDDO
!
  DO m=-mmx,mmx
!    do n=-nmx,nmx
    IF(m.GT.0) THEN
      nmin=MAX(-nmx,INT(aiota_min*m+1.d0))
      nmax=MIN(nmx,INT(aiota_max*m))
    ELSEIF(m.LT.0) THEN
      nmin=MAX(-nmx,INT(aiota_max*m))
      nmax=MIN(nmx,INT(aiota_min*m-1.d0))
    ELSE
      nmin=0
      nmax=0
    ENDIF
    DO n=nmin,nmax
      sin_mt_np_0 = sinmt(m)*cosnp(n) - cosmt(m)*sinnp(n)
      cos_mt_np_0 = cosmt(m)*cosnp(n) + sinmt(m)*sinnp(n)
      sin_mt_np = sin_mt_np_0*copsmn(m,n) + cos_mt_np_0*sipsmn(m,n)
      cos_mt_np = cos_mt_np_0*copsmn(m,n) - sin_mt_np_0*sipsmn(m,n)
!      br = br + m*cos_mt_np*expmn(m,n)
!      bt = bt + sin_mt_np*expmn(m,n)
!      dbr_dt = dbr_dt - m**2*sin_mt_np*expmn(m,n)
!      dbr_dp = dbr_dp + m*n*sin_mt_np*expmn(m,n)
!      dbt_dt = dbt_dt + m*cos_mt_np*expmn(m,n)
!      dbt_dp = dbt_dp - n*cos_mt_np*expmn(m,n)

      dpsi_dtfl=dpsi_dtfl+expmn(m,n)/sqs*sin_mt_np
      dpsi_dtheta=dpsi_dtheta+2.d0*tfl_e*expmn(m,n)*sqs*m*cos_mt_np

      ddpsi_dtfl2=ddpsi_dtfl2-0.5d0*expmn(m,n)/sqs/s/tfl_e*sin_mt_np
      ddpsi_dtheta2=ddpsi_dtheta2-2.d0*tfl_e*expmn(m,n)*sqs*m**2*sin_mt_np
      ddpsi_dtfl_dtheta=ddpsi_dtfl_dtheta+expmn(m,n)/sqs*m*cos_mt_np
      ddpsi_dtfl_dp=ddpsi_dtfl_dp-expmn(m,n)/sqs*n*cos_mt_np
      ddpsi_dtheta_dp=ddpsi_dtheta_dp+2.d0*tfl_e*expmn(m,n)*sqs*m*n*sin_mt_np
    ENDDO
  ENDDO

  ! big island:
  cos_mt_np=COS(m_isl*theta-n_isl*p)
  sin_mt_np=SIN(m_isl*theta-n_isl*p)

  m=m_isl
  n=n_isl

  dpsi_dtfl=dpsi_dtfl+eps_isl/sqs*sin_mt_np
  dpsi_dtheta=dpsi_dtheta+2.d0*tfl_e*eps_isl*sqs*m*cos_mt_np

  ddpsi_dtfl2=ddpsi_dtfl2-0.5d0*eps_isl/sqs/s/tfl_e*sin_mt_np
  ddpsi_dtheta2=ddpsi_dtheta2-2.d0*tfl_e*eps_isl*sqs*m**2*sin_mt_np
  ddpsi_dtfl_dtheta=ddpsi_dtfl_dtheta+eps_isl/sqs*m*cos_mt_np
  ddpsi_dtfl_dp=ddpsi_dtfl_dp-eps_isl/sqs*n*cos_mt_np

END SUBROUTINE hamilt_field
