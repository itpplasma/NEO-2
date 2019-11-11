MODULE plagrange_mod

  USE magnetics_mod

  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE :: dp = KIND(1.0d0)

  PUBLIC plagrange_coeff
  PRIVATE plag_coeff
  INTERFACE plagrange_coeff
     MODULE PROCEDURE plag_coeff
  END INTERFACE

  PUBLIC plagrange_test
  PRIVATE plag_test,plag_testperiod_2
  INTERFACE plagrange_test
     MODULE PROCEDURE plag_test,plag_testperiod,plag_testperiod_2
  END INTERFACE

  PUBLIC plagrange_stencel
  PRIVATE plag_stencel
  INTERFACE plagrange_stencel
     MODULE PROCEDURE plag_stencel
  END INTERFACE

  PUBLIC plagrange_value
  PRIVATE plag_value_1,plag_value_2,plag_value_all,plag_value_all2
  INTERFACE plagrange_value
     MODULE PROCEDURE plag_value_1,plag_value_2,plag_value_all,&
          plag_value_all2
  END INTERFACE

  PUBLIC plagrange_interp
  PRIVATE plag_interp_2,plag_interp_3,plag_interp_all,plag_interp_all2
  INTERFACE plagrange_interp
     MODULE PROCEDURE plag_interp_2,plag_interp_3,plag_interp_all,&
          plag_interp_all2
  END INTERFACE

CONTAINS
  !---------------------------------------------------------------------

  SUBROUTINE plag_coeff(npoi,nder,x,xp,coef)
    !
    ! npoi - number of points (determines the order of Lagrange polynomial
    ! which is equal npoi-1)
    ! nder - number of derivatives computed 0 - function only, 1 - first
    ! derivative
    ! x - actual point where function and derivatives are evaluated
    ! xp(npoi) - array of points where function is known
    ! coef(0:nder,npoi) - weights for computation of function and derivatives,
    ! f=sum(fun(1:npoi)*coef(0,1:npoi) gives the function value
    ! df=sum(fun(1:npoi)*coef(1,1:npoi) gives the derivative value value
    !
    !
    INTEGER, INTENT(in)                                :: npoi,nder
    REAL(kind=dp), INTENT(in)                          :: x
    REAL(kind=dp), DIMENSION(npoi), INTENT(in)         :: xp
    REAL(kind=dp), DIMENSION(0:nder,npoi), INTENT(out) :: coef
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: dummy
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: fak_i
    !
    INTEGER                                            :: i,k,j,l
    REAL(kind=dp)                                      :: fac
    REAL(kind=dp)                                      :: j_sum,l_sum,k_prod
    !
    DO i=1,npoi
       coef(0,i)=1.d0
       DO k=1,npoi
          IF(k.EQ.i) CYCLE
          coef(0,i)=coef(0,i)*(x-xp(k))/(xp(i)-xp(k))
       ENDDO
    ENDDO
    !
    IF(nder.EQ.0) RETURN
    !
    ALLOCATE(dummy(npoi))
    !
    DO i=1,npoi
       dummy=1.d0
       dummy(i)=0.d0
       DO k=1,npoi
          IF(k.EQ.i) CYCLE
          fac=(x-xp(k))/(xp(i)-xp(k))
          DO j=1,npoi
             IF(j.EQ.k) THEN
                dummy(j)=dummy(j)/(xp(i)-xp(k))
             ELSE
                dummy(j)=dummy(j)*fac
             ENDIF
          ENDDO
       ENDDO
       coef(1,i)=SUM(dummy)
    ENDDO
    !
    DEALLOCATE(dummy)
    !
    IF(nder.LE.1) RETURN
    !
    ALLOCATE(fak_i(npoi))
    do_i: DO i = 1,npoi
       fak_i = 0.0d0
       do_prep: DO k = 1,npoi
          IF (k .EQ. i) CYCLE
          fak_i(k) = (x-xp(k)) / (xp(i)-xp(k))
       END DO do_prep
       j_sum = 0.0d0
       do_j: DO j =1,npoi
          IF (j .EQ. i) CYCLE
          l_sum = 0.0d0
          do_l: DO l = 1,npoi
             IF (l .EQ. i .OR. l .EQ. j) CYCLE
             k_prod = 1.0d0
             do_k: DO k =1,npoi
                IF (k .EQ. i .OR. k .EQ. j .OR. k .EQ. l) CYCLE
                k_prod = k_prod * fak_i(k)
             END DO do_k
             l_sum = l_sum + k_prod / (xp(i)-xp(l))
          END DO do_l
          j_sum = j_sum + l_sum / (xp(i)-xp(j))
       END DO do_j
       coef(2,i)=j_sum
    END DO do_i
    DEALLOCATE(fak_i)

    RETURN
  END SUBROUTINE plag_coeff
  !--------------------------------------------------------------------
  SUBROUTINE plag_test(npoi,imax)
    INTEGER, INTENT(in) :: npoi
    INTEGER, INTENT(in) :: imax

    INTEGER, PARAMETER                         :: nder=1
    INTEGER                                    :: i
    REAL(kind=dp)                              :: u,umax
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE   :: up,fun
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: coeff
    !
    umax=3.1415926d0
    !
    ALLOCATE( up(npoi) )
    ALLOCATE( fun(npoi) )
    ALLOCATE( coeff(0:nder,npoi) )

    OPEN(1,file='input_plag.dat')
    DO i=1,npoi
       up(i)=umax*(float(i-1)/float(npoi-1))
       fun(i)=SIN(up(i))
       WRITE (1,*) up(i),fun(i)
    ENDDO
    CLOSE(1)
    !
    OPEN(1,file='sinus.dat')
    DO i=1,imax
       u=umax*(float(i-1)/float(imax-1))
       CALL plagrange_coeff(npoi,nder,u,up(:),coeff(:,:))
       WRITE (1,*) u,SIN(u),SUM(coeff(0,:)*fun),COS(u),SUM(coeff(1,:)*fun)
    ENDDO
    CLOSE(1)
    !
  END SUBROUTINE plag_test

  !---------------------------------------------------------------------
  !> \brief Output (to file) some (interpolated?) quantities for test purposes.
  SUBROUTINE plag_testperiod(fieldperiod,nlagrange,ndata)
    TYPE(fieldperiod_struct), POINTER :: fieldperiod
    INTEGER, INTENT(in)      :: nlagrange,ndata

    INTEGER, PARAMETER                         :: nder=1
    INTEGER :: nstep,i
    INTEGER :: u1 = 117

    REAL(kind=dp) :: phi_start,phi_end,phi_span,phi
    REAL(kind=dp) :: bhat,bhatder
    REAL(kind=dp) :: x1,x3,geodcu,h_phi,dlogbdphi
    REAL(kind=dp) :: dbcovar_s_hat_dphi,bcovar_s_hat,dlogbds

    nstep = fieldperiod%parent%parent%nstep
    phi_start = fieldperiod%coords%x2(0)
    phi_end   = fieldperiod%coords%x2(nstep)
    phi_span  = phi_end - phi_start

    OPEN(unit=u1,file='testperiod.dat')
    DO i = 0, nstep
       !! Modifications by Andreas F. Martitsch (11.06.2014)
       ! Optional output (necessary for modeling the magnetic rotation)
       IF ( ALLOCATED(fieldperiod%mdata%dbcovar_s_hat_dphi) .AND. &
            ALLOCATED(fieldperiod%mdata%bcovar_s_hat)       .AND. &
            ALLOCATED(fieldperiod%mdata%dlogbds) ) THEN
          WRITE(u1,*) &
               fieldperiod%coords%x1(i),fieldperiod%coords%x2(i), &
               fieldperiod%coords%x3(i), &
               fieldperiod%mdata%bhat(i),fieldperiod%mdata%geodcu(i), &
               fieldperiod%mdata%h_phi(i),fieldperiod%mdata%dlogbdphi(i), &
               fieldperiod%mdata%dbcovar_s_hat_dphi(i), &
               fieldperiod%mdata%bcovar_s_hat(i),fieldperiod%mdata%dlogbds(i)
       ELSE ! This is the old version:
          WRITE(u1,*) &
               fieldperiod%coords%x1(i),fieldperiod%coords%x2(i),&
               fieldperiod%coords%x3(i), &
               fieldperiod%mdata%bhat(i),fieldperiod%mdata%geodcu(i), &
               fieldperiod%mdata%h_phi(i),fieldperiod%mdata%dlogbdphi(i)
       END IF
       !! End Modifications by Andreas F. Martitsch (11.06.2014)
    END DO
    CLOSE(unit=u1)

    OPEN(unit=u1,file='testlagrange.dat')
    DO i = 0, ndata
       phi = phi_start + phi_span * DBLE(i) / DBLE(ndata)
       CALL plagrange_interp(fieldperiod,phi,nlagrange,bhat,bhatder)
       WRITE(u1,*) phi,bhat,bhatder
    END DO
    CLOSE(unit=u1)
    !PRINT *, 'Derivative'
    !PAUSE

    OPEN(unit=u1,file='testlagall.dat')
    DO i = 0, ndata
       phi = phi_start + phi_span * DBLE(i) / DBLE(ndata)
       !! Modifications by Andreas F. Martitsch (11.06.2014)
       ! Optional output (necessary for modeling the magnetic rotation)
       IF ( ALLOCATED(fieldperiod%mdata%dbcovar_s_hat_dphi) .AND. &
            ALLOCATED(fieldperiod%mdata%bcovar_s_hat)       .AND. &
            ALLOCATED(fieldperiod%mdata%dlogbds) ) THEN
          CALL plagrange_interp(fieldperiod,phi,nlagrange,x1,x3,bhat,geodcu,h_phi,&
               dlogbdphi,dbcovar_s_hat_dphi,bcovar_s_hat,dlogbds)
          WRITE(u1,*) x1,phi,x3,bhat,geodcu,h_phi,dlogbdphi,&
               dbcovar_s_hat_dphi,bcovar_s_hat,dlogbds
       ELSE ! This is the old version:
          CALL plagrange_interp(fieldperiod,phi,nlagrange,x1,x3,bhat,geodcu,h_phi,dlogbdphi)
          WRITE(u1,*) x1,phi,x3,bhat,geodcu,h_phi,dlogbdphi
       END IF
       !! End Modifications by Andreas F. Martitsch (11.06.2014)
    END DO
    CLOSE(unit=u1)




    PRINT *, '------ plag_testperiod ------'
    PRINT *, 'fieldperiod%tag', fieldperiod%tag
    PRINT *, 'nstep ',nstep

  END SUBROUTINE plag_testperiod
  !---------------------------------------------------------------------
  SUBROUTINE plag_testperiod_2(fieldperiod)
    TYPE(fieldperiod_struct), POINTER :: fieldperiod
    TYPE(fieldpropagator_struct), POINTER :: fieldpropagator
    TYPE(fieldripple_struct), POINTER :: fieldripple

    INTEGER :: u1 = 117
    INTEGER :: u2 = 118
    INTEGER :: u3 = 119
    INTEGER :: i,count
    INTEGER :: pend_tag
    INTEGER :: pstart_tag

    REAL(kind=dp) :: r0

    pend_tag = 15
    pstart_tag = -1

    fieldperiod => fieldperiod%parent%ch_fir
    fieldpropagator => fieldperiod%ch_ext
    fieldripple => fieldpropagator%ch_act

    r0 = fieldperiod%parent%parent%parent%r0


    OPEN(unit=u1,file='periods.dat')
    OPEN(unit=u2,file='limits.dat')
    count = 0
    DO
       IF (fieldperiod%tag .LT. pstart_tag) THEN
          IF (ASSOCIATED(fieldperiod%next)) THEN
             fieldperiod => fieldperiod%next
          ELSE
             EXIT
          END IF
          CYCLE
       END IF
       DO i = count, UBOUND(fieldperiod%coords%x2,1)
          WRITE(u1,'(1000f25.16)') &
               fieldperiod%coords%x1(i),fieldperiod%coords%x2(i),fieldperiod%coords%x3(i), &
               fieldperiod%mdata%bhat(i),fieldperiod%mdata%geodcu(i), &
               fieldperiod%mdata%h_phi(i),fieldperiod%mdata%dlogbdphi(i)
       END DO
       count = 1
       WRITE(u2,'(1000f25.16)') fieldperiod%phi_l
       WRITE(u2,'(1000f25.16)') fieldperiod%phi_r
       IF (ASSOCIATED(fieldperiod%next)) THEN
          fieldperiod => fieldperiod%next
       ELSE
          EXIT
       END IF
       IF (fieldperiod%tag .GE. pend_tag) EXIT
    END DO
    CLOSE(unit=u1)
    CLOSE(unit=u2)

    OPEN(unit=u1,file='props.dat')
    OPEN(unit=u2,file='extrema.dat')
    OPEN(unit=u3,file='proplim.dat')
    count = 0
    DO
       IF (fieldpropagator%parent%tag .LT. pstart_tag) THEN
          IF (ASSOCIATED(fieldpropagator%next)) THEN
             fieldpropagator => fieldpropagator%next
          ELSE
             EXIT
          END IF
          CYCLE
       END IF
       IF (fieldpropagator%parent%tag .GT. pend_tag) EXIT
       DO i = count, UBOUND(fieldpropagator%coords%x2,1)
          WRITE(u1,'(1000f25.16)') &
               fieldpropagator%coords%x1(i),fieldpropagator%coords%x2(i),fieldpropagator%coords%x3(i), &
               fieldpropagator%mdata%bhat(i),fieldpropagator%mdata%geodcu(i), &
               fieldpropagator%mdata%h_phi(i),fieldpropagator%mdata%dlogbdphi(i)
       END DO
       count = 1
       WRITE(u2,'(1000f25.16)') fieldpropagator%phi_l,fieldpropagator%b_l
       WRITE(u2,'(1000f25.16)') fieldpropagator%phi_min,fieldpropagator%b_min
       WRITE(u2,'(1000f25.16)') fieldpropagator%phi_r,fieldpropagator%b_r
       WRITE(u3,'(1000f25.16)') fieldpropagator%phi_l
       WRITE(u3,'(1000f25.16)') fieldpropagator%phi_r
       IF (ASSOCIATED(fieldpropagator%next)) THEN
          fieldpropagator => fieldpropagator%next
       ELSE
          EXIT
       END IF
    END DO
    CLOSE(unit=u1)
    CLOSE(unit=u2)
    CLOSE(unit=u3)

    OPEN(unit=u1,file='ripple.dat')
    DO
       IF (fieldripple%pa_fir%parent%tag .LT. pstart_tag) THEN
          IF (ASSOCIATED(fieldripple%next)) THEN
             fieldripple => fieldripple%next
          ELSE
             EXIT
          END IF
          CYCLE
       END IF
       IF (fieldripple%pa_fir%parent%tag .GT. pend_tag) EXIT
       WRITE(u1,'(1000f25.16)') fieldripple%pa_fir%phi_l,fieldripple%b_max_l
       WRITE(u1,'(1000f25.16)') fieldripple%pa_fir%phi_l + fieldripple%width_l / r0 ,fieldripple%b_min
       WRITE(u1,'(1000f25.16)') fieldripple%pa_las%phi_r,fieldripple%b_max_r
       IF (ASSOCIATED(fieldripple%next)) THEN
          fieldripple => fieldripple%next
       ELSE
          EXIT
       END IF
    END DO
    CLOSE(unit=u1)

  END SUBROUTINE plag_testperiod_2

  !---------------------------------------------------------------------
  SUBROUTINE plag_stencel(fieldperiod,phi,nlagrange,stencel)
    TYPE(fieldperiod_struct), POINTER                 :: fieldperiod, fp
    REAL(kind=dp), INTENT(in)                         :: phi
    INTEGER, INTENT(in)                               :: nlagrange
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(inout) :: stencel

    REAL(kind=dp) :: phi_start,phi_end,phi_span
    INTEGER       :: nstep,mid,i


    !INTEGER :: nstep
    REAL(kind=dp) ::  phi_l
    REAL(kind=dp) ::  phi_last,phi_first

    fp => fieldperiod

    nstep = UBOUND(fp%coords%x2,1)
    phi_l = phi
    ! fix phi_l locally to be within range
    phi_first = fp%parent%ch_fir%coords%x2(0)
    phi_last  = fp%parent%ch_las%coords%x2(nstep)
    phi_span  = phi_last - phi_first
    DO WHILE (phi_l .GT. phi_last)
       phi_l = phi_l - phi_span
    END DO
    DO WHILE (phi_l .LT. phi_first)
       phi_l = phi_l + phi_span
    END DO
    ! fix period if necessary
    phi_last  = fp%coords%x2(nstep)
    DO WHILE (phi_l .GT. phi_last)
       fp => fp%next
       phi_last  = fp%coords%x2(nstep)
    END DO
    phi_first = fp%coords%x2(0)
    DO WHILE (phi_l .LT. phi_first)
       fp => fp%prev
       phi_first  = fp%coords%x2(0)
    END DO

    nstep = fp%parent%parent%nstep
    phi_start = fp%coords%x2(0)
    phi_end   = fp%coords%x2(nstep)
    phi_span  = phi_end - phi_start

    IF (ALLOCATED(stencel)) DEALLOCATE(stencel)
    ALLOCATE(stencel(nlagrange+1))

    mid = INT( (phi_l - phi_start) / phi_span * DBLE(nstep) )
    DO i = 1,nlagrange+1
       stencel(i) = mid - (nlagrange+1)/2 + i
    END DO

  END SUBROUTINE plag_stencel
  !---------------------------------------------------------------------
  SUBROUTINE plag_value_1(fieldperiod,stencel,phi_arr)
    TYPE(fieldperiod_struct), POINTER                       :: fieldperiod
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(inout)       :: stencel
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: phi_arr

    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: hlp

    INTEGER :: nstep

    IF (ALLOCATED(phi_arr)) DEALLOCATE(phi_arr)
    ALLOCATE(phi_arr(LBOUND(stencel,1):UBOUND(stencel,1)))

    nstep = fieldperiod%parent%parent%nstep

    IF (MINVAL(stencel) .LT. 0) THEN
       ALLOCATE(hlp(-nstep:nstep))
       IF (ASSOCIATED(fieldperiod%prev)) THEN
          hlp(-nstep:0) = fieldperiod%prev%coords%x2(0:nstep)
          hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
       ELSE
          hlp(-nstep:0) = fieldperiod%parent%ch_las%coords%x2(0:nstep) &
               - fieldperiod%parent%ch_las%coords%x2(nstep) &
               + fieldperiod%coords%x2(0)
          hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
       END IF
    ELSEIF (MAXVAL(stencel) .GT. nstep) THEN
       ALLOCATE(hlp(0:2*nstep))
       IF (ASSOCIATED(fieldperiod%next)) THEN
          hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%next%coords%x2(0:nstep)
       ELSE
          hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%coords%x2(0:nstep) &
               - fieldperiod%parent%ch_fir%coords%x2(0) &
               + fieldperiod%coords%x2(nstep)
       END IF
    ELSE
       ALLOCATE(hlp(0:nstep))
       hlp(0:nstep) = fieldperiod%coords%x2(0:nstep)
    END IF
    phi_arr = hlp(stencel)
    DEALLOCATE(hlp)

  END SUBROUTINE plag_value_1
  !---------------------------------------------------------------------
  SUBROUTINE plag_value_2(fieldperiod,stencel,phi_arr,bhat_arr)
    TYPE(fieldperiod_struct), POINTER                       :: fieldperiod
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(inout)       :: stencel
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: phi_arr
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: bhat_arr

    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: hlp

    INTEGER :: nstep

    nstep = fieldperiod%parent%parent%nstep

    ! phi_arr
    IF (ALLOCATED(phi_arr)) DEALLOCATE(phi_arr)
    ALLOCATE(phi_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    IF (MINVAL(stencel) .LT. 0) THEN
       ALLOCATE(hlp(-nstep:nstep))
       IF (ASSOCIATED(fieldperiod%prev)) THEN
          hlp(-nstep:0) = fieldperiod%prev%coords%x2(0:nstep)
          hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
       ELSE
          hlp(-nstep:0) = fieldperiod%parent%ch_las%coords%x2(0:nstep) &
               - fieldperiod%parent%ch_las%coords%x2(nstep) &
               + fieldperiod%coords%x2(0)
          hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
       END IF
    ELSEIF (MAXVAL(stencel) .GT. nstep) THEN
       ALLOCATE(hlp(0:2*nstep))
       IF (ASSOCIATED(fieldperiod%next)) THEN
          hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%next%coords%x2(0:nstep)
       ELSE
          hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%coords%x2(0:nstep) &
               - fieldperiod%parent%ch_fir%coords%x2(0) &
               + fieldperiod%coords%x2(nstep)
       END IF
    ELSE
       ALLOCATE(hlp(0:nstep))
       hlp(0:nstep) = fieldperiod%coords%x2(0:nstep)
    END IF
    phi_arr = hlp(stencel)
    DEALLOCATE(hlp)

    ! bhat_arr
    IF (ALLOCATED(bhat_arr)) DEALLOCATE(bhat_arr)
    ALLOCATE(bhat_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    IF (MINVAL(stencel) .LT. 0) THEN
       ALLOCATE(hlp(-nstep:nstep))
       IF (ASSOCIATED(fieldperiod%prev)) THEN
          hlp(-nstep:0) = fieldperiod%prev%mdata%bhat(0:nstep)
          hlp(0:nstep)  = fieldperiod%mdata%bhat(0:nstep)
       ELSE
          hlp(-nstep:0) = fieldperiod%parent%ch_las%mdata%bhat(0:nstep)
          hlp(0:nstep)  = fieldperiod%mdata%bhat(0:nstep)
       END IF
    ELSEIF (MAXVAL(stencel) .GT. nstep) THEN
       ALLOCATE(hlp(0:2*nstep))
       IF (ASSOCIATED(fieldperiod%next)) THEN
          hlp(0:nstep)  = fieldperiod%mdata%bhat(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%next%mdata%bhat(0:nstep)
       ELSE
          hlp(0:nstep)  = fieldperiod%mdata%bhat(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%mdata%bhat(0:nstep)
       END IF
    ELSE
       ALLOCATE(hlp(0:nstep))
       hlp(0:nstep) = fieldperiod%mdata%bhat(0:nstep)
    END IF
    bhat_arr = hlp(stencel)
    DEALLOCATE(hlp)
    ! phi_arr, coords, x2
  END SUBROUTINE plag_value_2
  !---------------------------------------------------------------------
  SUBROUTINE plag_value_all(fieldperiod,stencel,phi_arr, &
       x1_arr,x3_arr,bhat_arr,geodcu_arr,h_phi_arr,dlogbdphi_arr)
    TYPE(fieldperiod_struct), POINTER                       :: fieldperiod
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(inout)       :: stencel
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: phi_arr
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: x1_arr,x3_arr,bhat_arr
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: geodcu_arr,h_phi_arr,dlogbdphi_arr

    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: hlp

    INTEGER :: nstep

    nstep = fieldperiod%parent%parent%nstep

    ! phi_arr
    IF (ALLOCATED(phi_arr)) DEALLOCATE(phi_arr)
    ALLOCATE(phi_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    IF (MINVAL(stencel) .LT. 0) THEN
       ALLOCATE(hlp(-nstep:nstep))
       IF (ASSOCIATED(fieldperiod%prev)) THEN
          hlp(-nstep:0) = fieldperiod%prev%coords%x2(0:nstep)
          hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
       ELSE
          hlp(-nstep:0) = fieldperiod%parent%ch_las%coords%x2(0:nstep) &
               - fieldperiod%parent%ch_las%coords%x2(nstep) &
               + fieldperiod%coords%x2(0)
          hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
       END IF
    ELSEIF (MAXVAL(stencel) .GT. nstep) THEN
       ALLOCATE(hlp(0:2*nstep))
       IF (ASSOCIATED(fieldperiod%next)) THEN
          hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%next%coords%x2(0:nstep)
       ELSE
          hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%coords%x2(0:nstep) &
               - fieldperiod%parent%ch_fir%coords%x2(0) &
               + fieldperiod%coords%x2(nstep)
       END IF
    ELSE
       ALLOCATE(hlp(0:nstep))
       hlp(0:nstep) = fieldperiod%coords%x2(0:nstep)
    END IF
    phi_arr = hlp(stencel)
    DEALLOCATE(hlp)

    ! x1_arr
    IF (ALLOCATED(x1_arr)) DEALLOCATE(x1_arr)
    ALLOCATE(x1_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    IF (MINVAL(stencel) .LT. 0) THEN
       ALLOCATE(hlp(-nstep:nstep))
       IF (ASSOCIATED(fieldperiod%prev)) THEN
          hlp(-nstep:0) = fieldperiod%prev%coords%x1(0:nstep)
          hlp(0:nstep)  = fieldperiod%coords%x1(0:nstep)
       ELSE
          hlp(-nstep:0) = fieldperiod%parent%ch_las%coords%x1(0:nstep)
          hlp(0:nstep)  = fieldperiod%coords%x1(0:nstep)
       END IF
    ELSEIF (MAXVAL(stencel) .GT. nstep) THEN
       ALLOCATE(hlp(0:2*nstep))
       IF (ASSOCIATED(fieldperiod%next)) THEN
          hlp(0:nstep)  = fieldperiod%coords%x1(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%next%coords%x1(0:nstep)
       ELSE
          hlp(0:nstep)  = fieldperiod%coords%x1(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%coords%x1(0:nstep)
       END IF
    ELSE
       ALLOCATE(hlp(0:nstep))
       hlp(0:nstep) = fieldperiod%coords%x1(0:nstep)
    END IF
    x1_arr = hlp(stencel)
    DEALLOCATE(hlp)
    ! x1

    ! x3_arr
    IF (ALLOCATED(x3_arr)) DEALLOCATE(x3_arr)
    ALLOCATE(x3_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    IF (MINVAL(stencel) .LT. 0) THEN
       ALLOCATE(hlp(-nstep:nstep))
       IF (ASSOCIATED(fieldperiod%prev)) THEN
          hlp(-nstep:0) = fieldperiod%prev%coords%x3(0:nstep)
          hlp(0:nstep)  = fieldperiod%coords%x3(0:nstep)
       ELSE
          hlp(-nstep:0) = fieldperiod%parent%ch_las%coords%x3(0:nstep)
          hlp(0:nstep)  = fieldperiod%coords%x3(0:nstep)
       END IF
    ELSEIF (MAXVAL(stencel) .GT. nstep) THEN
       ALLOCATE(hlp(0:2*nstep))
       IF (ASSOCIATED(fieldperiod%next)) THEN
          hlp(0:nstep)  = fieldperiod%coords%x3(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%next%coords%x3(0:nstep)
       ELSE
          hlp(0:nstep)  = fieldperiod%coords%x3(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%coords%x3(0:nstep)
       END IF
    ELSE
       ALLOCATE(hlp(0:nstep))
       hlp(0:nstep) = fieldperiod%coords%x3(0:nstep)
    END IF
    x3_arr = hlp(stencel)
    DEALLOCATE(hlp)
    ! x3

    ! bhat_arr
    IF (ALLOCATED(bhat_arr)) DEALLOCATE(bhat_arr)
    ALLOCATE(bhat_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    IF (MINVAL(stencel) .LT. 0) THEN
       ALLOCATE(hlp(-nstep:nstep))
       IF (ASSOCIATED(fieldperiod%prev)) THEN
          hlp(-nstep:0) = fieldperiod%prev%mdata%bhat(0:nstep)
          hlp(0:nstep)  = fieldperiod%mdata%bhat(0:nstep)
       ELSE
          hlp(-nstep:0) = fieldperiod%parent%ch_las%mdata%bhat(0:nstep)
          hlp(0:nstep)  = fieldperiod%mdata%bhat(0:nstep)
       END IF
    ELSEIF (MAXVAL(stencel) .GT. nstep) THEN
       ALLOCATE(hlp(0:2*nstep))
       IF (ASSOCIATED(fieldperiod%next)) THEN
          hlp(0:nstep)  = fieldperiod%mdata%bhat(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%next%mdata%bhat(0:nstep)
       ELSE
          hlp(0:nstep)  = fieldperiod%mdata%bhat(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%mdata%bhat(0:nstep)
       END IF
    ELSE
       ALLOCATE(hlp(0:nstep))
       hlp(0:nstep) = fieldperiod%mdata%bhat(0:nstep)
    END IF
    bhat_arr = hlp(stencel)
    DEALLOCATE(hlp)
    ! bhat_arr

    ! geodcu_arr
    IF (ALLOCATED(geodcu_arr)) DEALLOCATE(geodcu_arr)
    ALLOCATE(geodcu_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    IF (MINVAL(stencel) .LT. 0) THEN
       ALLOCATE(hlp(-nstep:nstep))
       IF (ASSOCIATED(fieldperiod%prev)) THEN
          hlp(-nstep:0) = fieldperiod%prev%mdata%geodcu(0:nstep)
          hlp(0:nstep)  = fieldperiod%mdata%geodcu(0:nstep)
       ELSE
          hlp(-nstep:0) = fieldperiod%parent%ch_las%mdata%geodcu(0:nstep)
          hlp(0:nstep)  = fieldperiod%mdata%geodcu(0:nstep)
       END IF
    ELSEIF (MAXVAL(stencel) .GT. nstep) THEN
       ALLOCATE(hlp(0:2*nstep))
       IF (ASSOCIATED(fieldperiod%next)) THEN
          hlp(0:nstep)  = fieldperiod%mdata%geodcu(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%next%mdata%geodcu(0:nstep)
       ELSE
          hlp(0:nstep)  = fieldperiod%mdata%geodcu(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%mdata%geodcu(0:nstep)
       END IF
    ELSE
       ALLOCATE(hlp(0:nstep))
       hlp(0:nstep) = fieldperiod%mdata%geodcu(0:nstep)
    END IF
    geodcu_arr = hlp(stencel)
    DEALLOCATE(hlp)
    ! geodcu_arr

    ! h_phi_arr
    IF (ALLOCATED(h_phi_arr)) DEALLOCATE(h_phi_arr)
    ALLOCATE(h_phi_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    IF (MINVAL(stencel) .LT. 0) THEN
       ALLOCATE(hlp(-nstep:nstep))
       IF (ASSOCIATED(fieldperiod%prev)) THEN
          hlp(-nstep:0) = fieldperiod%prev%mdata%h_phi(0:nstep)
          hlp(0:nstep)  = fieldperiod%mdata%h_phi(0:nstep)
       ELSE
          hlp(-nstep:0) = fieldperiod%parent%ch_las%mdata%h_phi(0:nstep)
          hlp(0:nstep)  = fieldperiod%mdata%h_phi(0:nstep)
       END IF
    ELSEIF (MAXVAL(stencel) .GT. nstep) THEN
       ALLOCATE(hlp(0:2*nstep))
       IF (ASSOCIATED(fieldperiod%next)) THEN
          hlp(0:nstep)  = fieldperiod%mdata%h_phi(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%next%mdata%h_phi(0:nstep)
       ELSE
          hlp(0:nstep)  = fieldperiod%mdata%h_phi(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%mdata%h_phi(0:nstep)
       END IF
    ELSE
       ALLOCATE(hlp(0:nstep))
       hlp(0:nstep) = fieldperiod%mdata%h_phi(0:nstep)
    END IF
    h_phi_arr = hlp(stencel)
    DEALLOCATE(hlp)
    ! h_phi_arr

    ! dlogbdphi_arr
    IF (ALLOCATED(dlogbdphi_arr)) DEALLOCATE(dlogbdphi_arr)
    ALLOCATE(dlogbdphi_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    IF (MINVAL(stencel) .LT. 0) THEN
       ALLOCATE(hlp(-nstep:nstep))
       IF (ASSOCIATED(fieldperiod%prev)) THEN
          hlp(-nstep:0) = fieldperiod%prev%mdata%dlogbdphi(0:nstep)
          hlp(0:nstep)  = fieldperiod%mdata%dlogbdphi(0:nstep)
       ELSE
          hlp(-nstep:0) = fieldperiod%parent%ch_las%mdata%dlogbdphi(0:nstep)
          hlp(0:nstep)  = fieldperiod%mdata%dlogbdphi(0:nstep)
       END IF
    ELSEIF (MAXVAL(stencel) .GT. nstep) THEN
       ALLOCATE(hlp(0:2*nstep))
       IF (ASSOCIATED(fieldperiod%next)) THEN
          hlp(0:nstep)  = fieldperiod%mdata%dlogbdphi(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%next%mdata%dlogbdphi(0:nstep)
       ELSE
          hlp(0:nstep)  = fieldperiod%mdata%dlogbdphi(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%mdata%dlogbdphi(0:nstep)
       END IF
    ELSE
       ALLOCATE(hlp(0:nstep))
       hlp(0:nstep) = fieldperiod%mdata%dlogbdphi(0:nstep)
    END IF
    dlogbdphi_arr = hlp(stencel)
    DEALLOCATE(hlp)
    ! dlogbdphi_arr

  END SUBROUTINE plag_value_all
  !---------------------------------------------------------------------
  !! Modifications by Andreas F. Martitsch (13.03.2014)
  ! Optional output (necessary for modeling the magnetic rotation)
  SUBROUTINE plag_value_all2(fieldperiod,stencel,phi_arr, &
       x1_arr,x3_arr,bhat_arr,geodcu_arr,h_phi_arr,dlogbdphi_arr,&
       dbcovar_s_hat_dphi_arr,bcovar_s_hat_arr,dlogbds_arr)
    TYPE(fieldperiod_struct), POINTER                       :: fieldperiod
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(inout)       :: stencel
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: phi_arr
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: x1_arr,x3_arr,bhat_arr
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: geodcu_arr,h_phi_arr
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: dlogbdphi_arr
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: dbcovar_s_hat_dphi_arr
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: bcovar_s_hat_arr
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: dlogbds_arr

    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: hlp

    INTEGER :: nstep

    nstep = fieldperiod%parent%parent%nstep

    ! Compute phi_arr, x1_arr, x3_arr, bhat_arr, geodcu_arr, h_phi_arr, dlogbdphi_arr,
    ! dbcovar_s_hat_dphi_arr, bcovar_s_hat_arr, dlogbds_arr
    CALL plag_value_all(fieldperiod,stencel,phi_arr, &
         x1_arr,x3_arr,bhat_arr,geodcu_arr,h_phi_arr,dlogbdphi_arr)

    ! dbcovar_s_hat_dphi_arr
    IF (ALLOCATED(dbcovar_s_hat_dphi_arr)) DEALLOCATE(dbcovar_s_hat_dphi_arr)
    ALLOCATE(dbcovar_s_hat_dphi_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    IF (MINVAL(stencel) .LT. 0) THEN
       ALLOCATE(hlp(-nstep:nstep))
       IF (ASSOCIATED(fieldperiod%prev)) THEN
          hlp(-nstep:0) = fieldperiod%prev%mdata%dbcovar_s_hat_dphi(0:nstep)
          hlp(0:nstep)  = fieldperiod%mdata%dbcovar_s_hat_dphi(0:nstep)
       ELSE
          hlp(-nstep:0) = fieldperiod%parent%ch_las%mdata%dbcovar_s_hat_dphi(0:nstep)
          hlp(0:nstep)  = fieldperiod%mdata%dbcovar_s_hat_dphi(0:nstep)
       END IF
    ELSEIF (MAXVAL(stencel) .GT. nstep) THEN
       ALLOCATE(hlp(0:2*nstep))
       IF (ASSOCIATED(fieldperiod%next)) THEN
          hlp(0:nstep)  = fieldperiod%mdata%dbcovar_s_hat_dphi(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%next%mdata%dbcovar_s_hat_dphi(0:nstep)
       ELSE
          hlp(0:nstep)  = fieldperiod%mdata%dbcovar_s_hat_dphi(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%mdata%dbcovar_s_hat_dphi(0:nstep)
       END IF
    ELSE
       ALLOCATE(hlp(0:nstep))
       hlp(0:nstep) = fieldperiod%mdata%dbcovar_s_hat_dphi(0:nstep)
    END IF
    dbcovar_s_hat_dphi_arr = hlp(stencel)
    DEALLOCATE(hlp)
    ! dbcovar_s_hat_dphi_arr

    ! bcovar_s_hat_arr
    IF (ALLOCATED(bcovar_s_hat_arr)) DEALLOCATE(bcovar_s_hat_arr)
    ALLOCATE(bcovar_s_hat_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    IF (MINVAL(stencel) .LT. 0) THEN
       ALLOCATE(hlp(-nstep:nstep))
       IF (ASSOCIATED(fieldperiod%prev)) THEN
          hlp(-nstep:0) = fieldperiod%prev%mdata%bcovar_s_hat(0:nstep)
          hlp(0:nstep)  = fieldperiod%mdata%bcovar_s_hat(0:nstep)
       ELSE
          hlp(-nstep:0) = fieldperiod%parent%ch_las%mdata%bcovar_s_hat(0:nstep)
          hlp(0:nstep)  = fieldperiod%mdata%bcovar_s_hat(0:nstep)
       END IF
    ELSEIF (MAXVAL(stencel) .GT. nstep) THEN
       ALLOCATE(hlp(0:2*nstep))
       IF (ASSOCIATED(fieldperiod%next)) THEN
          hlp(0:nstep)  = fieldperiod%mdata%bcovar_s_hat(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%next%mdata%bcovar_s_hat(0:nstep)
       ELSE
          hlp(0:nstep)  = fieldperiod%mdata%bcovar_s_hat(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%mdata%bcovar_s_hat(0:nstep)
       END IF
    ELSE
       ALLOCATE(hlp(0:nstep))
       hlp(0:nstep) = fieldperiod%mdata%bcovar_s_hat(0:nstep)
    END IF
    bcovar_s_hat_arr = hlp(stencel)
    DEALLOCATE(hlp)
    ! bcovar_s_hat_arr

    ! dlogbds_arr
    IF (ALLOCATED(dlogbds_arr)) DEALLOCATE(dlogbds_arr)
    ALLOCATE(dlogbds_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    IF (MINVAL(stencel) .LT. 0) THEN
       ALLOCATE(hlp(-nstep:nstep))
       IF (ASSOCIATED(fieldperiod%prev)) THEN
          hlp(-nstep:0) = fieldperiod%prev%mdata%dlogbds(0:nstep)
          hlp(0:nstep)  = fieldperiod%mdata%dlogbds(0:nstep)
       ELSE
          hlp(-nstep:0) = fieldperiod%parent%ch_las%mdata%dlogbds(0:nstep)
          hlp(0:nstep)  = fieldperiod%mdata%dlogbds(0:nstep)
       END IF
    ELSEIF (MAXVAL(stencel) .GT. nstep) THEN
       ALLOCATE(hlp(0:2*nstep))
       IF (ASSOCIATED(fieldperiod%next)) THEN
          hlp(0:nstep)  = fieldperiod%mdata%dlogbds(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%next%mdata%dlogbds(0:nstep)
       ELSE
          hlp(0:nstep)  = fieldperiod%mdata%dlogbds(0:nstep)
          hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%mdata%dlogbds(0:nstep)
       END IF
    ELSE
       ALLOCATE(hlp(0:nstep))
       hlp(0:nstep) = fieldperiod%mdata%dlogbds(0:nstep)
    END IF
    dlogbds_arr = hlp(stencel)
    DEALLOCATE(hlp)
    ! dlogbds_arr

  END SUBROUTINE plag_value_all2
  !! End Modifications by Andreas F. Martitsch (13.03.2014)
  !---------------------------------------------------------------------
  SUBROUTINE plag_interp_2(fieldperiod,phi,nlagrange,bhat,bhatder)
    TYPE(fieldperiod_struct), POINTER                 :: fieldperiod, fp
    REAL(kind=dp), INTENT(in)                         :: phi
    INTEGER, INTENT(in)                               :: nlagrange
    REAL(kind=dp), INTENT(out)                        :: bhat,bhatder

    INTEGER, PARAMETER :: nder = 1
    INTEGER, DIMENSION(:), ALLOCATABLE :: stencel
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: phi_arr,bhat_arr
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: coeff


    INTEGER :: nstep
    REAL(kind=dp) ::  phi_l
    REAL(kind=dp) ::  phi_span,phi_last,phi_first

    fp => fieldperiod

    nstep = UBOUND(fp%coords%x2,1)
    phi_l = phi
    ! fix phi_l locally to be within range
    phi_first = fp%parent%ch_fir%coords%x2(0)
    phi_last  = fp%parent%ch_las%coords%x2(nstep)
    phi_span  = phi_last - phi_first
    DO WHILE (phi_l .GT. phi_last)
       phi_l = phi_l - phi_span
    END DO
    DO WHILE (phi_l .LT. phi_first)
       phi_l = phi_l + phi_span
    END DO
    ! fix period if necessary
    phi_last  = fp%coords%x2(nstep)
    DO WHILE (phi_l .GT. phi_last)
       fp => fp%next
       phi_last  = fp%coords%x2(nstep)
    END DO
    phi_first = fp%coords%x2(0)
    DO WHILE (phi_l .LT. phi_first)
       fp => fp%prev
       phi_first  = fp%coords%x2(0)
    END DO

    CALL plagrange_stencel(fp,phi_l,nlagrange,stencel)
    CALL plagrange_value(fp,stencel,phi_arr,bhat_arr)
    ALLOCATE( coeff(0:nder,nlagrange+1) )
    CALL plagrange_coeff(nlagrange+1,nder,phi_l,phi_arr(:),coeff(:,:))
    bhat = SUM(coeff(0,:)*bhat_arr)
    bhatder = SUM(coeff(1,:)*bhat_arr)

!!$    !if (abs(phi_l) .lt. 0.2) then
!!$       print *, 'lagrange  ',stencel
!!$       print *, 'lagrange  ',phi_arr
!!$       print *, 'lagrange  ',bhat_arr
!!$       print *, 'lagrange  ',phi_l,bhat
!!$       print *, ' '
!!$    !end if

    IF (ALLOCATED(stencel)) DEALLOCATE(stencel)
    IF (ALLOCATED(phi_arr)) DEALLOCATE(phi_arr)
    IF (ALLOCATED(bhat_arr)) DEALLOCATE(bhat_arr)
    IF (ALLOCATED(coeff)) DEALLOCATE(coeff)

  END SUBROUTINE plag_interp_2
  !---------------------------------------------------------------------
  SUBROUTINE plag_interp_3(fieldperiod,phi,nlagrange,bhat,bhatder,bhatdder)
    TYPE(fieldperiod_struct), POINTER                 :: fieldperiod, fp
    REAL(kind=dp), INTENT(in)                         :: phi
    INTEGER, INTENT(in)                               :: nlagrange
    REAL(kind=dp), INTENT(out)                        :: bhat,bhatder,bhatdder

    INTEGER, PARAMETER :: nder = 2
    INTEGER, DIMENSION(:), ALLOCATABLE :: stencel
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: phi_arr,bhat_arr
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: coeff


    INTEGER :: nstep
    REAL(kind=dp) ::  phi_l
    REAL(kind=dp) ::  phi_span,phi_last,phi_first

    fp => fieldperiod

    nstep = UBOUND(fp%coords%x2,1)
    phi_l = phi
    ! fix phi_l locally to be within range
    phi_first = fp%parent%ch_fir%coords%x2(0)
    phi_last  = fp%parent%ch_las%coords%x2(nstep)
    phi_span  = phi_last - phi_first
    DO WHILE (phi_l .GT. phi_last)
       phi_l = phi_l - phi_span
    END DO
    DO WHILE (phi_l .LT. phi_first)
       phi_l = phi_l + phi_span
    END DO
    ! fix period if necessary
    phi_last  = fp%coords%x2(nstep)
    DO WHILE (phi_l .GT. phi_last)
       fp => fp%next
       phi_last  = fp%coords%x2(nstep)
    END DO
    phi_first = fp%coords%x2(0)
    DO WHILE (phi_l .LT. phi_first)
       fp => fp%prev
       phi_first  = fp%coords%x2(0)
    END DO

    CALL plagrange_stencel(fp,phi_l,nlagrange,stencel)
    CALL plagrange_value(fp,stencel,phi_arr,bhat_arr)
    ALLOCATE( coeff(0:nder,nlagrange+1) )
    CALL plagrange_coeff(nlagrange+1,nder,phi_l,phi_arr(:),coeff(:,:))
    bhat = SUM(coeff(0,:)*bhat_arr)
    bhatder = SUM(coeff(1,:)*bhat_arr)
    bhatdder = SUM(coeff(2,:)*bhat_arr)

!!$    !if (abs(phi_l) .lt. 0.2) then
!!$       print *, 'lagrange  ',stencel
!!$       print *, 'lagrange  ',phi_arr
!!$       print *, 'lagrange  ',bhat_arr
!!$       print *, 'lagrange  ',phi_l,bhat
!!$       print *, ' '
!!$    !end if

    IF (ALLOCATED(stencel)) DEALLOCATE(stencel)
    IF (ALLOCATED(phi_arr)) DEALLOCATE(phi_arr)
    IF (ALLOCATED(bhat_arr)) DEALLOCATE(bhat_arr)
    IF (ALLOCATED(coeff)) DEALLOCATE(coeff)

  END SUBROUTINE plag_interp_3
  !---------------------------------------------------------------------
  SUBROUTINE plag_interp_all(fieldperiod,phi,nlagrange,x1,x3,bhat,geodcu,h_phi,dlogbdphi)
    TYPE(fieldperiod_struct), POINTER                 :: fieldperiod, fp
    REAL(kind=dp), INTENT(in)                         :: phi
    INTEGER, INTENT(in)                               :: nlagrange
    REAL(kind=dp), INTENT(out)                        :: x1,x3,bhat,geodcu,h_phi,dlogbdphi

    INTEGER, PARAMETER :: nder = 0
    INTEGER, DIMENSION(:), ALLOCATABLE :: stencel
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: phi_arr
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: x1_arr,x3_arr,bhat_arr
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: geodcu_arr,h_phi_arr,dlogbdphi_arr
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: coeff

    INTEGER :: nstep
    REAL(kind=dp) ::  phi_l
    REAL(kind=dp) ::  phi_span,phi_last,phi_first
    fp => fieldperiod
    nstep = UBOUND(fp%coords%x2,1)
    phi_l = phi
    ! fix phi_l locally to be within range
    phi_first = fp%parent%ch_fir%coords%x2(0)
    phi_last  = fp%parent%ch_las%coords%x2(nstep)
    phi_span  = phi_last - phi_first
    DO WHILE (phi_l .GT. phi_last)
       phi_l = phi_l - phi_span
    END DO
    DO WHILE (phi_l .LT. phi_first)
       phi_l = phi_l + phi_span
    END DO
    ! fix period if necessary
    phi_last  = fp%coords%x2(nstep)
    DO WHILE (phi_l .GT. phi_last)
       fp => fp%next
       phi_last  = fp%coords%x2(nstep)
    END DO
    phi_first = fp%coords%x2(0)
    DO WHILE (phi_l .LT. phi_first)
       fp => fp%prev
       phi_first  = fp%coords%x2(0)
    END DO


    CALL plagrange_stencel(fp,phi_l,nlagrange,stencel)
    CALL plagrange_value(fp,stencel,phi_arr, &
         x1_arr,x3_arr,bhat_arr,geodcu_arr,h_phi_arr,dlogbdphi_arr)
    ALLOCATE( coeff(0:nder,nlagrange+1) )
    CALL plagrange_coeff(nlagrange+1,nder,phi_l,phi_arr(:),coeff(:,:))
    x1   = SUM(coeff(0,:)*x1_arr)
    x3   = SUM(coeff(0,:)*x3_arr)
    bhat = SUM(coeff(0,:)*bhat_arr)
    geodcu = SUM(coeff(0,:)*geodcu_arr)
    h_phi = SUM(coeff(0,:)*h_phi_arr)
    dlogbdphi = SUM(coeff(0,:)*dlogbdphi_arr)

    IF (ALLOCATED(stencel)) DEALLOCATE(stencel)
    IF (ALLOCATED(phi_arr)) DEALLOCATE(phi_arr)
    IF (ALLOCATED(x1_arr)) DEALLOCATE(x1_arr)
    IF (ALLOCATED(x3_arr)) DEALLOCATE(x3_arr)
    IF (ALLOCATED(bhat_arr)) DEALLOCATE(bhat_arr)
    IF (ALLOCATED(geodcu_arr)) DEALLOCATE(geodcu_arr)
    IF (ALLOCATED(h_phi_arr)) DEALLOCATE(h_phi_arr)
    IF (ALLOCATED(dlogbdphi_arr)) DEALLOCATE(dlogbdphi_arr)
    IF (ALLOCATED(coeff)) DEALLOCATE(coeff)
  END SUBROUTINE plag_interp_all
  !---------------------------------------------------------------------
  !! Modifications by Andreas F. Martitsch (13.03.2014)
  ! Optional output (necessary for modeling the magnetic rotation)
  SUBROUTINE plag_interp_all2(fieldperiod,phi,nlagrange,x1,x3,bhat,geodcu,&
       h_phi,dlogbdphi,dbcovar_s_hat_dphi,bcovar_s_hat,dlogbds)
    TYPE(fieldperiod_struct), POINTER                 :: fieldperiod, fp
    REAL(kind=dp), INTENT(in)                         :: phi
    INTEGER, INTENT(in)                               :: nlagrange
    REAL(kind=dp), INTENT(out)                        :: x1,x3,bhat,geodcu,h_phi,dlogbdphi
    REAL(kind=dp), INTENT(out)                        :: dbcovar_s_hat_dphi,bcovar_s_hat
    REAL(kind=dp), INTENT(out)                        :: dlogbds

    INTEGER, PARAMETER :: nder = 0
    INTEGER, DIMENSION(:), ALLOCATABLE :: stencel
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: phi_arr
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: x1_arr,x3_arr,bhat_arr
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: geodcu_arr,h_phi_arr,dlogbdphi_arr
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: dbcovar_s_hat_dphi_arr,bcovar_s_hat_arr
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: dlogbds_arr
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: coeff

    INTEGER :: nstep
    REAL(kind=dp) ::  phi_l
    REAL(kind=dp) ::  phi_span,phi_last,phi_first
    fp => fieldperiod
    nstep = UBOUND(fp%coords%x2,1)
    phi_l = phi
    ! fix phi_l locally to be within range
    phi_first = fp%parent%ch_fir%coords%x2(0)
    phi_last  = fp%parent%ch_las%coords%x2(nstep)
    phi_span  = phi_last - phi_first
    DO WHILE (phi_l .GT. phi_last)
       phi_l = phi_l - phi_span
    END DO
    DO WHILE (phi_l .LT. phi_first)
       phi_l = phi_l + phi_span
    END DO
    ! fix period if necessary
    phi_last  = fp%coords%x2(nstep)
    DO WHILE (phi_l .GT. phi_last)
       fp => fp%next
       phi_last  = fp%coords%x2(nstep)
    END DO
    phi_first = fp%coords%x2(0)
    DO WHILE (phi_l .LT. phi_first)
       fp => fp%prev
       phi_first  = fp%coords%x2(0)
    END DO


    CALL plagrange_stencel(fp,phi_l,nlagrange,stencel)
    CALL plagrange_value(fp,stencel,phi_arr, &
         x1_arr,x3_arr,bhat_arr,geodcu_arr,h_phi_arr,dlogbdphi_arr,&
         dbcovar_s_hat_dphi_arr,bcovar_s_hat_arr,dlogbds_arr)
    ALLOCATE( coeff(0:nder,nlagrange+1) )
    CALL plagrange_coeff(nlagrange+1,nder,phi_l,phi_arr(:),coeff(:,:))
    x1                 = SUM(coeff(0,:)*x1_arr)
    x3                 = SUM(coeff(0,:)*x3_arr)
    bhat               = SUM(coeff(0,:)*bhat_arr)
    geodcu             = SUM(coeff(0,:)*geodcu_arr)
    h_phi              = SUM(coeff(0,:)*h_phi_arr)
    dlogbdphi          = SUM(coeff(0,:)*dlogbdphi_arr)
    dbcovar_s_hat_dphi = SUM(coeff(0,:)*dbcovar_s_hat_dphi_arr)
    bcovar_s_hat       = SUM(coeff(0,:)*bcovar_s_hat_arr)
    dlogbds            = SUM(coeff(0,:)*dlogbds_arr)

    IF (ALLOCATED(stencel)) DEALLOCATE(stencel)
    IF (ALLOCATED(phi_arr)) DEALLOCATE(phi_arr)
    IF (ALLOCATED(x1_arr)) DEALLOCATE(x1_arr)
    IF (ALLOCATED(x3_arr)) DEALLOCATE(x3_arr)
    IF (ALLOCATED(bhat_arr)) DEALLOCATE(bhat_arr)
    IF (ALLOCATED(geodcu_arr)) DEALLOCATE(geodcu_arr)
    IF (ALLOCATED(h_phi_arr)) DEALLOCATE(h_phi_arr)
    IF (ALLOCATED(dlogbdphi_arr)) DEALLOCATE(dlogbdphi_arr)
    IF (ALLOCATED(dbcovar_s_hat_dphi_arr)) DEALLOCATE(dbcovar_s_hat_dphi_arr)
    IF (ALLOCATED(bcovar_s_hat_arr)) DEALLOCATE(bcovar_s_hat_arr)
    IF (ALLOCATED(dlogbds_arr)) DEALLOCATE(dlogbds_arr)
    IF (ALLOCATED(coeff)) DEALLOCATE(coeff)
  END SUBROUTINE plag_interp_all2
  !! End Modifications by Andreas F. Martitsch (13.03.2014)
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------

  subroutine fix_phiplacement_problem(ibeg,iend,npart,subsqmin,        &
                                      phi_mfl,bhat_mfl,eta)

    use device_mod, only : fieldpropagator

    implicit none

    integer :: i,ibeg,iend,npart,istep,ibmin,npassing,npassing_prev
    integer :: ncross_l,ncross_r

    double precision :: subsqmin

    integer, dimension(1)              :: idummy
    integer, dimension(:), allocatable :: icross_l,icross_r

    double precision, dimension(0:npart)        :: eta
    double precision, dimension(ibeg:iend)      :: phi_mfl,bhat_mfl
    double precision, dimension(:), allocatable :: eta_cross_l,eta_cross_r

    npassing = -1


    ! determine level crossings:

    idummy=minloc(bhat_mfl(ibeg:iend))
    ibmin=idummy(1)+ibeg-1

    ncross_l=0
    if (ibmin .gt. ibeg) then
      istep=ibmin
      do i=0,npart
        if(1.d0-bhat_mfl(istep)*eta(i) .gt. subsqmin) then
          npassing=i
        else
          exit
        end if
      end do
      npassing_prev=npassing
      do istep=ibmin-1,ibeg,-1
        do i=0,npart
          if (1.d0-bhat_mfl(istep)*eta(i) .gt. subsqmin) then
            npassing=i
          else
            exit
          end if
        end do
        if (npassing.lt.npassing_prev) then
          ncross_l=ncross_l+1
          npassing_prev=npassing
        end if
      end do
      if (ncross_l.gt.0) then
        allocate(icross_l(ncross_l),eta_cross_l(ncross_l))
        ncross_l=0
        istep=ibmin
        do i=0,npart
          if(1.d0-bhat_mfl(istep)*eta(i) .gt. subsqmin) then
            npassing=i
          else
            exit
          end if
        end do
        npassing_prev=npassing
        do istep=ibmin-1,ibeg,-1
          do i=0,npart
            if(1.d0-bhat_mfl(istep)*eta(i) .gt. subsqmin) then
              npassing=i
            else
              exit
            end if
          end do
          if (npassing.lt.npassing_prev) then
            ncross_l=ncross_l+1
            icross_l(ncross_l)=istep
            eta_cross_l(ncross_l)=eta(npassing_prev)
            npassing_prev=npassing
          end if
        end do
        do i=1,ncross_l
          istep=icross_l(i)
          if (abs(bhat_mfl(istep-1)*eta_cross_l(i)-1.d0).lt. &
              abs(bhat_mfl(istep)  *eta_cross_l(i)-1.d0)) then
            open(111,file='phi_placement_problem.dat',position='append')
            write(111,*) ' propagator tag = ',fieldpropagator%tag, &
                         ' step number = ',istep-1,                &
                         ' 1 / bhat = ',1.d0/bhat_mfl(istep-1),    &
                         ' eta = ',eta_cross_l(i)
            close(111)
            bhat_mfl(istep-1)=1/eta_cross_l(i)
          elseif (abs(bhat_mfl(istep+1)*eta_cross_l(i)-1.d0).lt. &
                  abs(bhat_mfl(istep)  *eta_cross_l(i)-1.d0)) then
            open(111,file='phi_placement_problem.dat',position='append')
            write(111,*) ' propagator tag = ',fieldpropagator%tag, &
                         ' step number = ',istep+1,                &
                         ' 1 / bhat = ',1.d0/bhat_mfl(istep+1),    &
                         ' eta = ',eta_cross_l(i)
            bhat_mfl(istep+1)=1/eta_cross_l(i)
            close(111)
          end if
        end do
        deallocate(icross_l,eta_cross_l)
      end if
    end if

    ncross_r=0
    if (ibmin.lt.iend) then
      istep=ibmin
      do i=0,npart
        if (1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) then
          npassing=i
        else
          exit
        end if
      end do
      npassing_prev=npassing
      do istep=ibmin+1,iend
        do i=0,npart
          if (1.d0-bhat_mfl(istep)*eta(i).gt.subsqmin) then
            npassing=i
          else
            exit
          end if
        end do
        if (npassing.lt.npassing_prev) then
          ncross_r=ncross_r+1
          npassing_prev=npassing
        end if
      end do
      if (ncross_r.gt.0) then
        allocate(icross_r(ncross_r),eta_cross_r(ncross_r))
        ncross_r=0
        istep=ibmin
        do i=0,npart
          if (1.d0-bhat_mfl(istep)*eta(i).gt.subsqmin) then
            npassing=i
          else
            exit
          end if
        end do
        npassing_prev=npassing
        do istep=ibmin+1,iend
          do i=0,npart
            if (1.d0-bhat_mfl(istep)*eta(i).gt.subsqmin) then
              npassing=i
            else
              exit
            end if
          end do
          if (npassing.lt.npassing_prev) then
            ncross_r=ncross_r+1
            icross_r(ncross_r)=istep
            eta_cross_r(ncross_r)=eta(npassing_prev)
            npassing_prev=npassing
          end if
        end do
        do i=1,ncross_r
          istep=icross_r(i)
          if (abs(bhat_mfl(istep-1)*eta_cross_r(i)-1.d0).lt. &
              abs(bhat_mfl(istep)  *eta_cross_r(i)-1.d0)) then
            open(111,file='phi_placement_problem.dat',position='append')
            write(111,*) ' propagator tag = ',fieldpropagator%tag, &
                         ' step number = ',istep-1,                &
                         ' 1 / bhat = ',1.d0/bhat_mfl(istep-1),    &
                         ' eta = ',eta_cross_r(i)
            close(111)
            bhat_mfl(istep-1)=1/eta_cross_r(i)
          elseif (abs(bhat_mfl(istep+1)*eta_cross_r(i)-1.d0).lt. &
                  abs(bhat_mfl(istep)  *eta_cross_r(i)-1.d0)) then
            open(111,file='phi_placement_problem.dat',position='append')
            write(111,*) ' propagator tag = ',fieldpropagator%tag, &
                         ' step number = ',istep+1,                &
                         ' 1 / bhat = ',1.d0/bhat_mfl(istep+1),    &
                         ' eta = ',eta_cross_r(i)
            close(111)
            bhat_mfl(istep+1)=1/eta_cross_r(i)
          end if
        end do
        deallocate(icross_r,eta_cross_r)
      end if
    end if

  end subroutine fix_phiplacement_problem

END MODULE plagrange_mod
