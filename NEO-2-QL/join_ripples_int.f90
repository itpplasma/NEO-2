! This routines handle the work between propagtaor and join_ripples
!
! at the moment
!   join_ripples_bsfitsplit
!   join_ripples_bsfitjoin
! cannot be used (not really implemeted)
!
! the join_ripple_nn (new notation) does a few internal things
! which winny needs and then calls the ripple_solver
! when things are settled with the joiner one can remove this call
! and just use the new notation (see naming convention below)
!
!
! ATTENTION change of relevant sizes of arrays affects the whole file
!
! naming convention
!
! one has to use e.g. (same for all relevant variables)
!  o%p%amat_p_p   old value of amat_plus_plus
!  n%p%amat_p_p   new value of amat_plus_plus
!  o%p%source_p_g old value of source_p_g
!  n%p%source_p_g new value of source_p_g
!
! one should be able to use the old join_ripples with a
!  translation of variables
!
!   old name              =>   new name
!   amat_plus_plus        =>   o%p%amat_p_p
!   amat_plus_minus       =>   o%p%amat_p_m
!   amat_plus_plus_new    =>   n%p%amat_p_p
!   amat_plus_minus_new   =>   n%p%amat_p_m
!   source_p              =>   o%p%source_p_g
!                         =>   o%p%source_p_e
!   source_p_new          =>   n%p%source_p_g
!                         =>   n%p%source_p_e
!   ............


SUBROUTINE join_ripples_bsfitsplit(eta,loc)
  ! this is the external routine which has to take care of
  !  splitting information because of a new level in eta
  !
  ! the location for the split is given in the variable loc
  ! eta contains this new value already at position loc

  USE propagator_mod, ONLY: binarysplit, propagator, prop_c, prop_diagnostic, &
         & prop_modifyold, prop_fluxsplitmode
  USE fluxsplit_mod

  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0d0)

  REAL(kind=dp), DIMENSION(0:), INTENT(in)    :: eta
  INTEGER,                      INTENT(in)    :: loc

  LOGICAL                              :: boundary_split

  REAL(kind=dp)                        :: eta_boundary
  REAL(kind=dp)                        :: val
  REAL(kind=dp)                        :: boundary_width,eta_width
  REAL(kind=dp)                        :: wsplit

  INTEGER                              :: mode,ishift,ierr
  INTEGER                              :: j1,j2,ub

  INTEGER                              :: n1,n2,npass,i
  TYPE(propagator), POINTER            :: a

  REAL(kind=dp), ALLOCATABLE           :: work(:,:)
  REAL(kind=dp), ALLOCATABLE           :: eta_splitmat(:)
  REAL(kind=dp), ALLOCATABLE           :: eta_mod(:)
  REAL(kind=dp), ALLOCATABLE           :: splitmatrix(:,:)

  ! pointer
  IF (prop_modifyold .EQ. 1) THEN
     a => prop_c%prev
     eta_boundary = a%p%eta_boundary_r
  ELSE
     a => prop_c
     eta_boundary = a%p%eta_boundary_l
  END IF
  ! sizes
  n1 = SIZE(a%p%cmat,1)
  n2 = SIZE(a%p%cmat,2)
  ub = UBOUND(eta,1)
  ! number of passing without boundary
  npass = n1 - 1
  ! mode for fluxsplitter
  mode = prop_fluxsplitmode

  ! WINNY
  ! here now splitting has to be handled
  !  one always works on a%p%cmat
  !  forward  is on the old
  !  backward is on the new (automatically taken care of by prop_modifyold

  ! compute different weight for splitting if boundary has to be splitted
  IF (loc .EQ. npass+1 .AND. eta(loc) .LT. eta_boundary) THEN
     boundary_split = .TRUE.
     boundary_width = eta_boundary - eta(loc-1)
     eta_width      = 2.0_dp * (eta(loc) - eta(loc-1))
     wsplit         = eta_width / boundary_width / 2.0_dp
  ELSE
     boundary_split = .FALSE.
     wsplit         = 0.5_dp
  END IF

  ! here comes the splitting
  IF (loc .LE. npass .OR. boundary_split) THEN
     ! working array and copy parts which are not effected
     ALLOCATE(work(n1+1,n2))
     work(1:loc-1,:) = a%p%cmat(1:loc-1,:)
     work(loc+2:,:)  = a%p%cmat(loc+1:,:)
     IF (mode .EQ. 0) THEN
        ! the simple mode 0
        work(loc,:)     = wsplit*a%p%cmat(loc,:)
        work(loc+1,:)   = (1.0_dp - wsplit) * a%p%cmat(loc,:)
     ELSE IF (mode .EQ. 1 .OR. mode .EQ. 2) THEN
        ! higher modes
        ! this is a modified eta vector without the already splitted eta
        ! makes things easier
        ALLOCATE(eta_mod(0:ub-1))
        eta_mod(0:loc-1) = eta(0:loc-1)
        eta_mod(loc:)    = eta(loc+1:)
        ! compute ishift and range of eta values to be taken
        ishift = 0
        j1 = loc - mode - 1
        j2 = j1 + 2*mode + 1
        IF (j1 .LT. 0) THEN
           ishift = j1
           j1 = j1 - ishift
           j2 = j2 - ishift
        END IF
        IF (j2 .GT. npass+1) THEN
           ishift = j2 - (npass+1)
           j1 = j1 - ishift
           j2 = j2 - ishift
        END IF
        ! now take eta values which should be given to the fluxsplitter
        ! and call the fluxsplitter
        ALLOCATE(eta_splitmat(2*mode+2))
        IF (j2 .EQ. npass+1) THEN
           eta_splitmat = (/eta_mod(j1:j2-1),eta_boundary/)
        ELSE
           eta_splitmat = eta_mod(j1:j2)
        END IF
        CALL fluxsplitter(mode,ishift,wsplit,eta_splitmat,splitmatrix,ierr)
        IF (ierr .LT. 0) THEN
           PRINT *, 'Error in fluxsplitter: ',ierr
           STOP
        END IF
        ! now where to place it
        j1 = loc - mode
        j2 = loc + mode
        IF (j1 .LT. 1) THEN
           ishift = j1 - 1
           j1 = j1 - ishift
           j2 = j2 - ishift
        END IF
        IF (j2 .GT. n1) THEN
           ishift = j2 - n1
           j1 = j1 - ishift
           j2 = j2 - ishift
        END IF
        ! and place it
        work(loc:loc+1,:) = MATMUL(splitmatrix,a%p%cmat(j1:j2,:))
        ! deallocate
        DEALLOCATE(splitmatrix)
        DEALLOCATE(eta_splitmat)
        DEALLOCATE(eta_mod)
     ELSE
        PRINT *, 'prop_fluxsplitmode not defined: ',prop_fluxsplitmode
        STOP
     END IF

     DEALLOCATE(a%p%cmat)
     ALLOCATE(a%p%cmat(n1+1,n2))
     a%p%cmat = work
     DEALLOCATE(work)
     ! diagnostic
     IF (prop_diagnostic .GE. 3) THEN
        PRINT *, ' '
        PRINT *, 'join_ripples_bsfitsplit ',loc,eta(loc)
     END IF
     IF (prop_diagnostic .GE. 3) THEN
        OPEN(unit=1000,file='cmat.dat')
        DO i = 1,n1+1
           WRITE(1000,'(1000f12.5)') a%p%cmat(i,:)
        END DO
        CLOSE(unit=1000)
        PRINT *, 'join_ripples_bsfitsplit - cmat.dat written'
     END IF
  END IF


  ! final cleaning
  NULLIFY(a)

END SUBROUTINE join_ripples_bsfitsplit

SUBROUTINE join_ripples_bsfitjoin(eta,loc)
  ! this is the external routine which has to take care of
  !  joining information because of a disappearing level in eta
  !
  ! in the array eta this value is already removed
  !
  ! in the test the next value is moved down
  ! in reality two levels have to be merged (sum?)

  USE propagator_mod, only: binarysplit, propagator, prop_c, prop_diagnostic, &
         & prop_modifyold

  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0d0)

  REAL(kind=dp), DIMENSION(0:), INTENT(in)    :: eta
  INTEGER,                      INTENT(in)    :: loc

  REAL(kind=dp)                        :: eta_boundary
  INTEGER                              :: i,n1,n2
  INTEGER                              :: npass

  TYPE(propagator), POINTER            :: a

  REAL(kind=dp), ALLOCATABLE           :: work(:,:)

  ! pointer
  IF (prop_modifyold .EQ. 1) THEN
     a => prop_c%prev
     eta_boundary = a%p%eta_boundary_r
  ELSE
     a => prop_c
     eta_boundary = a%p%eta_boundary_l
  END IF
  ! sizes
  n1 = SIZE(a%p%cmat,1)
  n2 = SIZE(a%p%cmat,2)
  ! number of passing without boundary
  npass = n1 - 1

  ! WINNY
  ! here now joining has to be handled
  !  one always works on a%pcmat
  !  forward  is on the old
  !  backward is on the new (automatically taken care of by prop_modifyold
  IF (loc .LE. npass) THEN
     ALLOCATE(work(n1-1,n2))
     work(1:loc-1,:) = a%p%cmat(1:loc-1,:)
     work(loc,:)     = a%p%cmat(loc,:) + a%p%cmat(loc+1,:)
     work(loc+1:,:)  = a%p%cmat(loc+2:,:)
     DEALLOCATE(a%p%cmat)
     ALLOCATE(a%p%cmat(n1-1,n2))
     a%p%cmat = work
     DEALLOCATE(work)
     ! diagnostic
     IF (prop_diagnostic .GE. 3) THEN
        PRINT *, ' '
        PRINT *, 'join_ripples_bsfitjoin  ',loc,eta(loc)
     END IF
     IF (prop_diagnostic .GE. 3) THEN
        OPEN(unit=1000,file='cmat.dat')
        DO i = 1,n1-1
           WRITE(1000,'(1000f12.5)') a%p%cmat(i,:)
        END DO
        CLOSE(unit=1000)
        PRINT *, 'join_ripples_bsfitjoin - cmat.dat written'
     END IF
  END IF


  ! final cleaning
  NULLIFY(a)

END SUBROUTINE join_ripples_bsfitjoin

SUBROUTINE join_ripples_nn(ierr,cstat)
  ! this is the external subroutine to join ripples
  !  at this point it is ensured that all arrays already have the
  !  same size. So, all differences from binary_split are already
  !  removed
  !
  ! the results of joining should be stored in the old propagator
  ! see below for the naming convention

  USE propagator_mod
  USE binarysplit_mod

  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0d0)

  INTEGER, INTENT(out) :: ierr
  CHARACTER(len=5), INTENT(in) :: cstat

  TYPE(propagator), POINTER            :: o
  TYPE(propagator), POINTER            :: n

  ! initialize
  ierr = 0
  o => prop_c%prev
  n => prop_c

  o%fieldpropagator_tag_e = n%fieldpropagator_tag_e

  ! here comes the real physical process of joining
  ! (see naming convention from top of the file)

  IF (cstat .EQ. 'inter') THEN
     CALL join_ripples(ierr)
  ELSE IF (cstat .EQ. 'final') THEN
     CALL join_ends(ierr)
  ELSE
     PRINT *, 'Error in Joiner - cstat: ',cstat,' not implemeted'
     STOP
  END IF
  ! WINNY
  !
  ! npart??? has one to do something
  !
  ! fix right side of old propagator (it should have now
  ! eta_boundary_r and npass_r from the new one
  o%p%eta_boundary_r   = n%p%eta_boundary_r
  ! o%p%npass_r = n%p%npass_r (does the joiner)

  ! fix also the y value (actually yend of propagator)
  ! and phi_value at the end
  o%y = n%y
  o%phi_r = n%phi_r
  ! fix also the binarysplit information
  o%eta_bs_r = n%eta_bs_r
  ! fix also eta information
  IF (ALLOCATED(o%p%eta_r)) DEALLOCATE(o%p%eta_r)
  ALLOCATE(o%p%eta_r(LBOUND(n%p%eta_r,1):UBOUND(n%p%eta_r,1)))
  o%p%eta_r = n%p%eta_r
  ! final cleaning
  NULLIFY(o)
  NULLIFY(n)

END SUBROUTINE join_ripples_nn
