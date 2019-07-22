MODULE partpa_mod
  ! Winny ipmin,iminr added
  INTEGER                                     :: ipmax,ipmin,istart,iminr
  INTEGER                                     :: npassing_start,npass_max
  DOUBLE PRECISION                            :: pardeb0,bmod0,dirint,hxeta
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eta
  ! Winny for talking between
  INTEGER :: nhalfstep
  ! Winny end
END MODULE
