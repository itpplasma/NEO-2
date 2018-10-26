!> Exchange between flint_bo and rhs_bo
module partpa_bo
  use nrtype, only : dp
  use sizey_bo

  integer                                  :: ipmax
  integer,       dimension(:), allocatable :: isw,ipa,icount
  real(kind=dp)                            :: pard0,bmod0
  real(kind=dp), dimension(:), allocatable :: eta
end module partpa_bo
