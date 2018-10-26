module inthecore_mod
  use nrtype, only : dp

  logical :: prop=.true.
  integer :: npoi,ijumpb,ibeg,iend
  real(kind=dp), parameter :: epssep=1.d-6
  real(kind=dp) :: rc,zc,sig,psi_sep,psi_cut,sigpsi,cutoff
  real(kind=dp), dimension(:), allocatable :: rho2i,theti
  integer :: incore
  real(kind=dp) :: vacf,dvacdr,dvacdz,d2vacdr2,d2vacdrdz,d2vacdz2
  real(kind=dp) :: plaf,dpladr,dpladz,d2pladr2,d2pladrdz,d2pladz2
end module inthecore_mod
