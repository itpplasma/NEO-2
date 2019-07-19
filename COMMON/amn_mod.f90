module amn_mod
  ! Fourier ampitudes of the original field:
  integer :: ntor_amn=1,mpol,ntor_ff,mpol_ff,nsqpsi,icall=0
  double precision :: sqpsimin,sqpsimax,hsqpsi
  double complex, dimension(:,:,:,:), allocatable :: splapsi,splatet
  double complex, dimension(:,:), allocatable :: amnpsi,   amntet,     &
                                                 amnpsi_s, amntet_s,   &
                                                 amnpsi_ss,amntet_ss
  double complex, dimension(:),   allocatable :: expthe,expphi

  ! Formfactors:
  integer :: nsqpsi_ff,nmodes_ff
  double precision :: sqpsimin_ff,sqpsimax_ff,hsqpsi_ff
  integer,        dimension(:,:), allocatable :: ipoi_ff
  double complex, dimension(:,:,:), allocatable :: splffp,splfft
  double complex, dimension(:),   allocatable :: fmnpsi,   fmntet,     &
                                                 fmnpsi_s, fmntet_s,   &
                                                 fmnpsi_ss,fmntet_ss
end module amn_mod
