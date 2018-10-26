module amn_mod
  use nrtype, only : dp, dpc

  ! Fourier ampitudes of the original field:
  integer :: ntor_amn=1,mpol,ntor_ff,mpol_ff,nsqpsi,icall=0
  real(kind=dp) :: sqpsimin,sqpsimax,hsqpsi
  complex(kind=dpc), dimension(:,:,:,:), allocatable :: splapsi,splatet
  complex(kind=dpc), dimension(:,:), allocatable :: amnpsi,   amntet,     &
                                                 amnpsi_s, amntet_s,   &
                                                 amnpsi_ss,amntet_ss
  double complex, dimension(:),   allocatable :: expthe,expphi
  ! Formfactors:
  integer :: nsqpsi_ff,nmodes_ff
  real(kind=dp) :: sqpsimin_ff,sqpsimax_ff,hsqpsi_ff
  integer,        dimension(:,:), allocatable :: ipoi_ff
  complex(kind=dpc), dimension(:,:,:), allocatable :: splffp,splfft
  complex(kind=dpc), dimension(:),   allocatable :: fmnpsi,   fmntet,     &
                                                 fmnpsi_s, fmntet_s,   &
                                                 fmnpsi_ss,fmntet_ss
end module amn_mod
