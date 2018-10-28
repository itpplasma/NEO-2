module arnoldi_mod
  use nrtype, only : dp, dpc

  integer :: ngrow,ierr
  integer :: ntol,mode=0
  real(kind=dp) :: tol
  complex(kind=dpc), dimension(:),   allocatable :: fzero,ritznum
  complex(kind=dpc), dimension(:,:), allocatable :: eigvecs
end module arnoldi_mod
