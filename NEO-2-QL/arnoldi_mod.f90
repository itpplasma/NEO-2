MODULE arnoldi_mod
  INTEGER :: ngrow,ierr
  INTEGER :: ntol,mode=0
  DOUBLE PRECISION :: tol
  DOUBLE COMPLEX, DIMENSION(:),   ALLOCATABLE :: fzero,ritznum
  DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE :: eigvecs
END MODULE arnoldi_mod
