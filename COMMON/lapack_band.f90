MODULE lapack_band

  IMPLICIT NONE

  PUBLIC mat2band,band2mat,gbsv,gbtrs,gesv

  INTERFACE mat2band
     MODULE PROCEDURE mat2band_d, mat2band_i
  END INTERFACE
  
  INTERFACE band2mat
     MODULE PROCEDURE band2mat_d, band2mat_i
  END INTERFACE
  
  INTERFACE gbsv
     MODULE PROCEDURE gbsv_b1, gbsv_b2, gbsv_b2_ind, gbsv_b2_ind_b
  END INTERFACE
  
  INTERFACE gbtrs
     MODULE PROCEDURE gbtrs_b1, gbtrs_b2, gbtrs_b2_ind, gbtrs_b2_ind_b
  END INTERFACE
  
  INTERFACE gesv
     MODULE PROCEDURE gesv_b2,gesv_b2_ind
  END INTERFACE
  
CONTAINS
  
  ! converts (n x n)-Matrix into band storage
  ! band must be allocatable and is allocated inside the routine
  SUBROUTINE mat2band_d(mat,kl,ku,band)
    use nrtype, only : dp

    REAL(kind=dp), DIMENSION(:,:), INTENT(in)                  :: mat
    INTEGER                                                    :: kl,ku
    REAL(kind=dp), DIMENSION(:,:), INTENT(out), ALLOCATABLE    :: band
    
    INTEGER                                                    :: n,ldab
    INTEGER                                                    :: m,i1,i2
    INTEGER                                                    :: icol
    INTEGER                                                    :: irow,irowb

    n = SIZE(mat,1)
    ldab = 2*kl + ku + 1
    IF (ALLOCATED(band)) DEALLOCATE(band)
    ALLOCATE(band(ldab,n))

    m = kl + ku + 1
    DO icol = 1, n
       i1 = MAX(1, icol - ku)
       i2 = MIN(n, icol + kl)
       DO irow = i1, i2
          irowb = irow - icol + m
          band(irowb,icol) = mat(irow,icol)
       END DO
    END DO
    
  END SUBROUTINE mat2band_d

  SUBROUTINE mat2band_i(mat,kl,ku,band)
    
    INTEGER,       DIMENSION(:,:), INTENT(in)                  :: mat
    INTEGER                                                    :: kl,ku
    INTEGER,       DIMENSION(:,:), INTENT(out), ALLOCATABLE    :: band
    
    INTEGER                                                    :: n,ldab
    INTEGER                                                    :: m,i1,i2
    INTEGER                                                    :: icol
    INTEGER                                                    :: irow,irowb

    n = SIZE(mat,1)
    ldab = 2*kl + ku + 1
    IF (ALLOCATED(band)) DEALLOCATE(band)
    ALLOCATE(band(ldab,n))

    m = kl + ku + 1
    DO icol = 1, n
       i1 = MAX(1, icol - ku)
       i2 = MIN(n, icol + kl)
       DO irow = i1, i2
          irowb = irow - icol + m
          band(irowb,icol) = mat(irow,icol)
       END DO
    END DO
    
  END SUBROUTINE mat2band_i

  ! converts band storage into (n x n)-matrix
  ! mat must be allocatable and is allocated inside the routine
  SUBROUTINE band2mat_d(band,kl,ku,mat)
    use nrtype, only : dp

    REAL(kind=dp), DIMENSION(:,:), INTENT(in)                  :: band
    INTEGER                                                    :: kl,ku
    REAL(kind=dp), DIMENSION(:,:), INTENT(out), ALLOCATABLE    :: mat
    
    INTEGER                                                    :: n,ldab
    INTEGER                                                    :: m,i1,i2
    INTEGER                                                    :: icol
    INTEGER                                                    :: irow,irowb

    n = SIZE(band,2)
    ldab = 2*kl + ku + 1
    IF (ALLOCATED(mat)) DEALLOCATE(mat)
    ALLOCATE(mat(n,n))

    m = kl + ku + 1
    DO icol = 1, n
       i1 = MAX(1, icol - ku)
       i2 = MIN(n, icol + kl)
       DO irow = i1, i2
          irowb = irow - icol + m
          mat(irow,icol) = band(irowb,icol)
       END DO
    END DO
    
  END SUBROUTINE band2mat_d

!  Purpose
  SUBROUTINE band2mat_i(band,kl,ku,mat)
    
    INTEGER,       DIMENSION(:,:), INTENT(in)                  :: band
    INTEGER                                                    :: kl,ku
    INTEGER,       DIMENSION(:,:), INTENT(out), ALLOCATABLE    :: mat
    
    INTEGER                                                    :: n,ldab
    INTEGER                                                    :: m,i1,i2
    INTEGER                                                    :: icol
    INTEGER                                                    :: irow,irowb

    n = SIZE(band,2)
    ldab = 2*kl + ku + 1
    IF (ALLOCATED(mat)) DEALLOCATE(mat)
    ALLOCATE(mat(n,n))

    m = kl + ku + 1
    DO icol = 1, n
       i1 = MAX(1, icol - ku)
       i2 = MIN(n, icol + kl)
       DO irow = i1, i2
          irowb = irow - icol + m
          mat(irow,icol) = band(irowb,icol)
       END DO
    END DO
    
  END SUBROUTINE band2mat_i

!  =======
!
!  DGBSV computes the solution to a real system of linear equations
!  A * X = B, where A is a band matrix of order N with KL subdiagonals
!  and KU superdiagonals, and X and B are N-by-NRHS matrices.
!
!  The LU decomposition with partial pivoting and row interchanges is
!  used to factor A as A = L * U, where L is a product of permutation
!  and unit lower triangular matrices with KL subdiagonals, and U is
!  upper triangular with KL+KU superdiagonals.  The factored form of A
!  is then used to solve the system of equations A * X = B.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of linear equations, i.e., the order of the
!          matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (N)
!          The pivot indices that define the permutation matrix P;
!          row i of the matrix was interchanged with row IPIV(i).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the N-by-NRHS right hand side matrix B.
!          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
!                has been completed, but the factor U is exactly
!                singular, and the solution has not been computed.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U because of fill-in resulting from the row interchanges.
!
!  =====================================================================


  ! Fortran 90 interface with b being a vector (1D)
  ! b is changed on output
  SUBROUTINE gbsv_b1(kl,ku,a,ipivot,b,info)
    use nrtype, only : dp

    INTEGER, INTENT(in)                          :: kl
    INTEGER, INTENT(in)                          :: ku
    REAL(kind=dp), DIMENSION(:,:), INTENT(in)    :: a
    INTEGER,       DIMENSION(:),   INTENT(out)   :: ipivot
    REAL(kind=dp), DIMENSION(:),   INTENT(inout) :: b
    INTEGER,                       INTENT(out)   :: info
    
    INTEGER,       PARAMETER                     :: nrhs = 1
    INTEGER                                      :: n,ldab,ldb
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: b2

    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: a_band
    CALL mat2band(a,kl,ku,a_band)


    n = SIZE(a_band,2)
    ALLOCATE(b2(n,1))
    b2(:,1) = b
    ldb = n
    ldab = 2*kl + ku + 1

    CALL dgbsv(n,kl,ku,nrhs,a_band,ldab,ipivot,b2,ldb,info)
    b = b2(:,1)
    DEALLOCATE(b2)
    DEALLOCATE(a_band)

  END SUBROUTINE gbsv_b1

  ! Fortran 90 interface with b being a matrix (2D)
  ! b is changed on output
  SUBROUTINE gbsv_b2(kl,ku,a,ipivot,b,info)
    use nrtype, only : dp

    INTEGER, INTENT(in)                          :: kl
    INTEGER, INTENT(in)                          :: ku
    REAL(kind=dp), DIMENSION(:,:), INTENT(in)    :: a
    INTEGER,       DIMENSION(:),   INTENT(out)   :: ipivot
    REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER,                       INTENT(out)   :: info
    
    INTEGER                                      :: n,nrhs,ldab,ldb

    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: a_band
    CALL mat2band(a,kl,ku,a_band)

    n = SIZE(a_band,2)
    ldb = n
    nrhs = SIZE(b,2)
    ldab = 2*kl + ku + 1

    CALL dgbsv(n,kl,ku,nrhs,a_band,ldab,ipivot,b,ldb,info)
    DEALLOCATE(a_band)

  END SUBROUTINE gbsv_b2

!!$  SUBROUTINE gbsv_b2_ind(kl,ku,a,ipivot,b,info,ind_s,ind_e)
!!$
!!$    INTEGER, INTENT(in)                          :: kl
!!$    INTEGER, INTENT(in)                          :: ku
!!$    REAL(kind=dp), DIMENSION(:,:), INTENT(in)    :: a
!!$    INTEGER,       DIMENSION(:),   INTENT(out)   :: ipivot
!!$    REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
!!$    INTEGER,                       INTENT(out)   :: info
!!$    INTEGER,       DIMENSION(:),   INTENT(in)    :: ind_s,ind_e
!!$    
!!$    INTEGER                                      :: n,nrhs,ldab,ldb
!!$
!!$    INTEGER                                      :: na,ns
!!$    INTEGER                                      :: nrhsr = 1
!!$    INTEGER                                      :: i_s,i_e
!!$    INTEGER                                      :: k
!!$
!!$    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: a_band
!!$    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: ar,br
!!$    INTEGER,       DIMENSION(:),   ALLOCATABLE   :: ipivotr
!!$
!!$    CALL mat2band(a,kl,ku,a_band)
!!$
!!$    n = SIZE(a,2)
!!$    ldb = n
!!$    nrhs = SIZE(b,2)
!!$    ldab = 2*kl + ku + 1
!!$
!!$    na = MAXVAL(ind_e-ind_s)+1
!!$    ALLOCATE( br(na,1) )
!!$
!!$    DO k = 1,nrhs
!!$       i_s = ind_s(k)
!!$       i_e = ind_e(k)
!!$       ns = i_e - i_s + 1
!!$
!!$       ALLOCATE( ar(ldab,ns) )
!!$       ar(:,1:ns) = a_band(:,i_s:i_e)
!!$       ALLOCATE( ipivotr(ns) )
!!$       br(1:ns,1) = b(i_s:i_e,k)
!!$
!!$       CALL dgbsv(ns,kl,ku,nrhsr,ar,ldab,ipivotr,br,na,info)
!!$
!!$       !b(1:i_s-1,k) = 0.0_dp
!!$       b(i_s:i_e,k) = br(1:ns,1)
!!$       !b(i_e+1:n,k) = 0.0_dp       
!!$       DEALLOCATE(ar)
!!$       DEALLOCATE(ipivotr)
!!$    END DO
!!$    info = info
!!$
!!$    DEALLOCATE(a_band)
!!$    DEALLOCATE(br)
!!$
!!$  END SUBROUTINE gbsv_b2_ind


  SUBROUTINE gbsv_b2_ind(kl,ku,a,ipivot,b,info,ind_s,ind_e)
    use nrtype, only : dp

    INTEGER, INTENT(in)                          :: kl
    INTEGER, INTENT(in)                          :: ku
    REAL(kind=dp), DIMENSION(:,:), INTENT(in)    :: a
    INTEGER,       DIMENSION(:),   INTENT(out)   :: ipivot
    REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER,                       INTENT(out)   :: info
    INTEGER,       DIMENSION(:),   INTENT(in)    :: ind_s,ind_e
    
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: a_band
    CALL mat2band(a,kl,ku,a_band)
    CALL gbsv(kl,ku,a_band,ipivot,b,info,ind_s,ind_e,'b')
    DEALLOCATE(a_band)

  END SUBROUTINE gbsv_b2_ind

  SUBROUTINE gbsv_b2_ind_b(kl,ku,a,ipivot,b,info,ind_s,ind_e,meth)
    use nrtype, only : dp

    INTEGER, INTENT(in)                          :: kl
    INTEGER, INTENT(in)                          :: ku
    REAL(kind=dp), DIMENSION(:,:), INTENT(in)    :: a
    INTEGER,       DIMENSION(:),   INTENT(out)   :: ipivot
    REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER,                       INTENT(out)   :: info
    INTEGER,       DIMENSION(:),   INTENT(in)    :: ind_s,ind_e
    CHARACTER(len=1),              INTENT(in)    :: meth 
    
    INTEGER                                      :: n,nrhs,ldab,ldb

    INTEGER                                      :: na,ns
    INTEGER                                      :: nrhsr = 1
    INTEGER                                      :: i_s,i_e
    INTEGER                                      :: k

    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: a_band
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: ar,br
    INTEGER,       DIMENSION(:),   ALLOCATABLE   :: ipivotr

    ipivot = 0

    IF (meth .NE. 'b' .AND. meth .NE. 'B') THEN
       PRINT *, 'Message from gbtrs'
       PRINT *, 'Method must be b or B'
       STOP
    END IF

    n = SIZE(a,2)
    ldb = n
    nrhs = SIZE(b,2)
    ldab = 2*kl + ku + 1

    na = MAXVAL(ind_e-ind_s)+1
    ALLOCATE( br(na,1) )

    DO k = 1,nrhs
       i_s = ind_s(k)
       i_e = ind_e(k)
       ns = i_e - i_s + 1

       ALLOCATE( ar(ldab,ns) )
       ar(:,1:ns) = a(:,i_s:i_e)
       ALLOCATE( ipivotr(ns) )
       br(1:ns,1) = b(i_s:i_e,k)

       CALL dgbsv(ns,kl,ku,nrhsr,ar,ldab,ipivotr,br,na,info)

       !b(1:i_s-1,k) = 0.0_dp
       b(i_s:i_e,k) = br(1:ns,1)
       !b(i_e+1:n,k) = 0.0_dp       
       DEALLOCATE(ar)
       DEALLOCATE(ipivotr)
    END DO
    info = info

    DEALLOCATE(br)

  END SUBROUTINE gbsv_b2_ind_b


!  Purpose
!  =======
!
!  DGBTRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general band matrix A using the LU factorization computed
!  by DGBTRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations.
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
!          Details of the LU factorization of the band matrix A, as
!          computed by DGBTRF.  U is stored as an upper triangular band
!          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
!          the multipliers used during the factorization are stored in
!          rows KL+KU+2 to 2*KL+KU+1.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices; for 1 <= i <= N, row i of the matrix was
!          interchanged with row IPIV(i).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!  Purpose
!  =======
!
!  DGBTRF computes an LU factorization of a real m-by-n band matrix A
!  using partial pivoting with row interchanges.
!
!  This is the blocked version of the algorithm, calling Level 3 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U because of fill-in resulting from the row interchanges.
!
!  =====================================================================
! 


  ! Fortran 90 interface with b being a vector (1D)
  ! lu-factorization is called within the routine
  ! a and b are changed on output
  SUBROUTINE gbtrs_b1(kl,ku,a,ipivot,b,info)
    use nrtype, only : dp

    INTEGER, INTENT(in)                          :: kl
    INTEGER, INTENT(in)                          :: ku
    REAL(kind=dp), DIMENSION(:,:), INTENT(in)    :: a
    INTEGER,       DIMENSION(:),   INTENT(out)   :: ipivot
    REAL(kind=dp), DIMENSION(:),   INTENT(inout) :: b
    INTEGER,                       INTENT(out)   :: info
    
    INTEGER,       PARAMETER                     :: nrhs = 1
    INTEGER                                      :: n,ldab,ldb
    INTEGER                                      :: info2
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: b2
    CHARACTER(len=1)                             :: trans = 'N'

    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: a_band
    CALL mat2band(a,kl,ku,a_band)

    n = SIZE(a_band,2)
    ALLOCATE(b2(n,1))
    b2(:,1) = b
    ldb = n
    ldab = 2*kl + ku + 1

    CALL dgbtrf(n,n,kl,ku,a_band,ldab,ipivot,info)
    CALL dgbtrs(trans,n,kl,ku,nrhs,a_band,ldab,ipivot,b2,ldb,info2)
    info = 100*info + info2
    b = b2(:,1)
    DEALLOCATE(b2)
    DEALLOCATE(a_band)

  END SUBROUTINE gbtrs_b1

  ! Fortran 90 interface with b being a matrix (2D)
  ! lu-factorization is called within the routine
  ! a and b are changed on output
  SUBROUTINE gbtrs_b2(kl,ku,a,ipivot,b,info)
    use nrtype, only : dp

    INTEGER, INTENT(in)                          :: kl
    INTEGER, INTENT(in)                          :: ku
    REAL(kind=dp), DIMENSION(:,:), INTENT(in)    :: a
    INTEGER,       DIMENSION(:),   INTENT(out)   :: ipivot
    REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER,                       INTENT(out)   :: info
    
    INTEGER                                      :: n,nrhs,ldab,ldb
    INTEGER                                      :: info2
    CHARACTER(len=1)                             :: trans = 'N'

    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: a_band
    CALL mat2band(a,kl,ku,a_band)

    n = SIZE(a_band,2)
    ldb = n
    nrhs = SIZE(b,2)
    ldab = 2*kl + ku + 1

    CALL dgbtrf(n,n,kl,ku,a_band,ldab,ipivot,info)
    CALL dgbtrs(trans,n,kl,ku,nrhs,a_band,ldab,ipivot,b,ldb,info2)
    info = 100*info + info2

    DEALLOCATE(a_band)
  END SUBROUTINE gbtrs_b2

  SUBROUTINE gbtrs_b2_ind(kl,ku,a,ipivot,b,info,ind_s,ind_e)
    use nrtype, only : dp

    INTEGER, INTENT(in)                          :: kl
    INTEGER, INTENT(in)                          :: ku
    REAL(kind=dp), DIMENSION(:,:), INTENT(in)    :: a
    INTEGER,       DIMENSION(:),   INTENT(out)   :: ipivot
    REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER,                       INTENT(out)   :: info
    INTEGER,       DIMENSION(:),   INTENT(in)    :: ind_s,ind_e
    
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: a_band
    CALL mat2band(a,kl,ku,a_band)
    CALL gbtrs(kl,ku,a_band,ipivot,b,info,ind_s,ind_e,'b')
    DEALLOCATE(a_band)
  END SUBROUTINE gbtrs_b2_ind

  SUBROUTINE gbtrs_b2_ind_b(kl,ku,a,ipivot,b,info,ind_s,ind_e,meth)
    use nrtype, only : dp

    INTEGER, INTENT(in)                          :: kl
    INTEGER, INTENT(in)                          :: ku
    REAL(kind=dp), DIMENSION(:,:), INTENT(in)    :: a
    INTEGER,       DIMENSION(:),   INTENT(out)   :: ipivot
    REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER,                       INTENT(out)   :: info
    INTEGER,       DIMENSION(:),   INTENT(in)    :: ind_s,ind_e
    CHARACTER(len=1),              INTENT(in)    :: meth 
   
    INTEGER                                      :: n,nrhs,ldab,ldb
    INTEGER                                      :: info2
    CHARACTER(len=1)                             :: trans = 'N'

    INTEGER                                      :: na,ns
    INTEGER                                      :: nrhsr = 1
    INTEGER                                      :: i_s,i_e
    INTEGER                                      :: k

    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: a_band
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: ar,br
    INTEGER,       DIMENSION(:),   ALLOCATABLE   :: ipivotr

    ipivot = 0

    IF (meth .NE. 'b' .AND. meth .NE. 'B') THEN
       PRINT *, 'Message from gbtrs'
       PRINT *, 'Method must be b or B'
       STOP
    END IF

    n = SIZE(a,2)
    ldb = n
    nrhs = SIZE(b,2)
    ldab = 2*kl + ku + 1

    na = MAXVAL(ind_e-ind_s)+1
    ALLOCATE( br(na,1) )

    DO k = 1,nrhs
       i_s = ind_s(k)
       i_e = ind_e(k)
       ns = i_e - i_s + 1

       ALLOCATE( ar(ldab,ns) )
       ar(:,1:ns) = a(:,i_s:i_e)
       ALLOCATE( ipivotr(ns) )
       br(1:ns,1) = b(i_s:i_e,k)

       CALL dgbtrf(ns,ns,kl,ku,ar,ldab,ipivotr,info)
       CALL dgbtrs(trans,ns,kl,ku,nrhsr,ar,ldab,ipivotr,br,na,info2)

       !b(1:i_s-1,k) = 0.0_dp
       b(i_s:i_e,k) = br(1:ns,1)
       !b(i_e+1:n,k) = 0.0_dp       
       DEALLOCATE(ar)
       DEALLOCATE(ipivotr)
    END DO
    info = 100*info + info2

    DEALLOCATE(br)

  END SUBROUTINE gbtrs_b2_ind_b

  ! conventional method (no band)
  SUBROUTINE gesv_b2(kl,ku,a,ipivot,b,info)
    use nrtype, only : dp

    INTEGER, INTENT(in)                          :: kl
    INTEGER, INTENT(in)                          :: ku
    REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: a
    INTEGER,       DIMENSION(:),   INTENT(out)   :: ipivot
    REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER,                       INTENT(out)   :: info
    
    INTEGER                                      :: n,nrhs


    n = SIZE(a,2)
    nrhs = SIZE(b,2)

    CALL dgesv(n,nrhs,a,n,ipivot,b,n,info)

  END SUBROUTINE gesv_b2

  SUBROUTINE gesv_b2_ind(kl,ku,a,ipivot,b,info,ind_s,ind_e)
    use nrtype, only : dp

    INTEGER, INTENT(in)                          :: kl
    INTEGER, INTENT(in)                          :: ku
    REAL(kind=dp), DIMENSION(:,:), INTENT(in)    :: a
    INTEGER,       DIMENSION(:),   INTENT(out)   :: ipivot
    REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER,                       INTENT(out)   :: info
    INTEGER,       DIMENSION(:),   INTENT(in)    :: ind_s,ind_e
    
    INTEGER                                      :: n,ns,nrhs,na
    INTEGER                                      :: nrhsr = 1
    INTEGER                                      :: i_s,i_e
    INTEGER                                      :: k

    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: ar,br
    
    n = SIZE(a,1)
    na = MAXVAL(ind_e-ind_s)+1
    nrhs = SIZE(b,2)
    ALLOCATE( ar(na,na) )
    ALLOCATE( br(na,1) )

    DO k = 1,nrhs
       i_s = ind_s(k)
       i_e = ind_e(k)
       ns = i_e - i_s + 1
       ar(1:ns,1:ns) = a(i_s:i_e,i_s:i_e)
       br(1:ns,1) = b(i_s:i_e,k)
       CALL dgesv(ns,nrhsr,ar,na,ipivot,br,na,info)
       !b(1:i_s-1,k) = 0.0_dp
       b(i_s:i_e,k) = br(1:ns,1)
       !b(i_e+1:n,k) = 0.0_dp
    END DO
    
    DEALLOCATE(ar)
    DEALLOCATE(br)

  END SUBROUTINE gesv_b2_ind


END MODULE lapack_band
