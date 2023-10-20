!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  SUBROUTINE pollag(n,coeflag)
! Computes coefficients of Laguerre polynomials of orders from 0 to n
!
!
! Input parameters:
!           Formal: n            - maximum order of Leguerre polynomials
! Output parameters:
!           Formal: coeflag(i,j) - j-th coefficient of re-normalized Leguerre 
!                                  polynomial (function $\varphi$)
!                                  of the order i
  IMPLICIT NONE

  INTEGER :: n,i,j,k

  DOUBLE PRECISION, PARAMETER  :: pi=3.14159265358979d0
  DOUBLE PRECISION :: bincoef,hm
  DOUBLE PRECISION, DIMENSION(0:n,0:n) :: coeflag

  coeflag=0.d0

  hm=3.d0/(8.d0*pi)
  coeflag(0,0)=1.d0/sqrt(hm)

  DO i=1,n
    hm=hm*(1.d0+1.5d0/i)
    DO k=0,i
      bincoef=1.d0
      do j=1,i-k
        bincoef=bincoef*(1.d0+(k+1.5d0)/j)
      enddo
      do j=1,k
        bincoef=-bincoef/j
      enddo
      coeflag(i,k)=bincoef
    ENDDO
    coeflag(i,0:i)=coeflag(i,0:i)/sqrt(hm)
  ENDDO

  END SUBROUTINE pollag


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine lagxmm(lag,x1mm,x2mm)

  implicit none

  double precision, parameter :: pi=3.14159265358979d0

  integer :: lag,m,mm,k,kk
  double precision :: cnorm
  double precision, dimension(0:lag,0:lag)    :: x1mm,x2mm,coeflag
  double precision, dimension(:), allocatable :: factorial

  cnorm=1.d0/(2.d0*pi*sqrt(pi))

  allocate(factorial(2*lag+2))

  factorial(1)=1.d0
  do k=2,2*lag+2
    factorial(k)=factorial(k-1)*real(k, kind=kind(0d0))
  enddo

  call pollag(lag,coeflag)

  do m=0,lag
    do mm=0,lag
      x1mm(m,mm)=0.d0
      x2mm(m,mm)=0.d0
      do k=0,m
        do kk=0,mm
          x1mm(m,mm)=x1mm(m,mm)+coeflag(m,k)*coeflag(mm,kk)*factorial(k+kk+1)
          x2mm(m,mm)=x2mm(m,mm)+coeflag(m,k)*coeflag(mm,kk)*factorial(k+kk+2)
        enddo
      enddo
    enddo
  enddo

  x1mm=x1mm*cnorm
  x2mm=x2mm*cnorm

  end subroutine lagxmm
