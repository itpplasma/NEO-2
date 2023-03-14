
! --------------------------------------------------------------------
!
!  inter_interfaces.f90
!
! --------------------------------------------------------------------

MODULE inter_interfaces


  INTERFACE lubksb
     SUBROUTINE lubksb_a(a,indx,b)
       USE nrtype
       REAL(DP),     DIMENSION(:,:), INTENT(IN)    :: a
       INTEGER(I4B), DIMENSION(:),   INTENT(IN)    :: indx
       REAL(DP),     DIMENSION(:),   INTENT(INOUT) :: b
     END SUBROUTINE lubksb_a
  END INTERFACE


  INTERFACE ludcmp
     SUBROUTINE ludcmp_a(a,indx,d)
       use nrtype, only : I4B, DP
       REAL(DP),     DIMENSION(:,:), INTENT(INOUT) :: a
       INTEGER(I4B), DIMENSION(:),   INTENT(OUT)   :: indx
       REAL(DP),                     INTENT(OUT)   :: d
     END SUBROUTINE ludcmp_a
  END INTERFACE


  INTERFACE splinecof3
     SUBROUTINE splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
          a, b, c, d, m, f)
       use nrtype, only : I4B, DP
       REAL(DP),                   INTENT(INOUT) :: c1, cn
       REAL(DP),     DIMENSION(:), INTENT(IN)    :: x
       REAL(DP),     DIMENSION(:), INTENT(IN)    :: y
       REAL(DP),     DIMENSION(:), INTENT(IN)    :: lambda1
       INTEGER(I4B), DIMENSION(:), INTENT(IN)    :: indx
       REAL(DP),     DIMENSION(:), INTENT(OUT)   :: a, b, c, d
       INTEGER(I4B),               INTENT(IN)    :: sw1, sw2
       REAL(DP),                   INTENT(IN)    :: m
       INTERFACE
          FUNCTION f(x,m)
            use nrtype, only : DP
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP), INTENT(IN) :: m
            REAL(DP)             :: f
          END FUNCTION f
       END INTERFACE
     END SUBROUTINE splinecof3_a
  END INTERFACE

  interface splinecof1
    subroutine splinecof1_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
        & a, b, c, d, m, f)
      use nrtype, only : I4B, DP

      real(DP),                   intent(inout) :: c1, cn
      real(DP),     dimension(:), intent(in)    :: x
      real(DP),     dimension(:), intent(in)    :: y
      real(DP),     dimension(:), intent(in)    :: lambda1
      integer(I4B), dimension(:), intent(in)    :: indx
      real(DP),     dimension(:), intent(out)   :: a, b, c, d
      integer(I4B),               intent(in)    :: sw1, sw2
      real(DP),                   intent(in)    :: m
      interface
        function f(x,m)
          use nrtype, only : DP
          implicit none
          real(DP), intent(in) :: x
          real(DP), intent(in) :: m
          real(DP)             :: f
        end function f
      end interface
    end subroutine splinecof1_a
  end interface


  INTERFACE reconstruction3
     SUBROUTINE reconstruction3_a(ai, bi, ci, di, h, a, b, c, d)
       use nrtype, only : DP
       REAL(DP), INTENT(IN)    :: ai, bi, ci, di
       REAL(DP), INTENT(IN)    :: h
       REAL(DP), INTENT(OUT)   :: a, b, c, d 
     END SUBROUTINE reconstruction3_a
  END INTERFACE

  interface reconstruction1
    subroutine reconstruction1_a(ai, bi, ci, di, h, a, b, c, d)
      use nrtype, only : DP
      real(DP), intent(in)    :: ai, bi, ci, di
      real(DP), intent(in)    :: h
      real(DP), intent(out)   :: a, b, c, d
    end subroutine reconstruction1_a
  end interface

  INTERFACE calc_opt_lambda3
     SUBROUTINE calc_opt_lambda3_a(x, y, lambda)
       use nrtype, only : DP
       REAL(DP), DIMENSION(:), INTENT(IN)  :: x, y
       REAL(DP), DIMENSION(:), INTENT(OUT) :: lambda
     END SUBROUTINE calc_opt_lambda3_a
  END INTERFACE


  INTERFACE dist_lin
     SUBROUTINE dist_lin_a(x, y, ymax, dist)
       use nrtype, only : DP
       REAL(DP), DIMENSION(:), INTENT(IN)  :: x, y
       REAL(DP),               INTENT(IN)  :: ymax
       REAL(DP),               INTENT(OUT) :: dist
     END SUBROUTINE dist_lin_a
  END INTERFACE


  INTERFACE splinecof3_lo_driv
     SUBROUTINE splinecof3_lo_driv_a(x, y, c1, cn, lambda, w, indx, &
          sw1, sw2, a, b, c, d, m, f)
       use nrtype, only : I4B, DP
       INTEGER(I4B), DIMENSION(:), INTENT(IN)    :: indx
       REAL(DP),                   INTENT(IN)    :: m
       REAL(DP),                   INTENT(INOUT) :: c1, cn
       REAL(DP),     DIMENSION(:), INTENT(IN)    :: x
       REAL(DP),     DIMENSION(:), INTENT(IN)    :: y
       REAL(DP),     DIMENSION(:), INTENT(IN)    :: lambda
       INTEGER(I4B), DIMENSION(:), INTENT(IN)    :: w
       REAL(DP),     DIMENSION(:), INTENT(OUT)   :: a, b, c, d
       INTEGER(I4B),               INTENT(IN)    :: sw1, sw2
       INTERFACE
          FUNCTION f(x,m)
            use nrtype, only : DP
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP), INTENT(IN) :: m
            REAL(DP)             :: f
          END FUNCTION f
       END INTERFACE
     END SUBROUTINE splinecof3_lo_driv_a
  END INTERFACE

  interface splinecof1_lo_driv
    subroutine splinecof1_lo_driv_a(x, y, c1, cn, lambda, w, indx, &
        & sw1, sw2, a, b, c, d, m, f)
      use nrtype, only : I4B, DP
      integer(I4B), dimension(:), intent(in)    :: indx
      real(DP),                   intent(in)    :: m
      real(DP),                   intent(inout) :: c1, cn
      real(DP),     dimension(:), intent(in)    :: x
      real(DP),     dimension(:), intent(in)    :: y
      real(DP),     dimension(:), intent(in)    :: lambda
      integer(I4B), dimension(:), intent(in)    :: w
      real(DP),     dimension(:), intent(out)   :: a, b, c, d
      integer(I4B),               intent(in)    :: sw1, sw2
      interface
        function f(x,m)
          use nrtype, only : DP
          implicit none
          real(DP), intent(in) :: x
          real(DP), intent(in) :: m
          real(DP)             :: f
        end function f
      end interface
    end subroutine splinecof1_lo_driv_a
  end interface


  INTERFACE splinecof3_hi_driv
     SUBROUTINE splinecof3_hi_driv_a(x, y, m, a, b, c, d, indx, f)
       use nrtype, only : I4B, DP
       INTEGER(I4B), DIMENSION(:),   INTENT(IN)  :: indx
       REAL(DP),     DIMENSION(:),   INTENT(IN)  :: m
       REAL(DP),     DIMENSION(:),   INTENT(IN)  :: x
       REAL(DP),     DIMENSION(:,:), INTENT(IN)  :: y
       REAL(DP),     DIMENSION(:,:), INTENT(OUT) :: a, b, c, d
       INTERFACE
          FUNCTION f(x,m)
            use nrtype, only : DP
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP), INTENT(IN) :: m
            REAL(DP)             :: f
          END FUNCTION f
       END INTERFACE
     END SUBROUTINE splinecof3_hi_driv_a
  END INTERFACE

  interface splinecof1_hi_driv
    subroutine splinecof1_hi_driv_a(x, y, m, a, b, c, d, indx, f)
      use nrtype, only : I4B, DP
      integer(I4B), dimension(:),   intent(in)  :: indx
      real(DP),     dimension(:),   intent(in)  :: m
      real(DP),     dimension(:),   intent(in)  :: x
      real(DP),     dimension(:,:), intent(in)  :: y
      real(DP),     dimension(:,:), intent(out) :: a, b, c, d
      interface
        function f(x,m)
          use nrtype, only : DP
          implicit none
          real(DP), intent(in) :: x
          real(DP), intent(in) :: m
          real(DP)             :: f
        end function f
      end interface
    end subroutine splinecof1_hi_driv_a
  end interface


  INTERFACE splint_horner3
     SUBROUTINE splint_horner3_a(xa,a,b,c,d,swd,m,x_in,f,fp,fpp,fppp,&
          y,yp,ypp,yppp)
       use nrtype, only : I4B, DP
       INTEGER(I4B),               INTENT(IN)  :: swd
       REAL(DP),                   INTENT(IN)  :: m
       REAL(DP),     DIMENSION(:), INTENT(IN)  :: xa, a, b, c, d
       REAL(DP),                   INTENT(IN)  :: x_in
       REAL(DP),                   INTENT(OUT) :: y, yp, ypp, yppp
       INTERFACE
          FUNCTION f(x,m)
            use nrtype, only : DP
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP), INTENT(IN) :: m
            REAL(DP)             :: f
          END FUNCTION f
       END INTERFACE
       INTERFACE
          FUNCTION fp(x,m)
            use nrtype, only : DP
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP), INTENT(IN) :: m
            REAL(DP)             :: fp
          END FUNCTION fp
       END INTERFACE
       INTERFACE
          FUNCTION fpp(x,m)
            use nrtype, only : DP
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP), INTENT(IN) :: m
            REAL(DP)             :: fpp
          END FUNCTION fpp
       END INTERFACE
       INTERFACE
          FUNCTION fppp(x,m)
            use nrtype, only : DP
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP), INTENT(IN) :: m
            REAL(DP)             :: fppp
          END FUNCTION fppp
       END INTERFACE
     END SUBROUTINE splint_horner3_a
  END INTERFACE

  interface splint_horner1
    subroutine splint_horner1_a(xa, a, b, c, d, swd, m, x_in, &
        & f, fp, fpp, fppp, y, yp, ypp, yppp)
      use nrtype, only : I4B, DP

      integer(I4B),               intent(in)  :: swd
      real(DP),                   intent(in)  :: m
      real(DP),     dimension(:), intent(in)  :: xa, a, b, c, d
      real(DP),                   intent(in)  :: x_in
      real(DP),                   intent(out) :: y, yp, ypp, yppp
      interface
        function f(x,m)
          use nrtype, only : DP
          implicit none
          real(DP), intent(in) :: x
          real(DP), intent(in) :: m
          real(DP)             :: f
        end function f
      end interface
      interface
        function fp(x,m)
          use nrtype, only : DP
          implicit none
          real(DP), intent(in) :: x
          real(DP), intent(in) :: m
          real(DP)             :: fp
        end function fp
      end interface
      interface
        function fpp(x,m)
          use nrtype, only : DP
          implicit none
          real(DP), intent(in) :: x
          real(DP), intent(in) :: m
          real(DP)             :: fpp
        end function fpp
      end interface
      interface
        function fppp(x,m)
          use nrtype, only : DP
          implicit none
          real(DP), intent(in) :: x
          real(DP), intent(in) :: m
          real(DP)             :: fppp
        end function fppp
      end interface
    end subroutine splint_horner1_a
  end interface

  INTERFACE splint_horner3_driv_s
     SUBROUTINE splint_horner3_driv_s_a(svec,a,b,c,d,swd,ixm,ixn,s,theta,phi,&
          f,fp,fpp,fppp,y,ys,yt,yp)
       use nrtype, only : I4B, DP
       REAL(DP),     DIMENSION(:),   INTENT(IN)  :: svec
       REAL(DP),     DIMENSION(:,:), INTENT(IN)  :: a, b, c, d
       INTEGER(I4B),                 INTENT(IN)  :: swd
       REAL(DP),     DIMENSION(:),   INTENT(IN)  :: ixm
       INTEGER(I4B), DIMENSION(:),   INTENT(IN)  :: ixn
       REAL(DP),                     INTENT(IN)  :: s, theta, phi
       REAL(DP),                     INTENT(OUT) :: y, ys, yt, yp
       INTERFACE
          FUNCTION f(x,m)
            use nrtype, only : DP
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP), INTENT(IN) :: m
            REAL(DP)             :: f
          END FUNCTION f
       END INTERFACE
       INTERFACE
          FUNCTION fp(x,m)
            use nrtype, only : DP
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP), INTENT(IN) :: m
            REAL(DP)             :: fp
          END FUNCTION fp
       END INTERFACE
       INTERFACE
          FUNCTION fpp(x,m)
            use nrtype, only : DP
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP), INTENT(IN) :: m
            REAL(DP)             :: fpp
          END FUNCTION fpp
       END INTERFACE
       INTERFACE
          FUNCTION fppp(x,m)
            use nrtype, only : DP
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP), INTENT(IN) :: m
            REAL(DP)             :: fppp
          END FUNCTION fppp
       END INTERFACE
     END SUBROUTINE splint_horner3_driv_s_a
  END INTERFACE


  INTERFACE splint_horner3_driv_c
     SUBROUTINE splint_horner3_driv_c_a(svec,a,b,c,d,swd,ixm,ixn,s,theta,phi,&
          f,fp,fpp,fppp,y,ys,yt,yp)
       use nrtype, only : I4B, DP
       REAL(DP),     DIMENSION(:),   INTENT(IN)  :: svec
       REAL(DP),     DIMENSION(:,:), INTENT(IN)  :: a, b, c, d
       INTEGER(I4B),                 INTENT(IN)  :: swd
       REAL(DP),     DIMENSION(:),   INTENT(IN)  :: ixm
       INTEGER(I4B), DIMENSION(:),   INTENT(IN)  :: ixn
       REAL(DP),                     INTENT(IN)  :: s, theta, phi
       REAL(DP),                     INTENT(OUT) :: y, ys, yt, yp
       INTERFACE
          FUNCTION f(x,m)
            use nrtype, only : DP
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP), INTENT(IN) :: m
            REAL(DP)             :: f
          END FUNCTION f
       END INTERFACE
       INTERFACE
          FUNCTION fp(x,m)
            use nrtype, only : DP
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP), INTENT(IN) :: m
            REAL(DP)             :: fp
          END FUNCTION fp
       END INTERFACE
       INTERFACE
          FUNCTION fpp(x,m)
            use nrtype, only : DP
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP), INTENT(IN) :: m
            REAL(DP)             :: fpp
          END FUNCTION fpp
       END INTERFACE
       INTERFACE
          FUNCTION fppp(x,m)
            use nrtype, only : DP
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP), INTENT(IN) :: m
            REAL(DP)             :: fppp
          END FUNCTION fppp
       END INTERFACE
     END SUBROUTINE splint_horner3_driv_c_a
  END INTERFACE


  INTERFACE
     FUNCTION tf(x,m)
       use nrtype, only : DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: x
       REAL(DP), INTENT(IN) :: m
       REAL(DP)             :: tf
     END FUNCTION tf
  END INTERFACE
  INTERFACE
     FUNCTION tfp(x,m)
       use nrtype, only : DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: x
       REAL(DP), INTENT(IN) :: m
       REAL(DP)             :: tfp
     END FUNCTION tfp
  END INTERFACE
  INTERFACE
     FUNCTION tfpp(x,m)
       use nrtype, only : DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: x
       REAL(DP), INTENT(IN) :: m
       REAL(DP)             :: tfpp
     END FUNCTION tfpp
  END INTERFACE
  INTERFACE
     FUNCTION tfppp(x,m)
       use nrtype, only : DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: x
       REAL(DP), INTENT(IN) :: m
       REAL(DP)             :: tfppp
     END FUNCTION tfppp
  END INTERFACE

  INTERFACE
     FUNCTION tfone(x,m)
       use nrtype, only : DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: x
       REAL(DP), INTENT(IN) :: m
       REAL(DP)             :: tfone
     END FUNCTION tfone
  END INTERFACE
  INTERFACE
     FUNCTION tfzero(x,m)
       use nrtype, only : DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: x
       REAL(DP), INTENT(IN) :: m
       REAL(DP)             :: tfzero
     END FUNCTION tfzero
  END INTERFACE

END MODULE inter_interfaces
