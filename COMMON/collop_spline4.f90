module collop_spline4
  use gsl_bspline_routines_mod
  use solve_systems
  use inter_precision, only: I4B, DP
  use collisionality_mod, only : phi_x_max

  implicit none
  !integer, parameter :: dp = 8

  logical :: initialized = .false.
  integer :: spl_order = 4
  real(fgsl_double), dimension(:), allocatable :: xknot
  real(kind=dp), dimension(:,:,:), allocatable :: ppcoefs

contains
  subroutine solve_eqsys2(a, b)
    
    implicit none
    
    real(kind=8), dimension(:,:), intent(INOUT) :: a
    real(kind=8), dimension(:,:), intent(INOUT) :: b
    integer(I4B)                                :: info
    integer(I4B) :: i_alloc
    integer(I4B) :: n, nrhs, lda, ldb
    integer(I4B), dimension(:), allocatable :: ipiv
    
    ! --------------------------------------------------------------------
    
    lda  = size(a,1)
    n    = size(a,2)
    ldb  = size(b,1)
    nrhs = size(b,2)
    
    allocate(ipiv(n),  stat = i_alloc)
    if(i_alloc /= 0) stop 'solve_eqsys: Allocation for array failed!'
    call dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
    if (info /= 0) stop 'Matrix singular!'
    
    deallocate(ipiv,  stat = i_alloc)
    if(i_alloc /= 0) stop 'solve_eqsys: Deallocation for array failed!'
    
  end subroutine solve_eqsys2
  
  subroutine binsrc(p,nmin,nmax,xi,i)
    integer                             :: n,nmin,nmax,i,imin,imax,k
    real(kind=dp)                       :: xi
    real(kind=dp), dimension(nmin:nmax) :: p

    !******************************************************************************
    ! Finds the index  i  of the array of increasing numbers   p  with dimension  n 
    ! which satisfies   p(i-1) <  xi  <  p(i) . Uses binary search algorithm.  
    !******************************************************************************
    imin=nmin
    imax=nmax
    n=nmax-nmin

    do k=1,n
       i=(imax-imin)/2+imin
       if(p(i).gt.xi) then
          imax=i
       else
          imin=i
       endif
       if(imax.eq.imin+1) exit
    enddo

    i=imax

    return
  end subroutine binsrc

  subroutine init_phi_spline4(lagmax, legmax)
    integer :: lagmax, legmax
    
    real(fgsl_double), dimension(:),  allocatable :: b
    real(fgsl_double), dimension(:,:),allocatable :: db, b_coef
    integer :: i, j, xi, k, o
    integer :: order, knots
    integer :: nder, nbf
    real(kind=dp) :: dx, xp
    real(kind=dp), dimension(:,:), allocatable :: yinhom
    real(kind=dp), dimension(:,:), allocatable :: Bmat
    real(kind=dp), dimension(:,:), allocatable :: alpha, Pmat
    real(kind=dp), dimension(100) :: work
    integer :: lwork, info

    if (.not. initialized) then
       initialized = .true.
       order = spl_order + 1
       knots = lagmax+1
       nder  = order-1
       nbf   = order + knots - 2

       allocate(Bmat(nbf-3, nbf))
       allocate(xknot(knots))
       allocate(b(knots+order-2))
       allocate(db(nder+1, knots+order-2))
       allocate(b_coef(knots+order-2, knots))
       allocate(alpha(nbf, knots))

       alpha = 0d0
       do i = 1, knots
          xknot(i) = phi_x_max/(knots-1) * (i-1)
          alpha(i,i) = 1!exp(xknot(i)**2)
       end do
       write (*,*) "Knots: ", xknot
       call init_bspline(order, knots)
       call set_bspline_knots(xknot)

       !**********************************************************
       ! Evaluation at knots
       !**********************************************************
       do j = 1, knots
          call bspline_eval(xknot(j), b)
          call bspline_deriv_eval(xknot(j), nder, db)
          Bmat(j,:) = b
       end do

       !**********************************************************
       ! Additional boundary conditions
       !**********************************************************
       !call bspline_deriv_eval(xknot(1), nder, db)
       !Bmat(knots+1,:) = db(3,:)
       !call bspline_deriv_eval(xknot(knots), nder, db)
       !Bmat(knots+2,:) = db(3,:)
       !call bspline_deriv_eval(xknot(knots), nder, db)
       !Bmat(knots+3,:) = db(4,:)
       !call bspline_deriv_eval(xknot(knots), nder, db)
       !Bmat(knots+4,:) = db(4,:)
       
       !**********************************************************
       ! Solve system for alpha
       !**********************************************************
       !call solve_eqsys2(Bmat, alpha)
       lwork = 100
       call dgels('N', size(Bmat,1), size(Bmat,2), size(alpha,2), Bmat, size(Bmat,1), alpha, size(alpha,1), &
            work, lwork, info)
       
       !**********************************************************
       ! Get pp-form for extrapolation
       !**********************************************************
       allocate(yinhom(order, knots))
       allocate(Pmat(order, order))
       allocate(ppcoefs(knots, order, knots))
       do k = 1, knots-1
          dx = (xknot(k+1) - xknot(k))/(order-1) -1d-12
          xp = xknot(k)
          do xi = 1, order
             write (*,*) xp
             call bspline_eval(xp, b)
             yinhom(xi,:) = matmul(b, alpha(1:nbf, :))
             do o = 1, order
                Pmat(xi, o) = (xp-xknot(k))**(order-o)
             end do
             xp = xp + dx
          end do
          call solve_eqsys2(Pmat, yinhom)
          ppcoefs(k,:,:) = yinhom(:,:)
       end do
       
    end if

    xp = 0d0
    do j = 1, 1000
       xp = 2*phi_x_max/999 * (j-1)
       write (545, *) xp, phi_spline4(0,xp), phi_spline4(1,xp), phi_spline4(2,xp), phi_spline4(3,xp), &
            phi_spline4(4,xp)
       write (546, *) xp, d_phi_spline4(0,xp), d_phi_spline4(1,xp), d_phi_spline4(2,xp), d_phi_spline4(3,xp), &
            d_phi_spline4(4,xp)
       write (547, *) xp, dd_phi_spline4(0,xp), dd_phi_spline4(1,xp), dd_phi_spline4(2,xp), dd_phi_spline4(3,xp), &
            dd_phi_spline4(4,xp)
    end do
    !stop
    
  end subroutine init_phi_spline4

  function phi_spline4(m, x) result(phi)
    integer :: m
    real(kind=dp) :: x, phi
    integer :: i, o

    call binsrc(xknot, lbound(xknot,1), ubound(xknot,1), x, i)
    phi = 0d0
    do o = 1, spl_order+1
       phi = phi + ppcoefs(i-1, o, m+1) * (x - xknot(i-1))**(spl_order+1-o)
    end do
    
  end function phi_spline4

  function d_phi_spline4(m, x) result(d_phi)
    integer :: m
    real(kind=dp) :: x, d_phi
    integer :: i, o

    call binsrc(xknot, lbound(xknot,1), ubound(xknot,1), x, i)
    d_phi = 0d0
    do o = 1, spl_order
       d_phi = d_phi + (spl_order+1-o) * ppcoefs(i-1, o, m+1) * (x - xknot(i-1))**(spl_order-o)
    end do
    
  end function d_phi_spline4

  function dd_phi_spline4(m, x) result(dd_phi)
    integer :: m
    real(kind=dp) :: x, dd_phi
    integer :: i, o

    call binsrc(xknot, lbound(xknot,1), ubound(xknot,1), x, i)
    dd_phi = 0d0
    do o = 1, spl_order-1
       dd_phi = dd_phi + (spl_order-o)*(spl_order+1-o) * ppcoefs(i-1, o, m+1) * (x - xknot(i-1))**(spl_order-o-1)
    end do
    
  end function dd_phi_spline4
  
end module collop_spline4
