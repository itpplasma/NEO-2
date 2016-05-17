!
module prodband_mod
contains
!
  subroutine prodband(nl1,nr1,nhb1,nl2,nr2,nhb2,a1,a2,nhb,a,ierr)
!
! Computes the product of band matrices a1 and a2 with halfband widths
! nhb1 and nhb2, respectively. Puts the result in band matrix a with halfband
! width nhb. 
!
! Input arguments :
!         Formal: nl1,nr1   - actual size of matrix a1
!                 nhb1      - half-band width of matrix a1
!                 nl2,nr2   - actual size of matrix a2
!                 nhb2      - half-band width of matrix a2
!                 a1(n1,n2) - matrix with sizes n1 and n2 .ge. nl1 and nr1, resp
!                 a2(n1,n2) - matrix with sizes n1 and n2 .ge. nl2 and nr2, resp
! Output arguments:
!         Formal: nhb       - half-band width of a product matrix
!                 a(nl1,nr2)- product matrix (allocated with size (nl1,nr2))    
!                 ierr      - error code: 0 - OK, 1 - actual matrix sizes 
!                             disagree
!
  implicit none
!
  integer :: nl1,nr1,nhb1,nl2,nr2,nhb2,nhb,i1,i2,imin2,imax2,imin,imax,ierr
  double precision, dimension(:,:), intent(in)               :: a1,a2
  double precision, dimension(:,:), intent(out), allocatable :: a
!
  if(nr1.eq.nl2) then
    ierr=0
  else
    ierr=1
    return
  endif
!
  nhb=nhb1+nhb2
!
  if(allocated(a)) deallocate(a)
  allocate(a(nl1,nr2))
!
  do i1=1,nl1
    imin2=max(1,i1-nhb)
    imax2=min(nr2,i1+nhb)
    do i2=imin2,imax2
      imin=max(1,i1-nhb1,i2-nhb2)
      imax=min(nr1,i1+nhb1,i2+nhb2)
      a(i1,i2)=sum(a1(i1,imin:imax)*a2(imin:imax,i2))
    enddo
  enddo
!
  return
  end subroutine
!
  subroutine prodband_var(nl1,nr1,nhb1,nl2,nr2,imi,ima,a1,a2,a,ierr)
!
! Computes the product of band matrix a1 with halfband width nhb1 with 
! a matrix a2 which has nonzero elements in columns between imi(i) and ima(i)
! where i is a raw number. Puts the result in matrix a.
!
! Input arguments :
!         Formal: nl1,nr1   - actual size of matrix a1
!                 nhb1      - half-band width of matrix a1
!                 nl2,nr2   - actual size of matrix a2
!                 imi(nr2),ima(nr2) - limits of the range of a2 column index
!                 a1(n1,n2) - matrix with sizes n1 and n2 .ge. nl1 and nr1, resp
!                 a2(n1,n2) - matrix with sizes n1 and n2 .ge. nl2 and nr2, resp
! Output arguments:
!         Formal: a(nl1,nr2)- product matrix (allocated with size (nl1,nr2))    
!                 ierr      - error code: 0 - OK, 1 - actual matrix sizes 
!                             disagree
!
  implicit none
!
  integer :: nl1,nr1,nhb1,nl2,nr2,i1,i2,imin1,imax1,imin,imax,ierr
  integer,          dimension(:),   intent(in)               :: imi,ima
  double precision, dimension(:,:), intent(in)               :: a1,a2
  double precision, dimension(:,:), intent(out), allocatable :: a
!
  if(nr1.eq.nl2) then
    ierr=0
  else
    ierr=1
    return
  endif
!
  if(allocated(a)) deallocate(a)
  allocate(a(nl1,nr2))
!
  do i2=1,nr2
    imin1=max(1,imi(i2)-nhb1)
    imax1=min(nl1,ima(i2)+nhb1)
    do i1=imin1,imax1
      imin=max(1,i1-nhb1,imi(i2))
      imax=min(nr1,i1+nhb1,ima(i2))
      a(i1,i2)=sum(a1(i1,imin:imax)*a2(imin:imax,i2))
    enddo
  enddo
!
  return
  end subroutine
end module
