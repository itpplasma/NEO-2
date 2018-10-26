module gfactor_mod
  use nrtype, only : dp

  implicit none

  integer :: ienter=1,npoia,npoib
  real(kind=dp) :: ha,hb
  real(kind=dp), dimension(:,:), allocatable :: garr

contains

  real(kind=dp) function gfactor(a,b)

    integer :: nistep,i,k,j
    real(kind=dp) :: a,b,atmp,btmp,wa,wb,x,hint
    real(kind=dp), dimension(:), allocatable :: arr

    if (ienter.EQ.1) then
      ienter=0
      npoia=100
      npoib=100
      nistep=100
      ha=1.d0/npoia
      hb=1.d0/(npoib+1)
      hint=1.d0/nistep
      allocate(garr(0:npoia,0:npoib),arr(0:nistep))
      do j=0,npoia
        atmp=j*ha
        do k=0,npoib
          btmp=k*hb
          do i=0,nistep
            x=i*hint
            arr(i)=SQRT(ABS(1.d0-atmp*(3.d0*x**2-2.d0*x**3)))          &
                    /(1.d0-btmp*(3.d0*(x-1.d0)**2+2.d0*(x-1.)**3))
          end do
          garr(j,k)=0.5d0*hint*(SUM(arr(0:nistep-1))+SUM(arr(1:nistep)))
        end do
      end do
      deallocate(arr)
    end if

    wa=a/ha
    j=INT(wa)
    j=MIN(npoia-1,MAX(0,j))
    wa=wa-j

    wb=b/hb
    k=INT(wb)
    k=MIN(npoib-1,MAX(0,k))
    wb=MIN(1.d0,wb-k)

    gfactor=(1.d0-wa)*(garr(j,k)*(1.d0-wb)+garr(j,k+1)*wb)                    &
         +       wa*(garr(j+1,k)*(1.d0-wb)+garr(j+1,k+1)*wb)

  end function gfactor
end module gfactor_mod
