module sparsevec_mod
  use nrtype, only : longint

  implicit none

  private dp
  integer, parameter    :: dp = kind(1.0d0)

  type :: sparsevec
     integer, dimension(:), allocatable :: idxvec
     integer(kind=longint), dimension(:), allocatable :: values
     integer                                          :: len = -1
     integer                                          :: len_sparse = 0

   contains

     procedure :: assign     => assign_sparsevec
     procedure :: modify     => modify_sparsevec
     procedure :: insert     => insert_sparsevec
     procedure :: get        => get_sparsevec
     procedure :: find       => find_sparsevec
     procedure :: reallocate => reallocate_sparsevec
     procedure :: clear      => clear_sparsevec
     procedure :: ibset      => ibset_sparsevec
     procedure :: ibclr      => ibclr_sparsevec
     
  end type sparsevec

contains

  subroutine insert_sparsevec(this, idx, val, idxnear)
    use nrtype, only : longint

    class(sparsevec) :: this
    integer          :: idx, idxnear
    integer(kind=longint) :: val
    logical          :: found
    integer(kind=longint) :: zero

    zero = 0
    if (val .ne. zero) then

    if (this%len_sparse .gt. 0) then

       call this%reallocate(this%len_sparse+1)

       if (idxnear .ne. this%len_sparse) then
          this%idxvec(idxnear+1:this%len_sparse) = this%idxvec(idxnear:this%len_sparse-1) 
          this%values(idxnear+1:this%len_sparse) = this%values(idxnear:this%len_sparse-1) 
       end if
       
       this%idxvec(idxnear) = idx
       this%values(idxnear) = val
    else
       call this%reallocate(this%len_sparse+1)
       this%idxvec(this%len_sparse) = idx
       this%values(this%len_sparse) = val
    end if

    end if
  end subroutine insert_sparsevec
  
  subroutine assign_sparsevec(this, obj)
    class(sparsevec) :: this
    class(sparsevec)  :: obj

    call this%clear()
    if (allocated(obj%idxvec)) then
       if (allocated(this%idxvec)) deallocate(this%idxvec)
       allocate(this%idxvec, source=obj%idxvec)
       !this%idxvec = obj%idxvec

    end if
    if (allocated(obj%values)) then
       if (allocated(this%values)) deallocate(this%values)
       allocate(this%values, source=obj%values)
       !this%values = obj%values

    end if
    this%len        = obj%len
    this%len_sparse = obj%len_sparse

  end subroutine assign_sparsevec
  
  subroutine reallocate_sparsevec(this, len_sparse)
    use nrtype, only : longint

    class(sparsevec)      :: this
    integer               :: len_sparse
    integer, dimension(:), allocatable :: idxvec_copy
    integer(kind=longint), dimension(:), allocatable :: values_copy
    
    if (this%len_sparse .ne. len_sparse) then
       
       if (this%len_sparse .ne. 0) then
          if ((.not. allocated(this%idxvec)) .or. (this%len_sparse .ge. ubound(this%idxvec,1))) then
             allocate(idxvec_copy(this%len_sparse))
             allocate(values_copy(this%len_sparse))

             idxvec_copy = this%idxvec(1:this%len_sparse)
             values_copy = this%values(1:this%len_sparse)

             !**********************************************************
             ! Reallocate more space than needed for performance reasons
             !**********************************************************
             deallocate(this%idxvec)
             deallocate(this%values)
             allocate(this%idxvec(4*len_sparse))
             allocate(this%values(4*len_sparse))
             this%idxvec = 0
             this%values = 0

             this%idxvec(1:this%len_sparse) = idxvec_copy
             this%values(1:this%len_sparse) = values_copy

             deallocate(idxvec_copy)
             deallocate(values_copy)

          end if
       else
          allocate(this%idxvec(len_sparse))
          allocate(this%values(len_sparse))
          !this%idxvec = 0
          !this%values = 0         
       end if
       this%len_sparse = len_sparse

    end if
  end subroutine reallocate_sparsevec

  subroutine clear_sparsevec(this)
    class(sparsevec) :: this

    if (allocated(this%idxvec)) deallocate(this%idxvec)
    if (allocated(this%values)) deallocate(this%values)
    this%len_sparse = 0
  end subroutine clear_sparsevec
  
  function get_sparsevec(this, idx) result(val)
    use nrtype, only : longint

    class(sparsevec)      :: this
    integer               :: idx, idxnear
    integer(kind=longint) :: val
    integer :: k, idx_sparse

    val = 0
    if (allocated(this%idxvec)) then
       !if (.not. any(this%idxvec .eq. idx)) then
       !   val = 0
       !   return
       !end if

       idx_sparse = this%find(idx, idxnear)
       if (idx_sparse .gt. 0) then
          val = this%values(idx_sparse)
       end if

    end if
       
  end function get_sparsevec
  
  subroutine modify_sparsevec(this, idx, val)
    use nrtype, only : longint

    class(sparsevec) :: this
    integer :: idx
    integer(kind=longint) :: val
    integer :: k, idx_sparse, idxnear
    logical :: found

    idx_sparse = this%find(idx, idxnear)
    if (idx_sparse .gt. 0) then
       this%values(idx_sparse) = val
    else
       call this%insert(idx, val, idxnear)
    end if
    
  end subroutine modify_sparsevec


  function find_sparsevec(this, idx, idxnear) result(idx_sparse)
    class(sparsevec) :: this
    integer          :: idx
    integer          :: idx_sparse
    integer, intent(out)         :: idxnear
    integer          :: k
    logical          :: found

    !**********************************************************
    ! Linear find
    !**********************************************************
    !idx_sparse = -1
    !do k = 1,this%len_sparse
    !   if (this%idxvec(k) .eq. idx) then
    !      idx_sparse = k
    !      return
    !   end if
    !end do

    !**********************************************************
    ! Binary find
    !**********************************************************
    found = .false.
    idx_sparse = -1
    if (this%len_sparse .gt. 0 ) then
       call binsrc(this%idxvec(1:this%len_sparse), 1, this%len_sparse, idx, idxnear, found)
       idx_sparse = idxnear - 1
    end if
    if (.not. found) idx_sparse = -1
        
  end function find_sparsevec
  
  subroutine ibset_sparsevec(this, idx, pos)
    use nrtype, only : longint

    class(sparsevec) :: this
    integer          :: idx
    integer          :: pos
    integer          :: idx_sparse, idxnear
    integer(kind=longint) :: zero

    ! This is very important, so that 0 has the correct bitsize
    zero = 0
    
    idx_sparse = this%find(idx, idxnear)
    if (idx_sparse .gt. 0) then
       this%values(idx_sparse) =  ibset(this%values(idx_sparse), pos)
    else
       call this%insert(idx, ibset(zero, pos), idxnear)
    end if
    
  end subroutine ibset_sparsevec

  subroutine ibclr_sparsevec(this, idx, pos)
    class(sparsevec) :: this
    integer          :: idx
    integer          :: pos
    integer          :: idx_sparse, idxnear
    
    idx_sparse = this%find(idx, idxnear)
    if (idx_sparse .gt. 0) then
       this%values(idx_sparse) =  ibclr(this%values(idx_sparse), pos)
    else
       write (*,*) "Warning in ibclr_sparsevec:", idx, "does not exist!"
       stop
    end if
    
  end subroutine ibclr_sparsevec

  subroutine binsrc(p,nmin,nmax,xi,i,found)
    integer                       :: n,nmin,nmax,i,imin,imax,k
    integer                       :: xi
    logical                       :: found
    integer, dimension(nmin:nmax) :: p

    !******************************************************************************
    ! Finds the index  i  of the array of increasing numbers   p  with dimension  n 
    ! which satisfies   p(i-1) <  xi  <  p(i) . Uses binary search algorithm.  
    !******************************************************************************
    imin=nmin
    imax=nmax
    n=nmax-nmin

    !**********************************************************
    ! Some special cases to increase performance
    !**********************************************************
    found = .false.
    if (xi .lt. p(nmin)) then
       i = nmin
       return
    elseif (xi .gt. p(nmax)) then
       i = nmax+1
       return
    elseif (xi .eq. p(nmax)) then
       i = nmax+1
       found = .true.
       return
    end if
    
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

    !**********************************************************
    ! Some special cases for the sparsevec to make insertion easier
    !**********************************************************
    if (nmax - nmin .ge. 1) then
       if ((p(i-1) .eq. xi) .or. (p(i) .eq. xi)) then
          found = .true.
       elseif (p(i) .eq. xi) then
          found = .true.
       end if
    end if
    return
  end subroutine binsrc

end module sparsevec_mod
