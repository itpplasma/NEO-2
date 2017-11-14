module sparsevec_mod

  implicit none
  
  private longint
  include 'longint.f90'

  private dp
  integer, parameter    :: dp = kind(1.0d0)

  type :: sparsevec
     integer, dimension(:), allocatable       :: idxvec
     integer(kind=longint), dimension(:), allocatable :: values
     integer                                  :: len = 0
     integer                                  :: len_sparse = 0

   contains

     procedure :: assign     => assign_sparsevec
     procedure :: modify     => modify_sparsevec
     procedure :: get        => get_sparsevec
     procedure :: find       => find_sparsevec
     procedure :: reallocate => reallocate_sparsevec
     procedure :: clear      => clear_sparsevec
     procedure :: ibset      => ibset_sparsevec
     procedure :: ibclr      => ibclr_sparsevec
     
  end type sparsevec

contains

  subroutine assign_sparsevec(this, obj)
    class(sparsevec) :: this
    class(sparsevec)  :: obj

    call this%clear()
    if (allocated(obj%idxvec)) this%idxvec = obj%idxvec
    if (allocated(obj%values)) this%values = obj%values
    this%len        = obj%len
    this%len_sparse = obj%len_sparse

  end subroutine assign_sparsevec
  
  subroutine reallocate_sparsevec(this, len_sparse)
    class(sparsevec)      :: this
    integer               :: len_sparse
    integer, dimension(:), allocatable :: idxvec_copy
    integer(kind=longint), dimension(:), allocatable :: values_copy
    
    if (this%len_sparse .ne. len_sparse) then
       
       if (this%len_sparse .ne. 0) then
          allocate(idxvec_copy(this%len_sparse))
          allocate(values_copy(this%len_sparse))
   
          idxvec_copy = this%idxvec
          values_copy = this%values

         ! write (*,*) "Reallocate", this%len_sparse, len_sparse, ubound(this%idxvec), ubound(idxvec_copy)
          
          deallocate(this%idxvec)
          deallocate(this%values)
          allocate(this%idxvec(len_sparse))
          allocate(this%values(len_sparse))
          this%idxvec = 0
          this%values = 0
         
          this%idxvec(1:this%len_sparse) = idxvec_copy
          this%values(1:this%len_sparse) = values_copy

          deallocate(idxvec_copy)
          deallocate(values_copy)
       else
          allocate(this%idxvec(len_sparse))
          allocate(this%values(len_sparse))
          this%idxvec = -1
          this%values = -1         
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
    class(sparsevec)      :: this
    integer               :: idx
    integer(kind=longint) :: val
    integer :: k

    val = 0
    do k = 1,this%len_sparse
       if (this%idxvec(k) .eq. idx) then
          val = this%values(k)
          return
       end if
    end do
       
  end function get_sparsevec
  
  subroutine modify_sparsevec(this, idx, val)
    class(sparsevec) :: this
    integer :: idx
    integer(kind=longint) :: val
    integer :: k
    logical :: found
    
    found = .false.
    do k = 1,this%len_sparse
       if (this%idxvec(k) .eq. idx) then
          this%values(k) = val
          found = .true.
          !write (*,*) "Found", k, idx, val
          !return
       end if
    end do

    if (.not. found) then
       !write (*,*) "Autoreallocate to", this%len_sparse+1, "because ", idx, "was not found in ", &
       !     this%idxvec, "and values are", this%values
       call this%reallocate(this%len_sparse+1)
       this%idxvec(this%len_sparse) = idx
       this%values(this%len_sparse) = val
       !read (*,*) 
    end if

    
    !if (this%len_sparse .ge. (ubound(this%values,1) - lbound(this%values,1))) then
    !   write (*,*) "Sparsevec: Not well initialized. Stopping."
    !   call abort()
    !end if
    
  end subroutine modify_sparsevec


  function find_sparsevec(this, idx) result(idx_sparse)
    class(sparsevec) :: this
    integer          :: idx
    integer          :: idx_sparse
    integer          :: k
    
    idx_sparse = -1
    do k = 1,this%len_sparse
       if (this%idxvec(k) .eq. idx) then
          idx_sparse = k
          return
       end if
    end do
  end function find_sparsevec
  
  subroutine ibset_sparsevec(this, idx, pos)
    class(sparsevec) :: this
    integer          :: idx
    integer          :: pos
    integer          :: idx_sparse
    
    idx_sparse = this%find(idx)
    if (idx_sparse .ge. 0) then
       this%values(idx_sparse) =  ibset(this%values(idx_sparse), pos)
    else
       call this%reallocate(this%len_sparse+1)
       this%idxvec(this%len_sparse) = idx
       this%values(this%len_sparse) = ibset(0, pos)
       !write (*,*) "Error in Sparsevec:", idx, "does not exist!"
       !call abort()
    end if
    
  end subroutine ibset_sparsevec

  subroutine ibclr_sparsevec(this, idx, pos)
    class(sparsevec) :: this
    integer          :: idx
    integer          :: pos
    integer          :: idx_sparse
    
    idx_sparse = this%find(idx)
    if (idx_sparse .ge. 0) then
       this%values(idx_sparse) =  ibclr(this%values(idx_sparse), pos)
    else
       !call this%reallocate(this%len_sparse+1)
       !this%idxvec(this%len_sparse) = idx
       !this%values(this%len_sparse) = ibset(0, pos)
       write (*,*) "Warning in ibclr_sparsevec:", idx, "does not exist!"
       !call abort()
    end if
    
  end subroutine ibclr_sparsevec
  
end module sparsevec_mod
