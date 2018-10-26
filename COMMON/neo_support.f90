module neo_support
  use nrtype, only : dp

  interface unit_check
     module procedure unit_check_1
  end interface

  interface strip_extension
     module procedure strip_extension_1
  end interface

  interface add_extension
     module procedure add_extension_1, add_extension_2
  end interface

contains
  !> checks for free unit number
  !> Deprecated due to fortran 2008 newunit specifier?
  subroutine unit_check_1(u)
    implicit none
    integer, intent(inout) :: u
    logical                :: lu
    checku: do
       inquire(unit=u,opened=lu)
       if (.NOT. lu) exit
       u = u + 1
    end do checku
  end subroutine unit_check_1

  subroutine strip_extension_1(str_in,ext,str_out)

    implicit none
    character(len=*), intent(in)    :: str_in
    character(len=*), intent(in)    :: ext
    character(len=*), intent(out)   :: str_out

    integer                         :: ind_ext

    ind_ext = INDEX(str_in,'.'//ext,back=.true.)
    if (ind_ext .NE. 0) then
       str_out = str_in(1:ind_ext-1)
    else
       str_out = str_in
    end if

  end subroutine strip_extension_1

  !> Add extension, in the form '.ext'.
  subroutine add_extension_1(str_in,ext,str_out)

    implicit none
    character(len=*), intent(in)    :: str_in
    character(len=*), intent(in)    :: ext
    character(len=*), intent(out)   :: str_out

    str_out = trim(adjustl(str_in))//'.'//ext

  end subroutine add_extension_1

  !> Add extension in the form '_integer', where integer is passed as an
  !> actual integer (not as a string).
  subroutine add_extension_2(str_in,int,str_out)

    implicit none
    character(len=*), intent(in)    :: str_in
    integer,          intent(in)    :: int
    character(len=*), intent(out)   :: str_out

    character(len=20)               :: ext

    write(ext,*) int
    str_out = trim(adjustl(str_in))//'_'//trim(adjustl(ext))

  end subroutine add_extension_2

end module neo_support
