module nctools_module

  use netcdf
  implicit none

  character(len=256), public :: nco_path = '/usr/bin/'

  interface nc_define
     module procedure nc_defineMatrix_double
     module procedure nc_defineMatrix_long
     module procedure nc_defineMatrix_int
     module procedure nc_defineArray_double
     module procedure nc_defineArray_int
  end interface nc_define

  interface nc_inquire
     module procedure nc_inquireMatrix_double
     module procedure nc_inquireArray_double
  end interface nc_inquire
  
  
contains

  subroutine nf90_check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       call abort
       stop
    end if
  end subroutine nf90_check

  subroutine nf90_createOrAppend(filename, ncid, exists)
    character(len=*) :: filename
    integer, intent(out) :: ncid
    logical, intent(out) :: exists
    integer :: ierr

    exists = .false.
    ierr = nf90_create(filename, ior(NF90_HDF5, NF90_NOCLOBBER), ncid)
    if (ierr .eq. NF90_EEXIST) then
       call nf90_check(nf90_open(filename, NF90_WRITE, ncid))
       exists = .true.
       !write (*,*) "APPENDING TO ", filename
    else
       call nf90_check(ierr)
    end if
  end subroutine nf90_createOrAppend

  subroutine nc_findGroup(ncid, name, grpid, found)
    integer :: ncid
    character(len=*) :: name
    integer, intent(out) :: grpid
    logical, intent(out) :: found
    character(len=256) :: filename

    found = .true.
    if (nf90_inq_ncid(ncid, name, grpid) /= NF90_NOERR) then
       write(filename,'(100A)') trim(adjustl(name)), '.nc'
       call nf90_check(nf90_open(filename, NF90_NOWRITE, grpid))
       found = .false.
    end if
    
  end subroutine nc_findGroup
  
  subroutine mergeNCFiles(regex, resultfilename)
    character(len=*) :: regex, resultfilename
    
    integer :: status
    character(len=300) :: command

    write (command, *) "find -regex '"// regex //"' -type f -print0 | xargs -r -0 "// &
         &trim(nco_path) // "/ncecat -A --gag -o "// resultfilename //" >> nc.log 2>&1"
    !write (*,*) command
    call system(command, status)

    if (status .eq. 0) then
       write (command, *) "find -regex '"// regex //"' -type f -delete"
       call system(command, status)
    else
       write (*,*) "An error occurred merging the propagator-files. Skipping deletion of files. &
            &See nc.log for more information."
    end if

  end subroutine mergeNCFiles

  subroutine nc_defineMatrix_int(ncid, name, var, varid)
    integer :: ncid
    character(len=*) :: name
    integer, dimension(:,:), allocatable :: var
    integer :: dimid(2)
    integer, intent(out) :: varid

    if (allocated(var)) then
       call nf90_check(nf90_def_dim(ncid, name // "_dim1", size(var,1), dimid(1)))
       call nf90_check(nf90_def_dim(ncid, name // "_dim2", size(var,2), dimid(2)))
       call nf90_check(nf90_def_var(ncid, name, NF90_INT, dimid, varid))

       call nf90_check(nf90_put_att(ncid, varid, 'lbound1', LBOUND(var,1)))
       call nf90_check(nf90_put_att(ncid, varid, 'ubound1', UBOUND(var,1)))
       call nf90_check(nf90_put_att(ncid, varid, 'lbound2', LBOUND(var,2)))
       call nf90_check(nf90_put_att(ncid, varid, 'ubound2', UBOUND(var,2)))
    end if
  end subroutine nc_defineMatrix_int
  
  subroutine nc_defineMatrix_long(ncid, name, var, varid)
    integer :: ncid
    character(len=*) :: name
    integer(kind=8), dimension(:,:), allocatable :: var
    integer :: dimid(2)
    integer, intent(out) :: varid

    if (allocated(var)) then
       call nf90_check(nf90_def_dim(ncid, name // "_dim1", size(var,1), dimid(1)))
       call nf90_check(nf90_def_dim(ncid, name // "_dim2", size(var,2), dimid(2)))
       call nf90_check(nf90_def_var(ncid, name, NF90_INT64, dimid, varid))

       call nf90_check(nf90_put_att(ncid, varid, 'lbound1', LBOUND(var,1)))
       call nf90_check(nf90_put_att(ncid, varid, 'ubound1', UBOUND(var,1)))
       call nf90_check(nf90_put_att(ncid, varid, 'lbound2', LBOUND(var,2)))
       call nf90_check(nf90_put_att(ncid, varid, 'ubound2', UBOUND(var,2)))
    end if
  end subroutine nc_defineMatrix_long

  subroutine nc_defineMatrix_double(ncid, name, var, varid)
    integer :: ncid
    character(len=*) :: name
    double precision, dimension(:,:), allocatable :: var
    integer :: dimid(2)
    integer, intent(out) :: varid

    if (allocated(var)) then
       call nf90_check(nf90_def_dim(ncid, name // "_dim1", size(var,1), dimid(1)))
       call nf90_check(nf90_def_dim(ncid, name // "_dim2", size(var,2), dimid(2)))
       call nf90_check(nf90_def_var(ncid, name, NF90_DOUBLE, dimid, varid))

       call nf90_check(nf90_put_att(ncid, varid, 'lbound1', LBOUND(var,1)))
       call nf90_check(nf90_put_att(ncid, varid, 'ubound1', UBOUND(var,1)))
       call nf90_check(nf90_put_att(ncid, varid, 'lbound2', LBOUND(var,2)))
       call nf90_check(nf90_put_att(ncid, varid, 'ubound2', UBOUND(var,2)))
    end if
    
  end subroutine nc_defineMatrix_double

   subroutine nc_defineArray_double(ncid, name, var, varid)
    integer :: ncid
    character(len=*) :: name
    double precision, dimension(:), allocatable :: var
    integer :: dimid
    integer, intent(out) :: varid

    if (allocated(var)) then
       call nf90_check(nf90_def_dim(ncid, name // "_dim", size(var,1), dimid))
       call nf90_check(nf90_def_var(ncid, name, NF90_DOUBLE, dimid, varid))

       call nf90_check(nf90_put_att(ncid, varid, 'lbound', LBOUND(var)))
       call nf90_check(nf90_put_att(ncid, varid, 'ubound', UBOUND(var)))
    end if
    
  end subroutine nc_defineArray_double

   subroutine nc_defineArray_int(ncid, name, var, varid)
    integer :: ncid
    character(len=*) :: name
    integer, dimension(:), allocatable :: var
    integer :: dimid
    integer, intent(out) :: varid

    if (allocated(var)) then
       call nf90_check(nf90_def_dim(ncid, name // "_dim", size(var,1), dimid))
       call nf90_check(nf90_def_var(ncid, name, NF90_INT, dimid, varid))

       call nf90_check(nf90_put_att(ncid, varid, 'lbound', LBOUND(var)))
       call nf90_check(nf90_put_att(ncid, varid, 'ubound', UBOUND(var)))
    end if
    
  end subroutine nc_defineArray_int  

  subroutine nc_inquireMatrix_double(ncid, name, varid, lb1, ub1, lb2, ub2)
    integer :: ncid
    character(len=*) :: name
    integer, intent(out) :: varid, lb1, ub1, lb2, ub2

    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_get_att(ncid, varid, "lbound1", lb1))
    call nf90_check(nf90_get_att(ncid, varid, "ubound1", ub1))
    call nf90_check(nf90_get_att(ncid, varid, "lbound2", lb2))
    call nf90_check(nf90_get_att(ncid, varid, "ubound2", ub2))
  end subroutine nc_inquireMatrix_double
  
  subroutine nc_inquireArray_double(ncid, name, varid, lb, ub)
    integer :: ncid
    character(len=*) :: name
    integer, intent(out) :: varid, lb, ub

    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_get_att(ncid, varid, "lbound", lb))
    call nf90_check(nf90_get_att(ncid, varid, "ubound", ub))
  end subroutine nc_inquireArray_double
  
end module nctools_module
