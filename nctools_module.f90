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
     module procedure nc_defineScalar_int
     module procedure nc_defineScalar_double
     module procedure nc_defineMultiDim
     module procedure nc_defineString
  end interface nc_define

  interface nc_defineUnlimited
     module procedure nc_defineUnlimitedArray
  end interface nc_defineUnlimited

  interface nc_inquire
     module procedure nc_inquireMatrix_double
     module procedure nc_inquireArray_double
  end interface nc_inquire

  interface nc_quickAdd
     module procedure nc_quickAddScalar_int
     module procedure nc_quickAddScalar_double
     module procedure nc_quickAddMatrix_double
     module procedure nc_quickAddArray_double
     module procedure nc_quickAddString
  end interface nc_quickAdd

  interface nc_quickGet
     module procedure nc_quickGetScalar_int
     module procedure nc_quickGetScalar_double
  end interface nc_quickGet

contains

  subroutine nc_quickGetScalar_int(ncid, name, var)
    integer :: ncid
    character(len=*) :: name
    integer, intent(out) :: var
    integer :: varid

    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_get_var(ncid, varid, var))
  end subroutine nc_quickGetScalar_int

  subroutine nc_quickGetScalar_double(ncid, name, var)
    integer :: ncid
    character(len=*) :: name
    double precision, intent(out) :: var
    integer :: varid

    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_get_var(ncid, varid, var))
  end subroutine nc_quickGetScalar_double

  subroutine nc_quickAddScalar_int(ncid, name, var, comment, unit)
    integer :: ncid
    character(len=*) :: name
    integer :: var
    character(len=*), optional :: comment, unit
    integer :: ierr, varid

    ierr = nf90_redef(ncid)
    call nc_define(ncid, name, var, varid, comment, unit)
    ierr = nf90_enddef(ncid)
    call nf90_check(nf90_put_var(ncid, varid, var))

  end subroutine nc_quickAddScalar_int

  subroutine nc_quickAddScalar_double(ncid, name, var, comment, unit)
    integer :: ncid
    character(len=*) :: name
    double precision :: var
    character(len=*), optional :: comment, unit
    integer :: ierr, varid

    ierr = nf90_redef(ncid)
    call nc_define(ncid, name, var, varid, comment, unit)
    ierr = nf90_enddef(ncid)
    call nf90_check(nf90_put_var(ncid, varid, var))

  end subroutine nc_quickAddScalar_double

  subroutine nc_quickAddMatrix_double(ncid, name, var, comment, unit)
    integer :: ncid
    character(len=*) :: name
    double precision, dimension(:,:), allocatable :: var
    character(len=*), optional :: comment, unit
    integer :: ierr, varid

    ierr = nf90_redef(ncid)
    call nc_define(ncid, name, var, varid, comment, unit)
    ierr = nf90_enddef(ncid)
    call nf90_check(nf90_put_var(ncid, varid, var))

  end subroutine nc_quickAddMatrix_double

  subroutine nc_quickAddArray_double(ncid, name, var, comment, unit)
    integer :: ncid
    character(len=*) :: name
    !double precision, dimension(:) :: var
    double precision, dimension(:), allocatable :: var
    character(len=*), optional :: comment, unit
    integer :: ierr, varid

    ierr = nf90_redef(ncid)
    call nc_define(ncid, name, var, varid, comment, unit)
    ierr = nf90_enddef(ncid)
    call nf90_check(nf90_put_var(ncid, varid, var))

  end subroutine nc_quickAddArray_double

  subroutine nc_quickAddArrayNonAlloc_double(ncid, name, var, comment, unit)
    integer :: ncid
    character(len=*) :: name
    double precision, dimension(:) :: var
    !double precision, dimension(:), allocatable :: var
    character(len=*), optional :: comment, unit
    integer :: ierr, varid
    
    ierr = nf90_redef(ncid)
    call nc_defineArrayNonAlloc_double(ncid, name, var, varid, comment, unit)
    ierr = nf90_enddef(ncid)
    call nf90_check(nf90_put_var(ncid, varid, var))
    
  end subroutine nc_quickAddArrayNonAlloc_double
  
  subroutine nc_quickAddString(ncid, name, var, comment, unit)
    integer :: ncid
    character(len=*) :: name
    character(len=*) :: var
    character(len=*), optional :: comment, unit
    integer :: ierr, varid

    ierr = nf90_redef(ncid)
    call nc_define(ncid, name, varid, comment, unit)
    ierr = nf90_enddef(ncid)
    call nf90_check(nf90_put_var(ncid, varid, var))

  end subroutine nc_quickAddString

  subroutine nf90_check(status, optException)
    integer, intent ( in) :: status
    logical, optional :: optException
    logical :: exception

    exception = .true.
    if (present(optException)) exception = optException

    if(status /= nf90_noerr) then
       if (exception) then
          print *, trim(nf90_strerror(status))
          call abort
          stop
       end if
    end if
  end subroutine nf90_check

  subroutine nc_create(filename, ncid, fileformat_version)
    character(len=*) :: filename
    integer, intent(out) :: ncid
    character(len=*), optional :: fileformat_version

    if (.not. present(fileformat_version)) then
       fileformat_version = '1.0'
    end if

    write (*,*) "Creating NetCDF-4 File: ", filename
    call nf90_check(nf90_create(filename, NF90_NETCDF4, ncid))
    call nf90_check(nf90_put_att(ncid, NF90_GLOBAL, 'Version', fileformat_version))
  end subroutine nc_create

  subroutine nc_close(ncid)
    integer :: ncid

    call nf90_check(nf90_close(ncid))
  end subroutine nc_close

  subroutine nc_open(filename, ncid, optException)
    character(len=*) :: filename
    integer :: ncid
    logical, optional :: optException
    logical :: exception

    exception = .true.
    if (present(optException)) exception = optException

    call nf90_check(nf90_open(filename, NF90_NOWRITE, ncid), exception)
  end subroutine nc_open

  subroutine nc_enddef(ncid)
    integer :: ncid

    call nf90_check(nf90_enddef(ncid))
  end subroutine nc_enddef

  subroutine nf90_createOrAppend(filename, ncid, exists)
    character(len=*) :: filename
    integer, intent(out) :: ncid
    logical, intent(out) :: exists
    integer :: ierr

    exists = .false.
    !write (*,*) "Create or append ", filename
    ierr = nf90_create(filename, ior(NF90_HDF5, NF90_NOCLOBBER), ncid)
    if (ierr .eq. NF90_EEXIST) then
       call nf90_check(nf90_open(filename, NF90_WRITE, ncid))
       exists = .true.
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
    integer :: ierr

    found = .true.
    ierr = nf90_inq_ncid(ncid, trim(name), grpid);
    if (ierr /= NF90_NOERR) then
       write(filename,'(100A)') trim(adjustl(name)), '.nc'
       call nf90_check(nf90_open(filename, NF90_NOWRITE, grpid))
       found = .false.
    end if

  end subroutine nc_findGroup

  subroutine nc_defineGroup(ncid, grpname, ncid_grp)
    integer :: ncid
    character(len=*) :: grpname
    integer, intent(out) :: ncid_grp

    call nf90_check(nf90_def_grp(ncid, trim(grpname), ncid_grp))

  end subroutine nc_defineGroup

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

  subroutine nc_defineMultiDim(ncid, name, type, dims, varid, comment, unit)
    integer :: ncid
    character(len=*) :: name
    integer :: type
    integer, dimension(:) :: dims
    integer, intent(out) :: varid
    character(len=*), optional :: comment, unit

    integer, dimension(:), allocatable :: dimid
    integer :: i
    character(len=32) :: dimname
    allocate(dimid(size(dims)))

    do i=1,size(dims)
       write (dimname, '(A,A,I1)') name, '_dim', i
       call nf90_check(nf90_def_dim(ncid, dimname, dims(i), dimid(i)))
    end do
    call nf90_check(nf90_def_var(ncid, name, type, dimid, varid))

    if (present(comment)) then
       call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
    end if

    if (present(unit)) then
       call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
    end if

    deallocate(dimid)

  end subroutine nc_defineMultiDim

  subroutine nc_defineMatrix_int(ncid, name, var, varid, comment, unit)
    integer :: ncid
    character(len=*) :: name
    integer, dimension(:,:), allocatable :: var
    integer :: dimid(2)
    integer, intent(out) :: varid
    character(len=*), optional :: comment, unit

    if (allocated(var)) then
       call nf90_check(nf90_def_dim(ncid, name // "_dim1", size(var,1), dimid(1)))
       call nf90_check(nf90_def_dim(ncid, name // "_dim2", size(var,2), dimid(2)))
       call nf90_check(nf90_def_var(ncid, name, NF90_INT, dimid, varid))

       call nf90_check(nf90_put_att(ncid, varid, 'lbound1', lbound(var,1)))
       call nf90_check(nf90_put_att(ncid, varid, 'ubound1', UBOUND(var,1)))
       call nf90_check(nf90_put_att(ncid, varid, 'lbound2', LBOUND(var,2)))
       call nf90_check(nf90_put_att(ncid, varid, 'ubound2', UBOUND(var,2)))
       if (present(comment)) then
          call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
       end if
       if (present(unit)) then
          call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
       end if
    end if
  end subroutine nc_defineMatrix_int

  subroutine nc_defineMatrix_long(ncid, name, var, varid, comment, unit)
    integer :: ncid
    character(len=*) :: name
    integer(kind=8), dimension(:,:), allocatable :: var
    integer :: dimid(2)
    integer, intent(out) :: varid
    character(len=*), optional :: comment, unit

    if (allocated(var)) then
       call nf90_check(nf90_def_dim(ncid, name // "_dim1", size(var,1), dimid(1)))
       call nf90_check(nf90_def_dim(ncid, name // "_dim2", size(var,2), dimid(2)))
       call nf90_check(nf90_def_var(ncid, name, NF90_INT64, dimid, varid))

       call nf90_check(nf90_put_att(ncid, varid, 'lbound1', LBOUND(var,1)))
       call nf90_check(nf90_put_att(ncid, varid, 'ubound1', UBOUND(var,1)))
       call nf90_check(nf90_put_att(ncid, varid, 'lbound2', LBOUND(var,2)))
       call nf90_check(nf90_put_att(ncid, varid, 'ubound2', UBOUND(var,2)))
       if (present(comment)) then
          call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
       end if
       if (present(unit)) then
          call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
       end if
    end if
  end subroutine nc_defineMatrix_long

  subroutine nc_defineMatrix_double(ncid, name, var, varid, comment, unit)
    integer :: ncid
    character(len=*) :: name
    double precision, dimension(:,:), allocatable :: var
    integer :: dimid(2)
    integer, intent(out) :: varid
    character(len=*), optional :: comment, unit

    if (allocated(var)) then
       call nf90_check(nf90_def_dim(ncid, name // "_dim1", size(var,1), dimid(1)))
       call nf90_check(nf90_def_dim(ncid, name // "_dim2", size(var,2), dimid(2)))
       call nf90_check(nf90_def_var(ncid, name, NF90_DOUBLE, dimid, varid))

       call nf90_check(nf90_put_att(ncid, varid, 'lbound1', LBOUND(var,1)))
       call nf90_check(nf90_put_att(ncid, varid, 'ubound1', UBOUND(var,1)))
       call nf90_check(nf90_put_att(ncid, varid, 'lbound2', LBOUND(var,2)))
       call nf90_check(nf90_put_att(ncid, varid, 'ubound2', UBOUND(var,2)))
       if (present(comment)) then
          call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
       end if
       if (present(unit)) then
          call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
       end if
    end if
  end subroutine nc_defineMatrix_double
  
  subroutine nc_defineArray_double(ncid, name, var, varid, comment, unit)
    integer :: ncid
    character(len=*) :: name
    double precision, dimension(:), allocatable :: var
    !double precision, dimension(:) :: var
    integer :: dimid
    integer, intent(out) :: varid
    character(len=*), optional :: comment, unit

    if (allocated(var)) then
       call nf90_check(nf90_def_dim(ncid, name // "_dim", size(var,1), dimid))
       call nf90_check(nf90_def_var(ncid, name, NF90_DOUBLE, dimid, varid))

       call nf90_check(nf90_put_att(ncid, varid, 'lbound', LBOUND(var)))
       call nf90_check(nf90_put_att(ncid, varid, 'ubound', UBOUND(var)))
       if (present(comment)) then
          call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
       end if
       if (present(unit)) then
          call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
       end if
    end if
  end subroutine nc_defineArray_double

  subroutine nc_defineArrayNonAlloc_double(ncid, name, var, varid, comment, unit)
    integer :: ncid
    character(len=*) :: name
    !double precision, dimension(:), allocatable :: var
    double precision, dimension(:) :: var
    integer :: dimid
    integer, intent(out) :: varid
    character(len=*), optional :: comment, unit
    
    !if associated(var) then
    call nf90_check(nf90_def_dim(ncid, name // "_dim", size(var,1), dimid))
    call nf90_check(nf90_def_var(ncid, name, NF90_DOUBLE, dimid, varid))
    
    call nf90_check(nf90_put_att(ncid, varid, 'lbound', lbound(var)))
    call nf90_check(nf90_put_att(ncid, varid, 'ubound', ubound(var)))
    if (present(comment)) then
       call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
    end if
    if (present(unit)) then
       call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
    end if
    !end if
  end subroutine nc_defineArrayNonAlloc_double
  
  subroutine nc_defineArray_int(ncid, name, var, varid, comment, unit)
    integer :: ncid
    character(len=*) :: name
    integer, dimension(:), allocatable :: var
    integer :: dimid
    integer, intent(out) :: varid
    character(len=*), optional :: comment, unit

    if (allocated(var)) then
       call nf90_check(nf90_def_dim(ncid, name // "_dim", size(var,1), dimid))
       call nf90_check(nf90_def_var(ncid, name, NF90_INT, dimid, varid))

       call nf90_check(nf90_put_att(ncid, varid, 'lbound', LBOUND(var)))
       call nf90_check(nf90_put_att(ncid, varid, 'ubound', UBOUND(var)))
       if (present(comment)) then
          call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
       end if
       if (present(unit)) then
          call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
       end if
    end if
  end subroutine nc_defineArray_int

   subroutine nc_defineScalar_int(ncid, name, var, varid, comment, unit)
    integer :: ncid
    character(len=*) :: name
    integer :: var
    integer, intent(out) :: varid
    character(len=*), optional :: comment, unit

    call nf90_check(nf90_def_var(ncid, name, NF90_INT, varid = varid))

    if (present(comment)) then
       call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
    end if
    if (present(unit)) then
       call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
    end if
  end subroutine nc_defineScalar_int

   subroutine nc_defineScalar_double(ncid, name, var, varid, comment, unit)
    integer :: ncid
    character(len=*) :: name
    double precision :: var
    integer, intent(out) :: varid
    character(len=*), optional :: comment, unit

    call nf90_check(nf90_def_var(ncid, name, NF90_DOUBLE, varid = varid))

    if (present(comment)) then
       call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
    end if
    if (present(unit)) then
       call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
    end if
  end subroutine nc_defineScalar_double

  subroutine nc_defineUnlimitedArray(ncid, name, type, varid, comment, unit)
    integer :: ncid
    character(len=*) :: name
    integer :: type
    integer, intent(out) :: varid
    character(len=*), optional :: comment, unit

    integer :: dimid

    call nf90_check(nf90_def_dim(ncid, name // '_dim', NF90_UNLIMITED, dimid))
    call nf90_check(nf90_def_var(ncid, name, type, dimid, varid))
    if (present(comment)) then
       call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
    end if
    if (present(unit)) then
       call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
    end if
  end subroutine nc_defineUnlimitedArray

  subroutine nc_defineString(ncid, name, varid, comment, unit)
    integer :: ncid
    character(len=*) :: name
    integer, intent(out) :: varid
    character(len=*), optional :: comment, unit

    call nf90_check(nf90_def_var(ncid, name, NF90_STRING, varid = varid))

    if (present(comment)) then
       call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
    end if
    if (present(unit)) then
       call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
    end if
  end subroutine nc_defineString

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
