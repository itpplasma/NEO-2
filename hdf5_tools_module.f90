module hdf5_tools_module

  use hdf5
  use h5lt
  use ISO_C_BINDING
  implicit none

  integer :: h5error

  interface h5_add
     module procedure h5_add_int
     module procedure h5_add_int_1
     module procedure h5_add_int8_2
     module procedure h5_add_double_0
     module procedure h5_add_double_1
     module procedure h5_add_double_2
     module procedure h5_add_string
  end interface h5_add

  interface h5_get
     module procedure h5_get_int
     module procedure h5_get_int_1
     module procedure h5_get_int8_2
     module procedure h5_get_double_0
     module procedure h5_get_double_1
     module procedure h5_get_double_2
     module procedure h5_get_double_4
  end interface h5_get

  interface h5_get_bounds
     module procedure h5_get_bounds_1
     module procedure h5_get_bounds_2
  end interface h5_get_bounds

  interface h5_append
     module procedure h5_append_int_0
     module procedure h5_append_double_0
     module procedure h5_append_double_4
  end interface h5_append

  interface h5_define_unlimited
     module procedure h5_define_unlimited_array
     module procedure h5_define_unlimited_matrix
  end interface h5_define_unlimited
  
contains   


  subroutine h5_init()
    call h5open_f(h5error)
    call h5eset_auto_f(1, h5error)
  end subroutine h5_init

  subroutine h5_deinit()
    call h5close_f(h5error)
  end subroutine h5_deinit

  subroutine h5_check()
    if (h5error < 0) then
       write (*,*) "HDF5 Error"
       call H5Eprint_f(h5error)
       call abort()
    end if
  end subroutine h5_check

  subroutine h5_create(filename, h5id, fileformat_version)
    character(len=*)           :: filename
    integer(HID_T)             :: h5id
    character(len=*), optional :: fileformat_version

    if (.not. present(fileformat_version)) then
       fileformat_version = '1.0'
    end if

    write (*,*) "Creating HDF-5 File: ", trim(filename)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, h5id, h5error)
    call h5ltmake_dataset_string_f(h5id, 'version', fileformat_version, h5error)

  end subroutine h5_create

  subroutine h5_close(h5id)
    integer(HID_T)             :: h5id 

    call h5fclose_f(h5id, h5error)
  end subroutine h5_close

  subroutine h5_open(filename, h5id)
    character(len=*)      :: filename
    integer(HID_T)        :: h5id
    
    write (*,*) "Opening HDF-5 File: ", trim(filename)
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, h5id, h5error)    
  end subroutine h5_open
  
  subroutine h5_open_rw(filename, h5id)
    character(len=*)      :: filename
    integer(HID_T)        :: h5id
    
    write (*,*) "Opening HDF-5 File: ", filename
    call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, h5id, h5error)    
  end subroutine h5_open_rw

  subroutine h5_define_group(h5id, grpname, h5grpid)
    integer(HID_T)       :: h5id
    character(len=*)     :: grpname
    integer(HID_T)       :: h5grpid

    call h5gcreate_f(h5id, grpname, h5grpid, h5error)
    call h5_check()
  end subroutine h5_define_group

  subroutine h5_open_group(h5id, grpname, h5id_grp)
    integer(HID_T)       :: h5id
    character(len=*)     :: grpname
    integer(HID_T)       :: h5id_grp

    !write (*,*) "Opening group ", grpname, "."
    call h5gopen_f(h5id, trim(grpname), h5id_grp, h5error)
    call h5_check()
  end subroutine h5_open_group

  subroutine h5_close_group(h5id_grp)
    integer(HID_T) :: h5id_grp

    call h5gclose_f(h5id_grp, h5error)
  end subroutine h5_close_group
  
  subroutine h5_get_nmembers(h5id, grpname, nmembers)
    integer(HID_T)       :: h5id
    character(len=*)     :: grpname
    integer              :: nmembers

    call h5gn_members_f(h5id, grpname, nmembers, h5error)
  end subroutine h5_get_nmembers

  subroutine h5_get_objinfo(h5id, grpname, idx, objname, type)
    integer(HID_T)       :: h5id
    character(len=*)     :: grpname, objname
    integer              :: idx
    integer              :: type
    
    call h5gget_obj_info_idx_f(h5id, grpname , idx, objname, type, &
         h5error)
  end subroutine h5_get_objinfo
  
  subroutine h5_define_unlimited_matrix(h5id, dataset, datatype, dims, dsetid)
    integer(HID_T)         :: h5id
    character(len=*)       :: dataset
    integer, dimension(:)  :: dims
    integer(HID_T)         :: datatype
    integer(HID_T)         :: dsetid

    integer                :: rank
    integer(SIZE_T), dimension(:), allocatable :: maxdims, startdims
    integer(HID_T)         :: dspaceid
    integer(HID_T)         :: crp_list 
    integer                :: k
    
    rank = size(dims,1)
    allocate(maxdims(rank), startdims(rank))
    maxdims = dims
    startdims = dims
    
    do k = lbound(maxdims,1), ubound(maxdims,1)
       if (maxdims(k) == -1) maxdims(k) = H5S_UNLIMITED_F
       if (maxdims(k) == -1) startdims(k)  = 1
    end do
    write (*,*) "Defining chunk: ", startdims
    call h5screate_simple_f(rank, startdims, dspaceid, h5error, maxdims)
    call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, h5error)
    call h5pset_chunk_f(crp_list, rank, startdims, h5error)
    call h5dcreate_f(h5id, dataset, datatype, dspaceid, dsetid, h5error, crp_list )
    call h5sclose_f(dspaceid, h5error)
    call h5pclose_f(crp_list, h5error)
    
    deallocate(maxdims)
    deallocate(startdims)
    
  end subroutine h5_define_unlimited_matrix
  
  subroutine h5_define_unlimited_array(h5id, dataset, datatype, dsetid)
    integer(HID_T)       :: h5id
    character(len=*)     :: dataset
    integer(HID_T)       :: datatype
    integer(HID_T)       :: dsetid

    integer(HID_T)        :: dspaceid
    integer(HID_T)        :: crp_list 
    integer               :: rank = 1
    integer(SIZE_T), dimension(1) :: dims = (/1/)
    integer(SIZE_T), dimension(1) :: maxdims

    maxdims = (/H5S_UNLIMITED_F/)

    call h5screate_simple_f(rank, dims, dspaceid, h5error, maxdims)
    call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, h5error)
    call h5pset_chunk_f(crp_list, rank, dims, h5error)
    CALL h5dcreate_f(h5id, dataset, datatype, dspaceid, &
         dsetid, h5error, crp_list )
    call h5sclose_f(dspaceid, h5error)
    call h5pclose_f(crp_list, h5error)
    
  end subroutine h5_define_unlimited_array

  subroutine h5_get_bounds_1(h5id, dataset, lb1, ub1, unlimited)
    integer(HID_T)         :: h5id
    character(len=*)       :: dataset
    integer, dimension(1)  :: lbounds, ubounds
    integer, intent(inout) :: lb1, ub1
    logical, optional      :: unlimited

    if (present(unlimited)) unlimited = .false.
    call h5eset_auto_f(0, h5error)
    call h5ltget_attribute_int_f(h5id, dataset,'lbounds', lbounds(1), h5error)    
    call h5ltget_attribute_int_f(h5id, dataset,'ubounds', ubounds(1), h5error)
    call h5eset_auto_f(1, h5error)

    if (h5error .gt. 0) then
       lb1 = lbounds(1)
       ub1 = ubounds(1)
    else
       if (present(unlimited)) unlimited = .true.
    end if
    !write (*,*) "get_bounds: ", dataset, lbounds, ubounds

  end subroutine h5_get_bounds_1
  
  subroutine h5_get_bounds_2(h5id, dataset, lb1, lb2, ub1, ub2)
    integer(HID_T)         :: h5id
    character(len=*)       :: dataset
    integer, dimension(1:2):: lbounds, ubounds
    integer, intent(out)   :: lb1, lb2, ub1, ub2

    call h5ltget_attribute_int_f(h5id, dataset,'lbounds', lbounds, h5error)
    call h5_check()

    call h5ltget_attribute_int_f(h5id, dataset,'ubounds', ubounds, h5error)
    call h5_check()

    lb1 = lbounds(1)
    lb2 = lbounds(2)
    ub1 = ubounds(1)
    ub2 = ubounds(2)
    
    !write (*,*) "get_bounds: ", dataset, lbounds, ubounds

  end subroutine h5_get_bounds_2

  subroutine h5_append_int_0(dsetid, value, offset)
    integer(HID_T)       :: dsetid
    integer              :: value
    integer              :: offset
    
    integer(SIZE_T), dimension(1) :: dims = (/1/)
    integer(SIZE_T), dimension(1) :: size
    integer(HID_T)                :: memspace
    integer                       :: rank = 1
    integer(HID_T)                :: dspaceid

    integer(HSIZE_T), dimension(1) :: offsetd
    integer(HSIZE_T), dimension(1) :: countd
    
    size    = (/offset/)
    offsetd = (/offset-1/)
    countd  = (/1/)
    
    call h5dset_extent_f(dsetid, size, h5error)
    call h5_check()
    call h5screate_simple_f(rank, dims, memspace, h5error)
    call h5_check()
    call h5dget_space_f(dsetid, dspaceid, h5error)
    call h5_check()
    call h5sselect_hyperslab_f(dspaceid, H5S_SELECT_SET_F, offsetd, countd, h5error)
    call h5_check()
    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, value, dims, h5error, memspace, dspaceid)
    call h5_check()

    call h5sclose_f(memspace, h5error)
    call h5sclose_f(dspaceid, h5error)
    
  end subroutine h5_append_int_0

 subroutine h5_append_double_0(dsetid, value, offset)
    integer(HID_T)       :: dsetid
    double precision     :: value
    integer              :: offset
    
    integer(SIZE_T), dimension(1) :: dims = (/1/)
    integer(SIZE_T), dimension(1) :: size
    integer(HID_T)                :: memspace
    integer                       :: rank = 1
    integer(HID_T)                :: dspaceid

    integer(HSIZE_T), dimension(1) :: offsetd
    integer(HSIZE_T), dimension(1) :: countd
    
    size    = (/offset/)
    offsetd = (/offset-1/)
    countd  = (/1/)
    
    call h5dset_extent_f(dsetid, size, h5error)
    call h5_check()
    call h5screate_simple_f(rank, dims, memspace, h5error)
    call h5_check()
    call h5dget_space_f(dsetid, dspaceid, h5error)
    call h5_check()
    call h5sselect_hyperslab_f(dspaceid, H5S_SELECT_SET_F, offsetd, countd, h5error)
    call h5_check()
    call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, value, dims, h5error, memspace, dspaceid)
    call h5_check()

    call h5sclose_f(memspace, h5error)
    call h5sclose_f(dspaceid, h5error)
    
  end subroutine h5_append_double_0

  subroutine h5_append_double_4(dsetid, value, offset)
    integer(HID_T)       :: dsetid
    double precision, dimension(:,:,:) :: value
    integer :: offset
    !integer, dimension(4)                :: offset

    integer(SIZE_T), dimension(4) :: dims
    integer(SIZE_T), dimension(4) :: size
    integer(HID_T)                :: memspace
    integer                       :: rank = 4
    integer(HID_T)                :: dspaceid

    integer(HSIZE_T), dimension(4) :: offsetd

    double precision :: start, end
    double precision, save :: diff = 0
    
    size    = (/shape(value), offset/)
    dims    = (/shape(value), 1/)
    offsetd = (/0, 0, 0, offset-1/)
    !offsetd = offset
    !countd  = shape(value)

    !write (*,*) "Size:", size
    !write (*,*) "Dims:", dims
    !write (*,*) "Offset:", offsetd
    call h5dset_extent_f(dsetid, size, h5error)
    call h5_check()
    call h5screate_simple_f(rank, dims, memspace, h5error)
    call h5_check()
    call h5dget_space_f(dsetid, dspaceid, h5error)
    call h5_check()
    call h5sselect_hyperslab_f(dspaceid, H5S_SELECT_SET_F, offsetd, dims, h5error)
    call h5_check()
    call cpu_time(start)

    call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, value, dims, h5error, memspace, dspaceid)
    call h5_check()
    call cpu_time(end)
    diff = diff + (end - start)
    !write (*,*) "TIME: ", diff

    call h5sclose_f(memspace, h5error)
    call h5sclose_f(dspaceid, h5error)
    
  end subroutine h5_append_double_4
    
  subroutine h5_add_int(h5id, dataset, value, comment, unit)
    integer(HID_T)                 :: h5id
    character(len=*)               :: dataset
    integer                        :: value
    character(len=*), optional     :: comment
    character(len=*), optional     :: unit
    integer(HSIZE_T)               :: dims(1) = (/0/)

    call h5ltmake_dataset_int_f(h5id, dataset, 0,dims, (/value/), h5error)
    if (present(comment)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'comment', comment, h5error)
    end if
    if (present(unit)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'unit', unit, h5error)
    end if
    call h5_check()

  end subroutine h5_add_int

  subroutine h5_add_int_1(h5id, dataset, value, lbounds, ubounds, comment, unit)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    integer, dimension(:)             :: value
    integer, dimension(:)             :: lbounds, ubounds
    character(len=*), optional        :: comment
    character(len=*), optional        :: unit
    integer(HSIZE_T), dimension(:), allocatable    :: dims
    integer(SIZE_T)                   :: size
    integer                           :: rank = 1

    allocate(dims(rank))
    dims = ubounds - lbounds + 1
    size = rank
    call h5ltmake_dataset_int_f(h5id, dataset, rank, dims, value, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbounds, size, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubounds, size, h5error)   
    
    if (present(comment)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'comment', comment, h5error)
    end if
    if (present(unit)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'unit', unit, h5error)
    end if

    deallocate(dims)
    call h5_check()

  end subroutine h5_add_int_1
  
  subroutine h5_add_int8_2(h5id, dataset, value, lbounds, ubounds, comment, unit)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    integer(kind=8), dimension(:,:),target :: value
    integer, dimension(:)             :: lbounds, ubounds
    character(len=*), optional        :: comment
    character(len=*), optional        :: unit
    integer(HSIZE_T), dimension(:), allocatable    :: dims
    integer(SIZE_T)                   :: size
    integer                           :: rank = 2
    integer(HID_T) :: dspace_id, dset_id
    integer(kind=8), dimension(:,:), pointer :: test
    type(C_PTR) :: f_ptr
    integer(HID_T) :: h5_kind_type_i ! HDF type corresponding to the specified KIND

    
    allocate(dims(rank))
    dims = ubounds - lbounds + 1
    size = rank
    h5_kind_type_i = h5kind_to_type(8,H5_INTEGER_KIND)

    call h5screate_simple_f(rank, dims, dspace_id, h5error)
    call h5dcreate_f(h5id, dataset, h5_kind_type_i, dspace_id, dset_id, h5error)

    test => value
    f_ptr = c_loc(test(1,1))
    call h5dwrite_f(dset_id, h5_kind_type_i, f_ptr, h5error)
  
    call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbounds, size, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubounds, size, h5error)   

    if (present(comment)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'comment', comment, h5error)
    end if
    if (present(unit)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'unit', unit, h5error)
    end if

    call h5dclose_f(dset_id, h5error)
    call h5sclose_f(dspace_id, h5error)
    
    deallocate(dims)

    call h5_check()
  end subroutine h5_add_int8_2

   subroutine h5_get_int8_2(h5id, dataset, value)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    integer(kind=8), dimension(:,:), target  :: value
    integer                           :: lb1, lb2, ub1, ub2
    integer(HSIZE_T), dimension(2)    :: dims
    integer(HID_T) :: dspace_id, dset_id
    integer(kind=8), dimension(:,:), pointer :: test
    integer(HID_T) :: h5_kind_type_i ! HDF type corresponding to the specified KIND
    type(C_PTR) :: f_ptr

    h5_kind_type_i = h5kind_to_type(8,H5_INTEGER_KIND)
    
    call h5_get_bounds(h5id, dataset, lb1, lb2, ub1, ub2)
    dims = (/ub1-lb1+1, ub2-lb2+1/)

    call h5dopen_f(h5id, dataset, dset_id, h5error)
    call h5dget_space_f(dset_id, dspace_id, h5error)
    test => value
    f_ptr = c_loc(test(1,1))

    call h5dread_f(dset_id, h5_kind_type_i, f_ptr, h5error)

    call h5dclose_f(dset_id, h5error)
    call h5sclose_f(dspace_id, h5error)

    call h5_check()

  end subroutine h5_get_int8_2
  
  subroutine h5_get_int(h5id, dataset, value)
    integer(HID_T)                 :: h5id
    character(len=*)               :: dataset
    integer, intent(out)           :: value
    integer, dimension(1)          :: buf
    integer(HSIZE_T)               :: dims(1) = (/0/)

    call h5ltread_dataset_int_f_1(h5id, dataset, buf, dims, h5error)
    value = buf(1)

    call h5_check()

  end subroutine h5_get_int

  subroutine h5_get_int_1(h5id, dataset, value)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    integer, dimension(:)             :: value
    integer                           :: lb1, ub1
    integer(HSIZE_T), dimension(1)    :: dims
    
    call h5_get_bounds(h5id, dataset, lb1, ub1)
    dims = (/ub1 - lb1 + 1/)
    !write (*,*) dims, 'UB', ub1, lb1
    call h5ltread_dataset_int_f(h5id, dataset, value, dims, h5error)

    call h5_check()
    
  end subroutine h5_get_int_1
  
  subroutine h5_add_double_0(h5id, dataset, value, comment, unit)
    integer(HID_T)                 :: h5id
    character(len=*)               :: dataset
    double precision               :: value
    character(len=*), optional     :: comment
    character(len=*), optional     :: unit
    integer(HSIZE_T) :: dims(1) = (/0/)

    call h5ltmake_dataset_double_f(h5id, dataset, 0, dims, (/value/), h5error)
    if (present(comment)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'comment', comment, h5error)
    end if
    if (present(unit)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'unit', unit, h5error)
    end if

    call h5_check()
    
  end subroutine h5_add_double_0

  subroutine h5_add_double_1(h5id, dataset, value, lbounds, ubounds, comment, unit)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    double precision, dimension(:)    :: value
    integer, dimension(:)             :: lbounds, ubounds
    character(len=*), optional        :: comment
    character(len=*), optional        :: unit
    integer(HSIZE_T), dimension(:), allocatable    :: dims
    integer(SIZE_T)                   :: size
    integer                           :: rank = 1

    allocate(dims(rank))
    dims = ubounds - lbounds + 1
    size = rank
    call h5ltmake_dataset_double_f(h5id, dataset, rank, dims, value, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbounds, size, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubounds, size, h5error)   
    
    if (present(comment)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'comment', comment, h5error)
    end if
    if (present(unit)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'unit', unit, h5error)
    end if

    deallocate(dims)
    call h5_check()

  end subroutine h5_add_double_1

  subroutine h5_get_double_0(h5id, dataset, value)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    double precision, intent(out)     :: value
    double precision, dimension(1)    :: buf
    integer(HSIZE_T), dimension(1)    :: dims = (/0/)
    
    call h5ltread_dataset_double_f(h5id, dataset, buf, dims, h5error)
    value = buf(1)

    call h5_check()

  end subroutine h5_get_double_0
  
  subroutine h5_get_double_1(h5id, dataset, value)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    double precision, dimension(:)    :: value
    integer                           :: lb1, ub1
    integer(HSIZE_T), dimension(1)    :: dims
    logical                           :: unlimited

    call h5_get_bounds(h5id, dataset, lb1, ub1, unlimited)
    if (.not. unlimited) then
       dims = (/ub1 - lb1 + 1/)
    else
       dims = shape(value)
       !write (*,*) "Unlimited dimension ", dims
    end if
    call h5ltread_dataset_double_f(h5id, dataset, value, dims, h5error)
    call h5_check()
        
  end subroutine h5_get_double_1

   subroutine h5_get_double_2(h5id, dataset, value)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    double precision, dimension(:,:)  :: value
    integer                           :: lb1, lb2, ub1, ub2
    integer(HSIZE_T), dimension(2)    :: dims

    call h5_get_bounds(h5id, dataset, lb1, lb2, ub1, ub2)
    dims = (/ub1-lb1+1, ub2-lb2+1/)
    call h5ltread_dataset_double_f(h5id, dataset, value, dims, h5error)

    call h5_check()
  end subroutine h5_get_double_2

  subroutine h5_get_double_4(h5id, dataset, value)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    double precision, dimension(:,:,:,:)  :: value
    integer(HSIZE_T), dimension(4)    :: dims

    dims = shape(value)
    call h5ltread_dataset_double_f(h5id, dataset, value, dims, h5error)

    call h5_check()
  end subroutine h5_get_double_4
  
  
  subroutine h5_add_double_2(h5id, dataset, value, lbounds, ubounds, comment, unit)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    double precision, dimension(:,:)  :: value
    integer, dimension(:)             :: lbounds, ubounds
    character(len=*), optional        :: comment
    character(len=*), optional        :: unit
    integer(HSIZE_T), dimension(:), allocatable    :: dims
    integer(SIZE_T)                   :: size
    integer                           :: rank = 2

    allocate(dims(rank))
    dims = ubounds - lbounds + 1
    size = rank
    call h5ltmake_dataset_double_f(h5id, dataset, rank, dims, value, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbounds, size, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubounds, size, h5error)   

    if (present(comment)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'comment', comment, h5error)
    end if
    if (present(unit)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'unit', unit, h5error)
    end if

    call h5_check()
  end subroutine h5_add_double_2


  subroutine h5_add_string(h5id, dataset, value, comment, unit)
    integer(HID_T)                 :: h5id
    character(len=*)               :: dataset
    character(len=*)               :: value
    character(len=*), optional     :: comment
    character(len=*), optional     :: unit

    call h5ltmake_dataset_string_f(h5id, dataset, value, h5error)
    if (present(comment)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'comment', comment, h5error)
    end if
    if (present(unit)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'unit', unit, h5error)
    end if

    call h5_check()
  end subroutine h5_add_string


end module hdf5_tools_module
