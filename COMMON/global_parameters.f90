module global_parameters

  implicit none

  ! Public parameters
  public ID_UNKNOWN, ID_PARALLEL, ID_QUASILINEAR
  public DEFAULT_STRING_LENGTH

  ! Public variables
  public  isw_code

  ! Public subroutines
  public init, string_to_id, id_to_string

  ! Public functions
  ! none so far

  integer, parameter :: ID_UNKNOWN = -1
  integer, parameter :: ID_PARALLEL = 0
  integer, parameter :: ID_QUASILINEAR = 1

  integer, parameter :: DEFAULT_STRING_LENGTH = 80

  integer, save :: isw_code

contains

  subroutine init(id_run)
    integer, intent(in) :: id_run

    isw_code = id_run

  end subroutine init

  function string_to_id(str) result(res)
    character(len=DEFAULT_STRING_LENGTH), intent(in) :: str

    integer :: res

    res = ID_UNKNOWN

    if ((trim(str) .eq. 'PARALLEL') .or. (trim(str) .eq. 'ID_PARALLEL')) then
      res = ID_PARALLEL
    elseif ((trim(str) .eq. 'QUASILINEAR') .or. (trim(str) .eq. 'ID_QUASILINEAR')) then
      res = ID_QUASILINEAR
    end if

  end function string_to_id

  function id_to_string(givenid) result(res)
    integer, intent(in) :: givenid

    character(len=DEFAULT_STRING_LENGTH) :: res

    select case (givenid)
    case (ID_PARALLEL)
      res = 'ID_PARALLEL'
    case (ID_QUASILINEAR)
      res = 'ID_QUASILINEAR'
    case (ID_UNKNOWN)
      res = 'ID_UNKNOWN'
    case default
      res = 'ID_UNKNOWN'
    end select
  end function id_to_string

end module global_parameters
