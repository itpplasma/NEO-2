module collop_nbi
  use nrtype, only : dp, pi

  implicit none

  ! --------------------------------------------------------------------
  ! -- namelist variables ----------------------------------------------
  ! --------------------------------------------------------------------
  character(len=80) :: name_nbi_data_file

  integer :: legmax_nbi

  namelist /nbi/ legmax_nbi, name_nbi_data_file
  ! --------------------------------------------------------------------

  logical :: lsw_nbi ! already in namelist collision

contains

  subroutine set_default_values_namelist()
    legmax_nbi = 4

    lsw_nbi = .false.

    name_nbi_data_file = 'nbi_data.dat'
  end subroutine set_default_values_namelist

  subroutine init_phi_nbi(namelist_file_unit)
    integer :: namelist_file_unit

    read(namelist_file_unit,nml=nbi)
    ! Load data

  end subroutine init_phi_nbi

  function phi_nbi(m, x) result(phi)
    integer       :: m
    real(kind=dp) :: x, phi

    ! Return interpolated NBI source data
    
  end function phi_nbi

  function d_phi_nbi(m, x) result(d_phi)
    integer       :: m
    real(kind=dp) :: x, d_phi

    ! Return first derivative of phi_nbi
    
  end function d_phi_nbi

  function dd_phi_nbi(m, x) result(dd_phi)
    integer       :: m
    real(kind=dp) :: x, dd_phi
    
  end function dd_phi_nbi

end module collop_nbi
