! Winny module for ripple added
module ripple_mod
  use nrtype, only : dp

  integer :: ripple_counter
  integer, dimension(:), allocatable :: imax1_ripple
  integer, dimension(:), allocatable :: imax2_ripple
  integer, dimension(:), allocatable :: imin_ripple
  real(kind=dp), dimension(:), allocatable :: col_ripple
  real(kind=dp), dimension(:), allocatable :: eta_x0
  real(kind=dp), dimension(:), allocatable :: eta_s
  real(kind=dp)                            :: bmax_abs,bmin_abs
  integer :: propagator_counter
  real(kind=dp), dimension(:), allocatable :: phi1_prop
  real(kind=dp), dimension(:), allocatable :: phi2_prop
  integer, dimension(:), allocatable :: ibeg_prop,iend_prop
  integer, dimension(:), allocatable :: ripplenumber_prop
  integer :: ibeg
  integer :: iend
  integer :: ibeg_prop_act
  integer :: iend_prop_act
  integer :: imax_eta
end module ripple_mod
