MODULE device_mod
  ! This module contains only type declarations for the following quantities:
  !
  !   device, surface, fieldline, fieldperiod, fieldpropagator, fieldripple
  !
  ! They are used later to store information in these pointers which in most 
  ! cases are linked to a previous and a next pointer
  !
  ! With this the whole magnetic structure is defined.
  USE magnetics_mod

  IMPLICIT NONE

  TYPE(device_struct),          POINTER, PUBLIC :: device
  TYPE(surface_struct),         POINTER, PUBLIC :: surface
  TYPE(fieldline_struct),       POINTER, PUBLIC :: fieldline
  TYPE(fieldperiod_struct),     POINTER, PUBLIC :: fieldperiod
  TYPE(fieldpropagator_struct), POINTER, PUBLIC :: fieldpropagator
  TYPE(fieldripple_struct),     POINTER, PUBLIC :: fieldripple

END MODULE device_mod
