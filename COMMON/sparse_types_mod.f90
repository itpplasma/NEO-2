MODULE sparse_types_mod
  ! Module containing basic type definitions and parameters
  ! Extracted from sparse_mod.f90 for better modularity
  
  IMPLICIT NONE
  
  ! Kind parameters
  INTEGER, PARAMETER :: dp = KIND(1.0d0)
  INTEGER, PARAMETER :: long = 8
  
END MODULE sparse_types_mod