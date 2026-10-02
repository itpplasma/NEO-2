SUBROUTINE rk4_kin(x,h)

  USE rk4_kin_mod
  USE fieldline_step_mod, ONLY: fieldline_step
  USE mag_interface_mod, ONLY: mag_coordinates
  USE rhs_kin_sub, only: rhs_kin

  IMPLICIT NONE

  INTEGER          :: ierr
  DOUBLE PRECISION :: x,h,length_scale

  IF (mag_coordinates .EQ. 0) THEN
    length_scale = ABS(y(1)) ! R
  ELSE
    length_scale = 1.d0      ! Boozer theta (rad)
  END IF

  call fieldline_step(x,h,y,rhs_kin,length_scale,ierr)
  if (ierr .ne. 0) then
    print *, 'rk4_kin: field line integration failed from phi = ', x, &
         ' to phi = ', x + h
    print *, 'rk4_kin: y at phi = ', y
    error stop 'rk4_kin: field line integration failed'
  end if

END SUBROUTINE rk4_kin
