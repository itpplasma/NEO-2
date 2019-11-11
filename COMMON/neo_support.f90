MODULE neo_support
  USE neo_precision

  INTERFACE unit_check
     MODULE PROCEDURE unit_check_1
  END INTERFACE

  INTERFACE strip_extension
     MODULE PROCEDURE strip_extension_1
  END INTERFACE

  INTERFACE add_extension
     MODULE PROCEDURE add_extension_1, add_extension_2
  END INTERFACE

CONTAINS
  ! checks for free unit number
  SUBROUTINE unit_check_1(u)
    IMPLICIT NONE
    INTEGER, INTENT(inout) :: u
    LOGICAL                :: lu
    checku: DO
       INQUIRE(unit=u,opened=lu)
       IF (.NOT. lu) EXIT
       u = u + 1
    END DO checku
  END SUBROUTINE unit_check_1

  SUBROUTINE strip_extension_1(str_in,ext,str_out)
    USE neo_precision
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in)    :: str_in
    CHARACTER(len=*), INTENT(in)    :: ext
    CHARACTER(len=*), INTENT(out)   :: str_out

    INTEGER                         :: ind_ext

    ind_ext = INDEX(str_in,'.'//ext,back=.TRUE.)
    IF (ind_ext .NE. 0) THEN
       str_out = str_in(1:ind_ext-1)
    ELSE
       str_out = str_in
    END IF

  END SUBROUTINE strip_extension_1

  SUBROUTINE add_extension_1(str_in,ext,str_out)
    USE neo_precision
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in)    :: str_in
    CHARACTER(len=*), INTENT(in)    :: ext
    CHARACTER(len=*), INTENT(out)   :: str_out

    str_out = TRIM(ADJUSTL(str_in))//'.'//ext

  END SUBROUTINE add_extension_1

  SUBROUTINE add_extension_2(str_in,int,str_out)
    USE neo_precision
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in)    :: str_in
    INTEGER,          INTENT(in)    :: int
    CHARACTER(len=*), INTENT(out)   :: str_out

    CHARACTER(len=20)               :: ext

    WRITE(ext,*) int
    str_out = TRIM(ADJUSTL(str_in))//'_'//TRIM(ADJUSTL(ext))

  END SUBROUTINE add_extension_2

END MODULE neo_support
