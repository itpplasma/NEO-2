MODULE multispecies_iteration_sync_mod
    USE mpiprovider_module, ONLY : mpro
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: gather_iteration_convergence

CONTAINS

    ! Every integral_part call contains a world allgather.  Species may have
    ! different local convergence rates, so all ranks must continue the outer
    ! iteration until every species has converged; otherwise a fast species can
    ! enter the next propagator collective while a slow species is still here.
    SUBROUTINE gather_iteration_convergence(local_change,local_threshold, &
            ispec,changes,thresholds,converged)
        DOUBLE PRECISION, INTENT(IN) :: local_change,local_threshold
        INTEGER, INTENT(IN) :: ispec
        DOUBLE PRECISION, DIMENSION(0:), INTENT(INOUT) :: changes,thresholds
        LOGICAL, INTENT(OUT) :: converged

        IF (SIZE(changes) .NE. mpro%getNumProcs() .OR. &
            SIZE(thresholds) .NE. mpro%getNumProcs()) THEN
            ERROR STOP 'multispecies convergence arrays do not match MPI size'
        END IF
        IF (ispec .LT. 0 .OR. ispec .GE. mpro%getNumProcs()) THEN
            ERROR STOP 'multispecies convergence species index is out of range'
        END IF

        changes(ispec)=local_change
        thresholds(ispec)=local_threshold
        CALL mpro%allgather_inplace(changes)
        CALL mpro%allgather_inplace(thresholds)
        converged=ALL(changes .LT. thresholds)
    END SUBROUTINE gather_iteration_convergence

END MODULE multispecies_iteration_sync_mod
