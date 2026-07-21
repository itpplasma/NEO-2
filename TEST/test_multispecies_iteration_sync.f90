PROGRAM test_multispecies_iteration_sync
    USE mpiprovider_module, ONLY : mpro
    USE multispecies_iteration_sync_mod, ONLY : gather_iteration_convergence
    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: changes,thresholds
    DOUBLE PRECISION :: local_change
    INTEGER :: ispec,iter,exit_iteration
    LOGICAL :: converged

    CALL mpro%init()
    IF (mpro%getNumProcs() .NE. 2) &
        ERROR STOP 'FAIL: test requires exactly two MPI ranks'
    ispec=mpro%getRank()
    ALLOCATE(changes(0:mpro%getNumProcs()-1))
    ALLOCATE(thresholds(0:mpro%getNumProcs()-1))
    changes=0.d0
    thresholds=0.d0
    exit_iteration=0

    DO iter=1,5
        IF (iter .GE. ispec+2) THEN
            local_change=0.25d0
        ELSE
            local_change=1.d0
        END IF
        CALL gather_iteration_convergence(local_change,0.5d0,ispec, &
            changes,thresholds,converged)
        IF (converged) THEN
            exit_iteration=iter
            EXIT
        END IF
    END DO

    IF (exit_iteration .NE. 3) &
        ERROR STOP 'FAIL: ranks did not wait for the slow species'
    IF (ANY(changes .GT. thresholds)) &
        ERROR STOP 'FAIL: gathered convergence state is inconsistent'
    IF (mpro%isMaster()) PRINT *, 'All tests passed!'

    DEALLOCATE(changes,thresholds)
    CALL mpro%deinit(.FALSE.)
END PROGRAM test_multispecies_iteration_sync
