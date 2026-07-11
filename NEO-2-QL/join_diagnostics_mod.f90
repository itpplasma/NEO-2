module join_diagnostics_mod
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    public :: report_join_failure, validate_join_normalization

contains

    subroutine report_join_failure(stage, info, old_tags, new_tags, ndim, ndim1, &
            matrix, rhs, ierr)
        integer, intent(in) :: stage, info, old_tags(2), new_tags(2), ndim, ndim1
        real(kind(1.0d0)), intent(in) :: matrix(:, :), rhs(:, :)
        integer, intent(out) :: ierr

        write(*, '(a,i0)') 'join_ripples: DGBSV stage=', stage
        write(*, '(a,i0)') 'join_ripples: DGBSV info=', info
        write(*, '(a,2(1x,i0))') 'join_ripples: old propagator tags=', old_tags
        write(*, '(a,2(1x,i0))') 'join_ripples: new propagator tags=', new_tags
        write(*, '(a,2(1x,i0))') 'join_ripples: dimensions=', ndim, ndim1
        write(*, '(a,i0,a,i0)') 'join_ripples: finite matrix entries=', &
            count(ieee_is_finite(matrix)), '/', size(matrix)
        write(*, '(a,i0,a,i0)') 'join_ripples: finite rhs entries=', &
            count(ieee_is_finite(rhs)), '/', size(rhs)
        ierr = stage
    end subroutine report_join_failure

    subroutine validate_join_normalization(facnorm, direction, column, old_tags, &
            new_tags, npass_l, npass_r, ierr)
        real(kind(1.0d0)), intent(in) :: facnorm
        character(len=1), intent(in) :: direction
        integer, intent(in) :: column, old_tags(2), new_tags(2), npass_l, npass_r
        integer, intent(out) :: ierr

        if (ieee_is_finite(facnorm) .and. facnorm /= 0.0d0) then
            ierr = 0
            return
        end if

        write(*, '(a,a)') 'join_ripples: normalization direction=', direction
        write(*, '(a,i0)') 'join_ripples: normalization column=', column
        write(*, '(a,es23.15)') 'join_ripples: normalization factor=', facnorm
        write(*, '(a,2(1x,i0))') 'join_ripples: old propagator tags=', old_tags
        write(*, '(a,2(1x,i0))') 'join_ripples: new propagator tags=', new_tags
        write(*, '(a,2(1x,i0))') 'join_ripples: pass dimensions=', npass_l, npass_r
        ierr = 5
    end subroutine validate_join_normalization
end module join_diagnostics_mod
