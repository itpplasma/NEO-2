module ntv_work_export
    use, intrinsic :: iso_fortran_env, only: dp => real64, int32, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: ntv_work_export_t, work_begin, work_basis, work_vector, work_end
    public :: work_node_deltas

    type :: ntv_work_export_t
        logical :: enabled = .false.
        integer :: text_unit = 0, basis_unit = 0
        integer(int64) :: basis_expected = 0, basis_written = 0
        character(128) :: prefix = ''
    end type

contains

    subroutine work_begin(export, surface, species, nvector, nbasis)
        type(ntv_work_export_t), intent(out) :: export
        real(dp), intent(in) :: surface
        integer, intent(in) :: species, nvector, nbasis
        logical :: exists
        integer :: unit, ios, selected_species
        integer, save :: serial = 0
        real(dp) :: selected_surface

        inquire(file='NEO2_WORK.REQUEST', exist=exists)
        if (.not. exists) return
        open(newunit=unit, file='NEO2_WORK.REQUEST', status='old', &
             action='read', iostat=ios)
        if (ios /= 0) error stop 'Cannot open work request'
        read(unit, *, iostat=ios) selected_surface, selected_species
        close(unit)
        if (ios /= 0) error stop 'Invalid work request'
        if (abs(surface - selected_surface) > 1.0e-10_dp) return
        if (selected_species >= 0) then
            if (species /= selected_species) return
        end if
        if (storage_size(0.0_dp) /= 64) error stop 'Work export requires real64'
        if (storage_size(cmplx(0.0_dp, 0.0_dp, dp)) /= 128) &
            error stop 'Work export requires complex128'
        serial = serial + 1
        write(export%prefix, '(A,I0,A,I0)') 'NEO2_WORK_s', species, '_j', serial
        open(newunit=export%text_unit, &
             file=trim(export%prefix)//'.txt', status='new', action='write')
        write(export%text_unit, '(A)') 'NEO2_COMPLEX_WORK_V1'
        write(export%text_unit, *) 'VECTOR_LENGTH', nvector
        export%basis_expected = int(nbasis, int64)
        open(newunit=export%basis_unit, &
             file=trim(export%prefix)//'_basis.bin', status='new', &
             action='write', access='stream', form='unformatted')
        write(export%basis_unit) 'N2BAS001', int(z'01020304', int32), &
            164_int32, export%basis_expected
        write(export%text_unit, *) 'BASIS_FILE ', trim(export%prefix)//'_basis.bin'
        export%enabled = .true.
    end subroutine

    subroutine work_node_deltas(node, lower, upper, positive, negative, values)
        integer, intent(in) :: node, lower, upper
        real(dp), intent(in) :: positive(lower:upper), negative(lower:upper)
        real(dp), intent(out) :: values(2)

        if (node < lower .or. node > upper) error stop 'Invalid work node'
        ! The native rearrangement never consumes positive(lower)/negative(upper).
        ! Select before reading; do not initialise or mask any solver coefficient.
        values = 0.0_dp
        if (node > lower) values(1) = positive(node)
        if (node < upper) values(2) = negative(node)
    end subroutine

    subroutine work_basis(export, node, band, stencil_start, lambda, eta0, &
                          amin2ovb, coefficients)
        type(ntv_work_export_t), intent(inout) :: export
        integer, intent(in) :: node, band, stencil_start
        real(dp), intent(in) :: lambda, eta0, amin2ovb, coefficients(4, 4)

        if (.not. export%enabled) return
        if (.not. all(ieee_is_finite(coefficients))) &
            error stop 'Nonfinite pitch reconstruction coefficients'
        write(export%basis_unit) int(node, int32), int(band, int32), &
            int(stencil_start, int32), lambda, eta0, amin2ovb, coefficients
        export%basis_written = export%basis_written + 1
    end subroutine

    subroutine work_vector(export, phase, channel, species, values)
        type(ntv_work_export_t), intent(in) :: export
        character(*), intent(in) :: phase
        integer, intent(in) :: channel, species
        complex(dp), intent(in) :: values(:)
        integer :: unit, first, last
        character(64) :: role
        character(256) :: path

        if (.not. export%enabled) return
        if (channel /= 1) then
            if (channel /= 3) return
        end if
        write(role, '(A,A,I0,A,I0)') phase, '_k', channel, '_b', species
        path = trim(export%prefix)//'_'//trim(role)//'.bin'
        open(newunit=unit, file=trim(path), status='new', action='write', &
             access='stream', form='unformatted')
        write(unit) 'N2WORK01', int(z'01020304', int32), 64_int32, 128_int32, &
            int(size(values), int64)
        ! Small chunks also handle a noncontiguous native flux_vector row.
        do first = 1, size(values), 4096
            last = min(size(values), first + 4095)
            if (.not. all(ieee_is_finite(real(values(first:last), dp)))) &
                error stop 'Nonfinite complex work real part'
            if (.not. all(ieee_is_finite(aimag(values(first:last))))) &
                error stop 'Nonfinite complex work imaginary part'
            write(unit) values(first:last)
        end do
        close(unit)
        write(export%text_unit, *) 'VECTOR ', trim(role), size(values), trim(path)
    end subroutine

    subroutine work_end(export)
        type(ntv_work_export_t), intent(inout) :: export

        if (.not. export%enabled) return
        if (export%basis_written /= export%basis_expected) &
            error stop 'Incomplete pitch reconstruction export'
        close(export%basis_unit)
        write(export%text_unit, '(A)') 'END_NEO2_COMPLEX_WORK_V1'
        close(export%text_unit)
        export%enabled = .false.
    end subroutine

end module
