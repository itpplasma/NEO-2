! Opt-in capture of small scalar systems for an independent solver oracle.
module linear_system_diagnostic_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: diagnostic_system, diagnostic_array
contains
    function diagnostic_path(axisymmetric, species, label) result(path)
        logical, intent(in) :: axisymmetric
        integer, intent(in) :: species
        character(*), intent(in) :: label
        character(:), allocatable :: path
        character(4096) :: directory
        character(32) :: species_text
        integer :: length, status

        call get_environment_variable('NEO2_DUMP_LINEAR_SYSTEM', directory, &
                                      length=length, status=status)
        path = ''
        if (status == 1) return
        if (length == 0) return
        if (status /= 0) error stop 'NEO2_DUMP_LINEAR_SYSTEM path is too long'
        write(species_text, '(i0)') species
        path = trim(directory)//'/'
        if (axisymmetric) then
            path = path//'axi_species'
        else
            path = path//'na_species'
        end if
        path = path//trim(species_text)//'_'//label//'.mtx'
    end function diagnostic_path

    subroutine diagnostic_system(axisymmetric, species, n, rows, cols, values, bra)
        logical, intent(in) :: axisymmetric
        integer, intent(in) :: species, n, rows(:), cols(:)
        complex(dp), intent(in) :: values(:), bra(:, :)
        character(:), allocatable :: path
        integer :: unit, j

        path = diagnostic_path(axisymmetric, species, 'matrix')
        if (len(path) == 0) return
        ! Deliberately bounded: this diagnostic is for small direct systems,
        ! not production-sized matrix exports or a new solver backend.
        if (n > 4096) error stop 'Linear-system diagnostic exceeds 4096 unknowns'
        if (size(rows) /= size(values)) error stop 'Diagnostic row size mismatch'
        if (size(cols) /= size(values)) error stop 'Diagnostic column size mismatch'
        if (size(bra, 2) /= n) error stop 'Diagnostic bra size mismatch'
        open(newunit=unit, file=path, status='new', action='write')
        write(unit, '(a)') '%%MatrixMarket matrix coordinate complex general'
        write(unit, '(3(i0,1x))') n, n, size(values)
        do j = 1, size(values)
            write(unit, '(2(i0,1x),2(es25.17e3,1x))') &
                rows(j), cols(j), real(values(j),dp), aimag(values(j))
        end do
        close(unit)
        call diagnostic_array(axisymmetric, species, 'flux_bra', transpose(bra))
    end subroutine diagnostic_system

    subroutine diagnostic_array(axisymmetric, species, label, values)
        logical, intent(in) :: axisymmetric
        integer, intent(in) :: species
        character(*), intent(in) :: label
        complex(dp), intent(in) :: values(:, :)
        character(:), allocatable :: path
        integer :: unit, row, col

        path = diagnostic_path(axisymmetric, species, label)
        if (len(path) == 0) return
        if (size(values,1) > 4096) &
            error stop 'Linear-system diagnostic exceeds 4096 unknowns'
        open(newunit=unit, file=path, status='new', action='write')
        write(unit, '(a)') '%%MatrixMarket matrix array complex general'
        write(unit, '(2(i0,1x))') size(values,1), size(values,2)
        do col = 1, size(values,2)
            do row = 1, size(values,1)
                write(unit, '(2(es25.17e3,1x))') &
                    real(values(row,col),dp), aimag(values(row,col))
            end do
        end do
        close(unit)
    end subroutine diagnostic_array
end module linear_system_diagnostic_mod
