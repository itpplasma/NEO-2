! The routine itself is extracted at configure time from the production source.
! Only its host data and one collective are supplied here. The mock collective
! checks this species' V projection, then gathers the other species' moments.
! This verifies multispecies arithmetic, NOT MPI communication or synchronization.
module native_integral_oracle
    use, intrinsic :: iso_fortran_env, only: real64
    use moment_test_support, only: check
    implicit none
    integer, parameter :: dp = real64
    integer :: lag, leg, ibeg, iend, npart, num_spec, ispec, n_2d_size
    integer, allocatable :: npl(:), ind_start(:)
    real(dp), allocatable :: pleg_bra(:, :, :), pleg_ket(:, :, :)
    real(dp), allocatable :: ailmm_aa(:, :, :, :, :), fact_pos_e(:), fact_neg_e(:)
    complex(dp), allocatable :: gathered(:, :, :, :)
    logical :: addboucol
    type :: mock_collective
        integer :: calls = 0
    contains
        procedure :: allgather_inplace
    end type
    type(mock_collective) :: mpro
contains
    subroutine allgather_inplace(this, local)
        class(mock_collective), intent(inout) :: this
        complex(dp), intent(inout) :: local(:, :, :, :)
        real(dp) :: scale
        scale = max(1.0_dp, maxval(abs(gathered(:, :, :, ispec))))
        call check(maxval(abs(local(:, :, :, ispec + 1) - gathered(:, :, :, ispec))) &
                   < 1.0e-13_dp*scale, 'native projection matches V before gather')
        local = gathered
        this%calls = this%calls + 1
    end subroutine
    include 'native_integral_part.inc'
end module
