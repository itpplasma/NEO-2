module helical_source_mod
    use nrtype, only: dp
    implicit none
    private

    public :: add_helical_source

contains

    subroutine add_helical_source(source, profile, moments, column, sigma_sign, &
            ibeg, iend, lag, npl, ind_start, fact_pos_b, fact_pos_e, &
            fact_neg_b, fact_neg_e)
        complex(dp), intent(inout) :: source(:, :)
        complex(dp), intent(in) :: profile(:, ibeg:)
        real(dp), intent(in) :: moments(0:)
        integer, intent(in) :: column, sigma_sign, ibeg, iend, lag
        integer, intent(in) :: npl(ibeg:iend), ind_start(ibeg:iend)
        real(dp), intent(in) :: fact_pos_b(ibeg:iend), fact_pos_e(ibeg:iend)
        real(dp), intent(in) :: fact_neg_b(ibeg:iend), fact_neg_e(ibeg:iend)

        call add_forward_source(source, profile, moments, column, ibeg, iend, &
            lag, npl, ind_start, fact_pos_b, fact_pos_e)
        call add_backward_source(source, profile, moments, column, sigma_sign, &
            ibeg, iend, lag, npl, ind_start, fact_neg_b, fact_neg_e)
    end subroutine add_helical_source

    subroutine add_forward_source(source, profile, moments, column, ibeg, iend, &
            lag, npl, ind_start, fact_b, fact_e)
        complex(dp), intent(inout) :: source(:, :)
        complex(dp), intent(in) :: profile(:, ibeg:)
        real(dp), intent(in) :: moments(0:), fact_b(ibeg:iend), fact_e(ibeg:iend)
        integer, intent(in) :: column, ibeg, iend, lag
        integer, intent(in) :: npl(ibeg:iend), ind_start(ibeg:iend)
        integer :: i, m, k, k_prev, np, np_prev, np_next

        do i = ibeg, iend
            np = npl(i)
            do m = 0, lag
                k = ind_start(i) + 2*(np + 1)*m
                if (i <= ibeg) cycle
                np_prev = npl(i - 1)
                k_prev = ind_start(i - 1) + 2*(np_prev + 1)*m
                if (mod(i - ibeg, 2) == 1) then
                    np_next = npl(i + 1)
                    source(k + 1:k + np + 1, column) = source(k + 1:k + np + 1, column) &
                        + moments(m)/1.5d0*profile(1:np + 1, i)*fact_e(i)
                    source(k + 1:k + np_prev + 1, column) = source(k + 1:k + np_prev + 1, column) &
                        + moments(m)/2.4d0*profile(1:np_prev + 1, i - 1)*fact_b(i - 1)
                    source(k + 1:k + np_next + 1, column) = source(k + 1:k + np_next + 1, column) &
                        - moments(m)/12d0*profile(1:np_next + 1, i + 1)*fact_e(i + 1)
                else
                    np_next = npl(i - 2)
                    source(k + 1:k + np + 1, column) = source(k + 1:k + np + 1, column) &
                        + moments(m)/2.4d0*profile(1:np + 1, i)*fact_e(i)
                    call add_forward_neighbor(source, profile, moments(m), column, k, k_prev, &
                        np, np_prev, i - 1, ibeg, 1.5d0, fact_b(i - 1))
                    call add_forward_neighbor(source, profile, -moments(m), column, k, k_prev, &
                        np, np_next, i - 2, ibeg, 12d0, fact_b(i - 2))
                end if
            end do
        end do
    end subroutine add_forward_source

    subroutine add_forward_neighbor(source, profile, moment, column, k, k_prev, &
            np, np_neighbor, i_neighbor, ibeg, divisor, factor)
        complex(dp), intent(inout) :: source(:, :)
        complex(dp), intent(in) :: profile(:, ibeg:)
        real(dp), intent(in) :: moment, divisor, factor
        integer, intent(in) :: column, k, k_prev, np, np_neighbor, i_neighbor, ibeg

        if (np_neighbor <= np) then
            source(k + 1:k + np_neighbor + 1, column) = &
                source(k + 1:k + np_neighbor + 1, column) &
                + moment/divisor*profile(1:np_neighbor + 1, i_neighbor)*factor
        else
            source(k + 1:k + np + 1, column) = source(k + 1:k + np + 1, column) &
                + moment/divisor*profile(1:np + 1, i_neighbor)*factor
            source(k_prev + np_neighbor + 2, column) = source(k_prev + np_neighbor + 2, column) &
                + moment/divisor*profile(np_neighbor + 1, i_neighbor)*factor
        end if
    end subroutine add_forward_neighbor

    subroutine add_backward_source(source, profile, moments, column, sigma_sign, &
            ibeg, iend, lag, npl, ind_start, fact_b, fact_e)
        complex(dp), intent(inout) :: source(:, :)
        complex(dp), intent(in) :: profile(:, ibeg:)
        real(dp), intent(in) :: moments(0:), fact_b(ibeg:iend), fact_e(ibeg:iend)
        integer, intent(in) :: column, sigma_sign, ibeg, iend, lag
        integer, intent(in) :: npl(ibeg:iend), ind_start(ibeg:iend)
        integer :: i, m, k, k_prev, np, np_prev, np_next
        real(dp) :: signed_moment

        do i = ibeg, iend
            np = npl(i)
            do m = 0, lag
                k = ind_start(i) + 2*(np + 1)*m
                if (i >= iend) cycle
                signed_moment = real(sigma_sign, dp)*moments(m)
                np_prev = npl(i + 1)
                k_prev = ind_start(i + 1) + 2*(np_prev + 1)*m
                if (mod(i - ibeg, 2) == 1) then
                    np_next = npl(i - 1)
                    source(k + np + 2:k + 2*np + 2, column) = &
                        source(k + np + 2:k + 2*np + 2, column) &
                        + signed_moment/1.5d0*profile(np + 1:1:-1, i)*fact_e(i)
                    source(k + 2*np + 2 - np_prev:k + 2*np + 2, column) = &
                        source(k + 2*np + 2 - np_prev:k + 2*np + 2, column) &
                        + signed_moment/2.4d0*profile(np_prev + 1:1:-1, i + 1)*fact_b(i + 1)
                    source(k + 2*np + 2 - np_next:k + 2*np + 2, column) = &
                        source(k + 2*np + 2 - np_next:k + 2*np + 2, column) &
                        - signed_moment/12d0*profile(np_next + 1:1:-1, i - 1)*fact_e(i - 1)
                else
                    np_next = npl(i + 2)
                    source(k + np + 2:k + 2*np + 2, column) = &
                        source(k + np + 2:k + 2*np + 2, column) &
                        + signed_moment/2.4d0*profile(np + 1:1:-1, i)*fact_e(i)
                    call add_backward_neighbor(source, profile, signed_moment, column, k, k_prev, &
                        np, np_prev, i + 1, ibeg, 1.5d0, fact_b(i + 1))
                    call add_backward_neighbor(source, profile, -signed_moment, column, k, k_prev, &
                        np, np_next, i + 2, ibeg, 12d0, fact_b(i + 2))
                end if
            end do
        end do
    end subroutine add_backward_source

    subroutine add_backward_neighbor(source, profile, moment, column, k, k_prev, &
            np, np_neighbor, i_neighbor, ibeg, divisor, factor)
        complex(dp), intent(inout) :: source(:, :)
        complex(dp), intent(in) :: profile(:, ibeg:)
        real(dp), intent(in) :: moment, divisor, factor
        integer, intent(in) :: column, k, k_prev, np, np_neighbor, i_neighbor, ibeg

        if (np_neighbor <= np) then
            source(k + 2*np + 2 - np_neighbor:k + 2*np + 2, column) = &
                source(k + 2*np + 2 - np_neighbor:k + 2*np + 2, column) &
                + moment/divisor*profile(np_neighbor + 1:1:-1, i_neighbor)*factor
        else
            source(k + np + 2:k + 2*np + 2, column) = &
                source(k + np + 2:k + 2*np + 2, column) &
                + moment/divisor*profile(np + 1:1:-1, i_neighbor)*factor
            source(k_prev + np_neighbor + 1, column) = source(k_prev + np_neighbor + 1, column) &
                + moment/divisor*profile(np_neighbor + 1, i_neighbor)*factor
        end if
    end subroutine add_backward_neighbor
end module helical_source_mod
