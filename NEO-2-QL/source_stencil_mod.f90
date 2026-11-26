!> Choice of the three-point stencil in the QL source and flux quadrature.
!>
!> Each step of the forward sweep integrates over [istep-1, istep], each step of
!> the backward sweep over [istep, istep+1], with a quadratic through three
!> neighbouring points. The default choice alternates with the parity of
!> istep - ibeg. When iend - ibeg is odd (possible after a phi_divide re-split),
!> that choice needs a point beyond iend, so the stencil is shifted inward.
module source_stencil_mod
    implicit none
    private

    public :: forward_stencil_centred, backward_stencil_centred

contains

    !> .true.: points istep-1, istep, istep+1; .false.: istep-2, istep-1, istep.
    pure logical function forward_stencil_centred(istep, ibeg, iend)
        integer, intent(in) :: istep, ibeg, iend

        forward_stencil_centred = mod(istep - ibeg, 2) == 1 .and. istep < iend
    end function forward_stencil_centred

    !> .true.: points istep-1, istep, istep+1; .false.: istep, istep+1, istep+2.
    pure logical function backward_stencil_centred(istep, ibeg, iend)
        integer, intent(in) :: istep, ibeg, iend

        backward_stencil_centred = mod(istep - ibeg, 2) == 1 .or. istep + 2 > iend
    end function backward_stencil_centred

end module source_stencil_mod
