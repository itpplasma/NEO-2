! Exact factorization of the retained QL integral collision action.
! This module does not change the production solver or its collision truncation.
module collision_moment_map
    use, intrinsic :: iso_fortran_env, only: real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    integer, parameter :: dp = real64
    integer, parameter, public :: moment_map_ok = 0, moment_map_invalid = 1

    type, public :: collision_moment_map_t
        private
        integer :: nf = 0, nr = 0, nv = 0, nl = 0, np = 0, width = 0
        logical :: add_boundary = .false.
        integer, allocatable :: passing(:), offset(:)
        real(dp), allocatable :: bra(:, :, :), ket(:, :, :), forward(:), backward(:)
    contains
        procedure :: field_size
        procedure :: moment_size
        procedure :: project => project_moments
        procedure :: lift => lift_moments
    end type
    public :: init_collision_moment_map, mix_collision_moments

contains
    integer function field_size(this)
        class(collision_moment_map_t), intent(in) :: this
        field_size = this%nf
    end function

    integer function moment_size(this)
        class(collision_moment_map_t), intent(in) :: this
        moment_size = this%nr
    end function

    ! Arrays use one-based storage: bra(l+1,p,phi), with physical l starting at 0.
    ! The distribution contains co/counter blocks for each energy basis and phi.
    subroutine init_collision_moment_map(map, nvelocity, passing, bra, ket, &
                                        forward, backward, add_boundary, info)
        type(collision_moment_map_t), intent(out) :: map
        integer, intent(in) :: nvelocity, passing(:)
        real(dp), intent(in) :: bra(:, :, :), ket(:, :, :), forward(:), backward(:)
        logical, intent(in) :: add_boundary
        integer, intent(out) :: info
        integer :: i, p, l, stat
        integer(int64) :: nf64, nr64

        info = moment_map_invalid
        p = size(passing)
        l = size(bra, 1)
        if (nvelocity < 1 .or. p < 1 .or. l < 1) return
        if (any(passing < 0) .or. any(passing == huge(0))) return
        if (size(bra, 2) < maxval(passing) + 1 .or. size(bra, 3) /= p) return
        if (any(shape(ket) /= shape(bra))) return
        if (size(forward) /= p .or. size(backward) /= p) return
        if (.not. all(ieee_is_finite(bra)) .or. .not. all(ieee_is_finite(ket))) return
        if (.not. all(ieee_is_finite(forward)) .or. .not. all(ieee_is_finite(backward))) return
        nf64 = 2_int64*int(nvelocity, int64)*sum(int(passing, int64) + 1_int64)
        nr64 = int(nvelocity, int64)*int(l, int64)*int(p, int64)
        if (nf64 > huge(0) .or. nr64 > huge(0)) return
        allocate (map%passing(p), map%offset(p), map%bra(l, size(bra, 2), p), &
                  map%ket(l, size(bra, 2), p), map%forward(p), map%backward(p), stat=stat)
        if (stat /= 0) return
        map%nf = int(nf64)
        map%nr = int(nr64)
        map%nv = nvelocity
        map%nl = l
        map%np = p
        map%width = 2*(maxval(passing) + 1)
        map%add_boundary = add_boundary
        map%passing = passing
        map%bra = bra
        map%ket = ket
        map%forward = forward
        map%backward = backward
        map%offset(1) = 0
        do i = 2, p
            map%offset(i) = map%offset(i - 1) + 2*(passing(i - 1) + 1)*nvelocity
        end do
        info = moment_map_ok
    end subroutine

    ! V: parity-aware Legendre projection. m is ordered (energy, Legendre, phi).
    subroutine project_moments(this, f, m, info)
        class(collision_moment_map_t), intent(in) :: this
        complex(dp), intent(in) :: f(:)
        complex(dp), intent(out) :: m(:)
        integer, intent(out) :: info
        integer :: i, v, l, k, n, j

        info = moment_map_invalid
        m = (0.0_dp, 0.0_dp)
        if (this%nf == 0 .or. size(f) /= this%nf .or. size(m) /= this%nr) return
        if (.not. finite_complex(f)) return
        do i = 1, this%np
            n = this%passing(i) + 1
            do v = 1, this%nv
                k = this%offset(i) + 2*n*(v - 1)
                do l = 1, this%nl
                    j = v + this%nv*((l - 1) + this%nl*(i - 1))
                    m(j) = sum(this%bra(l, 1:n, i)*f(k + 1:k + n)) &
                         + (-1.0_dp)**(l - 1)*sum(this%bra(l, 1:n, i)*f(k + 2*n:k + n + 1:-1))
                end do
            end do
        end do
        if (finite_complex(m)) info = moment_map_ok
    end subroutine

    ! C: exact energy/species matrix already retained by the collision model.
    ! coeff(mout,min,l,sout,sin), moments(energy-fastest/Legendre/phi,species).
    ! No inversion of C, rank assumption, or extra moment truncation is made.
    subroutine mix_collision_moments(coeff, moments, mixed, info)
        real(dp), intent(in) :: coeff(:, :, :, :, :)
        complex(dp), intent(in) :: moments(:, :)
        complex(dp), intent(out) :: mixed(:, :)
        integer, intent(out) :: info
        integer :: nv, nl, ns, np, i, l, a, b, j

        info = moment_map_invalid
        mixed = (0.0_dp, 0.0_dp)
        nv = size(coeff, 1); nl = size(coeff, 3); ns = size(coeff, 4)
        if (nv < 1 .or. nl < 1 .or. ns < 1) return
        if (size(coeff, 2) /= nv .or. size(coeff, 5) /= ns) return
        if (size(moments, 2) /= ns .or. any(shape(mixed) /= shape(moments))) return
        if (size(moments, 1) == 0) return
        if (int(nv, int64)*int(nl, int64) > size(moments, 1)) return
        if (mod(size(moments, 1), nv*nl) /= 0) return
        if (.not. all(ieee_is_finite(coeff))) return
        do a = 1, ns
            if (.not. finite_complex(moments(:, a))) return
        end do
        np = size(moments, 1)/(nv*nl)
        do i = 1, np
            do l = 1, nl
                j = nv*((l - 1) + nl*(i - 1))
                do a = 1, ns
                    do b = 1, ns
                        mixed(j + 1:j + nv, a) = mixed(j + 1:j + nv, a) &
                            + matmul(coeff(:, :, l, a, b), moments(j + 1:j + nv, b))
                    end do
                end do
            end do
        end do
        do a = 1, ns
            if (.not. finite_complex(mixed(:, a))) return
        end do
        info = moment_map_ok
    end subroutine

    ! U: Legendre reconstruction AND the original spatial deposition, including
    ! the 1/2 factor, both directions, ragged pitch support, and bounce terms.
    subroutine lift_moments(this, m, f, info)
        class(collision_moment_map_t), intent(in) :: this
        complex(dp), intent(in) :: m(:)
        complex(dp), intent(out) :: f(:)
        integer, intent(out) :: info
        complex(dp), allocatable :: tmp(:, :, :)
        integer :: i, v, l, j, k, n, nsame, knext, stat

        info = moment_map_invalid
        f = (0.0_dp, 0.0_dp)
        if (this%nf == 0 .or. size(f) /= this%nf .or. size(m) /= this%nr) return
        if (.not. finite_complex(m)) return
        allocate (tmp(this%nv, this%width, this%np), stat=stat)
        if (stat /= 0) return
        tmp = (0.0_dp, 0.0_dp)
        do i = 1, this%np
            n = this%passing(i) + 1
            do v = 1, this%nv
                do l = 1, this%nl
                    j = v + this%nv*((l - 1) + this%nl*(i - 1))
                    tmp(v, 1:n, i) = tmp(v, 1:n, i) + m(j)*this%ket(l, 1:n, i)
                    tmp(v, 2*n:n + 1:-1, i) = tmp(v, 2*n:n + 1:-1, i) &
                        + (-1.0_dp)**(l - 1)*m(j)*this%ket(l, 1:n, i)
                end do
            end do
        end do
        tmp = 0.5_dp*tmp
        do i = 2, this%np
            n = this%passing(i) + 1
            nsame = min(this%passing(i), this%passing(i - 1)) + 1
            do v = 1, this%nv
                k = this%offset(i) + 2*n*(v - 1)
                f(k + 1:k + n) = f(k + 1:k + n) + tmp(v, 1:n, i)*this%forward(i)
                f(k + 1:k + nsame) = f(k + 1:k + nsame) + tmp(v, 1:nsame, i - 1)*this%forward(i)
                if (this%passing(i) > this%passing(i - 1) .and. this%add_boundary) &
                    f(k + n) = f(k + n) + tmp(v, n + 1, i)*this%forward(i)
            end do
        end do
        do i = 1, this%np - 1
            n = this%passing(i) + 1
            nsame = min(this%passing(i), this%passing(i + 1)) + 1
            knext = 2*(this%passing(i + 1) + 1)
            do v = 1, this%nv
                k = this%offset(i) + 2*n*v
                f(k - n + 1:k) = f(k - n + 1:k) + tmp(v, n + 1:2*n, i)*this%backward(i)
                f(k - nsame + 1:k) = f(k - nsame + 1:k) &
                    + tmp(v, knext - nsame + 1:knext, i + 1)*this%backward(i)
                if (this%passing(i) > this%passing(i + 1) .and. this%add_boundary) &
                    f(k - n + 1) = f(k - n + 1) + tmp(v, n, i)*this%backward(i)
            end do
        end do
        if (finite_complex(f)) info = moment_map_ok
    end subroutine

    logical function finite_complex(x)
        complex(dp), intent(in) :: x(:)
        finite_complex = all(ieee_is_finite(real(x, dp))) .and. all(ieee_is_finite(aimag(x)))
    end function
end module collision_moment_map
