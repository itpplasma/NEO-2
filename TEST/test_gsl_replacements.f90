! Value-equivalence test for the FGSL/GSL replacements (issue #88).
! Links the system GSL as oracle; GSL is used here only, never in the
! shipped NEO-2 targets.
module gsl_oracle
    use, intrinsic :: iso_c_binding
    implicit none

    type, bind(c) :: gsl_function_t
        type(c_funptr) :: f
        type(c_ptr) :: params
    end type gsl_function_t

    interface
        function gsl_sf_bessel_kn(n, x) bind(c, name='gsl_sf_bessel_Kn')
            import :: c_int, c_double
            integer(c_int), value :: n
            real(c_double), value :: x
            real(c_double) :: gsl_sf_bessel_kn
        end function gsl_sf_bessel_kn

        function gsl_set_error_handler_off() &
            bind(c, name='gsl_set_error_handler_off')
            import :: c_ptr
            type(c_ptr) :: gsl_set_error_handler_off
        end function gsl_set_error_handler_off

        function gsl_integration_cquad_workspace_alloc(n) &
            bind(c, name='gsl_integration_cquad_workspace_alloc')
            import :: c_size_t, c_ptr
            integer(c_size_t), value :: n
            type(c_ptr) :: gsl_integration_cquad_workspace_alloc
        end function gsl_integration_cquad_workspace_alloc

        subroutine gsl_integration_cquad_workspace_free(w) &
            bind(c, name='gsl_integration_cquad_workspace_free')
            import :: c_ptr
            type(c_ptr), value :: w
        end subroutine gsl_integration_cquad_workspace_free

        function gsl_integration_cquad(f, a, b, epsabs, epsrel, ws, &
                                       result, abserr, nevals) &
            bind(c, name='gsl_integration_cquad')
            import :: c_int, c_double, c_size_t, c_ptr, gsl_function_t
            type(gsl_function_t) :: f
            real(c_double), value :: a, b, epsabs, epsrel
            type(c_ptr), value :: ws
            real(c_double) :: result, abserr
            integer(c_size_t) :: nevals
            integer(c_int) :: gsl_integration_cquad
        end function gsl_integration_cquad

        function gsl_integration_workspace_alloc(n) &
            bind(c, name='gsl_integration_workspace_alloc')
            import :: c_size_t, c_ptr
            integer(c_size_t), value :: n
            type(c_ptr) :: gsl_integration_workspace_alloc
        end function gsl_integration_workspace_alloc

        subroutine gsl_integration_workspace_free(w) &
            bind(c, name='gsl_integration_workspace_free')
            import :: c_ptr
            type(c_ptr), value :: w
        end subroutine gsl_integration_workspace_free

        function gsl_integration_qagiu(f, a, epsabs, epsrel, limit, ws, &
                                       result, abserr) &
            bind(c, name='gsl_integration_qagiu')
            import :: c_int, c_double, c_size_t, c_ptr, gsl_function_t
            type(gsl_function_t) :: f
            real(c_double), value :: a, epsabs, epsrel
            integer(c_size_t), value :: limit
            type(c_ptr), value :: ws
            real(c_double) :: result, abserr
            integer(c_int) :: gsl_integration_qagiu
        end function gsl_integration_qagiu

        function gsl_bspline_alloc(k, nbreak) bind(c, name='gsl_bspline_alloc')
            import :: c_size_t, c_ptr
            integer(c_size_t), value :: k, nbreak
            type(c_ptr) :: gsl_bspline_alloc
        end function gsl_bspline_alloc

        subroutine gsl_bspline_free(w) bind(c, name='gsl_bspline_free')
            import :: c_ptr
            type(c_ptr), value :: w
        end subroutine gsl_bspline_free

        function gsl_bspline_knots(breakpts, w) bind(c, name='gsl_bspline_knots')
            import :: c_int, c_ptr
            type(c_ptr), value :: breakpts, w
            integer(c_int) :: gsl_bspline_knots
        end function gsl_bspline_knots

        function gsl_bspline_eval(x, b, w) bind(c, name='gsl_bspline_eval')
            import :: c_int, c_double, c_ptr
            real(c_double), value :: x
            type(c_ptr), value :: b, w
            integer(c_int) :: gsl_bspline_eval
        end function gsl_bspline_eval

        function gsl_bspline_deriv_eval(x, nderiv, db, w) &
            bind(c, name='gsl_bspline_deriv_eval')
            import :: c_int, c_double, c_size_t, c_ptr
            real(c_double), value :: x
            integer(c_size_t), value :: nderiv
            type(c_ptr), value :: db, w
            integer(c_int) :: gsl_bspline_deriv_eval
        end function gsl_bspline_deriv_eval

        function gsl_vector_alloc(n) bind(c, name='gsl_vector_alloc')
            import :: c_size_t, c_ptr
            integer(c_size_t), value :: n
            type(c_ptr) :: gsl_vector_alloc
        end function gsl_vector_alloc

        subroutine gsl_vector_free(v) bind(c, name='gsl_vector_free')
            import :: c_ptr
            type(c_ptr), value :: v
        end subroutine gsl_vector_free

        subroutine gsl_vector_set(v, i, x) bind(c, name='gsl_vector_set')
            import :: c_size_t, c_double, c_ptr
            type(c_ptr), value :: v
            integer(c_size_t), value :: i
            real(c_double), value :: x
        end subroutine gsl_vector_set

        function gsl_vector_get(v, i) bind(c, name='gsl_vector_get')
            import :: c_size_t, c_double, c_ptr
            type(c_ptr), value :: v
            integer(c_size_t), value :: i
            real(c_double) :: gsl_vector_get
        end function gsl_vector_get

        function gsl_matrix_alloc(n1, n2) bind(c, name='gsl_matrix_alloc')
            import :: c_size_t, c_ptr
            integer(c_size_t), value :: n1, n2
            type(c_ptr) :: gsl_matrix_alloc
        end function gsl_matrix_alloc

        subroutine gsl_matrix_free(m) bind(c, name='gsl_matrix_free')
            import :: c_ptr
            type(c_ptr), value :: m
        end subroutine gsl_matrix_free

        function gsl_matrix_get(m, i, j) bind(c, name='gsl_matrix_get')
            import :: c_size_t, c_double, c_ptr
            type(c_ptr), value :: m
            integer(c_size_t), value :: i, j
            real(c_double) :: gsl_matrix_get
        end function gsl_matrix_get
    end interface

end module gsl_oracle


module test_integrands
    use, intrinsic :: iso_c_binding
    use nrtype, only : dp
    implicit none

contains

    function f_integrand(x)
        real(dp) :: f_integrand
        real(dp) :: x

        f_integrand = exp(-x*x)*cos(3.0d0*x)
    end function f_integrand

    function f_integrand_param(x, param1)
        real(dp) :: f_integrand_param
        real(dp) :: x, param1

        f_integrand_param = x**param1*exp(-x*x)
    end function f_integrand_param

    function c_integrand(x, params) bind(c) result(fx)
        real(c_double), value :: x
        type(c_ptr), value :: params
        real(c_double) :: fx

        real(dp), pointer :: p

        if (c_associated(params)) then
            call c_f_pointer(params, p)
            fx = x**p*exp(-x*x)
        else
            fx = exp(-x*x)*cos(3.0d0*x)
        end if
    end function c_integrand

end module test_integrands


program test_gsl_replacements
    use, intrinsic :: iso_c_binding
    use nrtype, only : dp
    use gsl_oracle
    use test_integrands
    use neo_bessel_k, only : bessel_kn
    use integration_routines_mod, only : fint1d_cquad, fint1d_qagiu
    use bspline_routines_mod, only : init_bspline, set_bspline_knots, &
                                     bspline_eval, bspline_deriv_eval

    implicit none

    integer :: nfail
    type(c_ptr) :: old_handler

    nfail = 0
    old_handler = gsl_set_error_handler_off()

    call test_bessel_kn
    call test_integration
    call test_bspline

    if (nfail .eq. 0) then
        print *, 'All tests passed!'
    else
        print *, 'FAIL: ', nfail, ' comparisons out of tolerance'
        stop 1
    end if

contains

    subroutine check(label, got, want, tol)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: got, want, tol

        if (abs(got - want) .gt. tol*max(1.0d0, abs(want))) then
            print *, 'FAIL ', label, ': got ', got, ' want ', want
            nfail = nfail + 1
        end if
    end subroutine check

    subroutine test_bessel_kn
        real(dp), dimension(7), parameter :: xs = &
            [0.05d0, 0.5d0, 1.0d0, 2.5d0, 10.0d0, 50.0d0, 400.0d0]
        integer :: n, i

        do n = 0, 2
            do i = 1, size(xs)
                call check('bessel_kn', bessel_kn(n, xs(i)), &
                           real(gsl_sf_bessel_kn(n, xs(i)), dp), 1.0d-13)
            end do
        end do
    end subroutine test_bessel_kn

    subroutine test_integration
        real(dp), parameter :: epsabs = 1.0d-10, epsrel = 1.0d-10
        real(dp), dimension(2) :: res
        real(dp) :: ra, rda
        real(dp), target :: param
        integer(c_size_t) :: neval
        integer(c_int) :: status
        type(c_ptr) :: ws
        type(gsl_function_t) :: f

        f%f = c_funloc(c_integrand)
        f%params = c_null_ptr

        ws = gsl_integration_cquad_workspace_alloc(100_c_size_t)
        status = gsl_integration_cquad(f, 0.0d0, 2.0d0, epsabs, epsrel, &
                                       ws, ra, rda, neval)
        call gsl_integration_cquad_workspace_free(ws)
        res = fint1d_cquad(f_integrand, 0.0d0, 2.0d0, epsabs, epsrel)
        call check('fint1d_cquad', res(1), ra, 1.0d-9)

        param = 1.5d0
        f%params = c_loc(param)
        ws = gsl_integration_cquad_workspace_alloc(100_c_size_t)
        status = gsl_integration_cquad(f, 0.0d0, 2.0d0, epsabs, epsrel, &
                                       ws, ra, rda, neval)
        call gsl_integration_cquad_workspace_free(ws)
        res = fint1d_cquad(f_integrand_param, param, 0.0d0, 2.0d0, &
                           epsabs, epsrel)
        call check('fint1d_cquad param', res(1), ra, 1.0d-9)

        f%params = c_null_ptr
        ws = gsl_integration_workspace_alloc(1000_c_size_t)
        status = gsl_integration_qagiu(f, 0.7d0, epsabs, epsrel, &
                                       1000_c_size_t, ws, ra, rda)
        call gsl_integration_workspace_free(ws)
        res = fint1d_qagiu(f_integrand, 0.7d0, epsabs, epsrel)
        call check('fint1d_qagiu', res(1), ra, 1.0d-9)
    end subroutine test_integration

    subroutine test_bspline
        integer, parameter :: order = 4, nbreak = 8
        integer, parameter :: ncbf = nbreak + order - 2
        real(dp), dimension(nbreak), target :: xknots
        real(dp), dimension(ncbf) :: b_new
        real(dp), dimension(11), parameter :: xs = &
            [0.0d0, 0.13d0, 0.5d0, 0.94d0, 1.0d0, 1.7d0, 2.3d0, 2.9d0, &
             3.3d0, 3.99d0, 4.0d0]
        real(dp), parameter :: dist = 1.4d0
        type(c_ptr) :: w, bv, brk, dbm
        integer(c_size_t) :: i, m
        integer(c_int) :: status
        integer :: k, nder
        real(dp) :: gam_all, x_del
        character(len=32) :: label

        gam_all = 0.0d0
        do k = 1, nbreak - 1
            gam_all = gam_all + dist**k
        end do
        xknots(1) = 0.0d0
        x_del = 4.0d0/gam_all
        do k = 1, nbreak - 1
            xknots(k + 1) = xknots(k) + x_del*dist**k
        end do
        xknots(nbreak) = 4.0d0

        call init_bspline(order, nbreak)
        call set_bspline_knots(xknots)

        w = gsl_bspline_alloc(int(order, c_size_t), int(nbreak, c_size_t))
        brk = gsl_vector_alloc(int(nbreak, c_size_t))
        do k = 1, nbreak
            call gsl_vector_set(brk, int(k - 1, c_size_t), xknots(k))
        end do
        status = gsl_bspline_knots(brk, w)
        bv = gsl_vector_alloc(int(ncbf, c_size_t))
        dbm = gsl_matrix_alloc(int(ncbf, c_size_t), int(order + 1, c_size_t))

        do k = 1, size(xs)
            status = gsl_bspline_eval(xs(k), bv, w)
            call bspline_eval(xs(k), b_new)
            do i = 0, int(ncbf - 1, c_size_t)
                write (label, '(A,F5.2,A,I2)') 'bspline x=', xs(k), ' m=', i
                call check(trim(label), b_new(i + 1), &
                           real(gsl_vector_get(bv, i), dp), 1.0d-12)
            end do

            do nder = 1, 2
                status = gsl_bspline_deriv_eval(xs(k), &
                                                int(order, c_size_t), dbm, w)
                call bspline_deriv_eval(xs(k), nder, b_new)
                do m = 0, int(ncbf - 1, c_size_t)
                    write (label, '(A,I1,A,F5.2,A,I2)') &
                        'bspline d', nder, ' x=', xs(k), ' m=', m
                    call check(trim(label), b_new(m + 1), &
                               real(gsl_matrix_get(dbm, m, &
                                                   int(nder, c_size_t)), dp), &
                               1.0d-10)
                end do
            end do
        end do

        call gsl_matrix_free(dbm)
        call gsl_vector_free(bv)
        call gsl_vector_free(brk)
        call gsl_bspline_free(w)
    end subroutine test_bspline

end program test_gsl_replacements
