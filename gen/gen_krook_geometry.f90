program gen_krook_geometry
    use, intrinsic :: iso_fortran_env, only: real64
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(*), operator(/), &
        operator(+), operator(-), operator(**)
    use fortsym_kernel_emit, only: kernel_emit_spec_t, emit_fortran_kernel_ir
    use fortsym_kernel_ir, only: kernel_ir_t, lower_kernel_ir
    use fortsym_string, only: chars, str, str_t
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: width, trapped_width, g_trapped, g_passing
    type(expr_t) :: geometry(4)
    type(kernel_ir_t) :: ir
    type(kernel_emit_spec_t) :: spec
    type(str_t) :: source
    character(len=2048) :: output, revision
    character(:), allocatable :: message, emitted
    integer :: unit, n, factorial
    logical :: ok

    call get_command_argument(1, output)
    call get_command_argument(2, revision)
    if (len_trim(output) == 0 .or. len_trim(revision) == 0) &
        error stop 'usage: gen_krook_geometry OUTPUT FORTSYM_REVISION'
    call arena%init()
    spec%generator = str('gen/gen_krook_geometry.f90')
    spec%generator_revision = str(trim(revision))
    spec%regenerate_command = &
        str('gen_krook_geometry COMMON/generated/krook_geometry_leaf.inc FORTSYM_REVISION')
    open (newunit=unit, file=trim(output), status='replace')
    width = sym(arena, 'width')
    trapped_width = sym(arena, 'trapped_width')
    g_trapped = sym(arena, 'g_trapped')
    g_passing = sym(arena, 'g_passing')
    ! Integral of g(lambda)*(lambda-lower)**(n-1)/(n-1)! d lambda,
    ! split at lower+trapped_width. The caller clips to a crossing band.
    factorial = 1
    do n = 1, 4
        factorial = factorial*n
        geometry(n) = (g_passing*width**n + &
            (g_trapped - g_passing)*trapped_width**n)/real(factorial, real64)
    end do
    call lower_kernel_ir(geometry, ir, ok, message)
    if (.not. ok) error stop message
    spec%name = str('krook_geometry_leaf')
    spec%args = [str('width'), str('trapped_width'), str('g_trapped'), &
        str('g_passing')]
    spec%outputs = [str('moment0'), str('moment1'), str('moment2'), str('moment3')]
    source = emit_fortran_kernel_ir(ir, spec, ok, message)
    if (.not. ok) error stop message
    emitted = chars(source)
    write (unit, '(a)') emitted(:verify(emitted, new_line('a'), back=.true.))
    close (unit)
end program gen_krook_geometry
