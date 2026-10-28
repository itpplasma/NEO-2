program gen_krook
    use, intrinsic :: iso_fortran_env, only: real64
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, sqrt, exp, operator(*), operator(/), &
        operator(+), operator(-), operator(**)
    use fortsym_kernel_emit, only: kernel_emit_spec_t, emit_fortran_kernel_ir
    use fortsym_kernel_ir, only: kernel_ir_t, lower_kernel_ir
    use fortsym_string, only: chars, str, str_t
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: nu, kappa, temperature, mass, roots(1)
    type(expr_t) :: speed, alpha, beta, prj, expansion, weights(2)
    type(expr_t) :: h_phi, bhat
    type(kernel_ir_t) :: ir
    type(kernel_emit_spec_t) :: spec
    type(str_t) :: source
    character(len=2048) :: output, revision
    character(:), allocatable :: message
    integer :: unit
    logical :: ok

    call get_command_argument(1, output)
    call get_command_argument(2, revision)
    if (len_trim(output) == 0 .or. len_trim(revision) == 0) &
        error stop 'usage: gen_krook OUTPUT FORTSYM_REVISION'
    call arena%init()
    nu = sym(arena, 'nu')
    kappa = sym(arena, 'kappa')
    temperature = sym(arena, 'temperature')
    mass = sym(arena, 'mass')
    ! Native rhs_mat_energ is minus one half of kappa times the pitch mass.
    roots(1) = 2.0_real64*nu/(kappa*sqrt(2.0_real64*temperature/mass))
    call lower_kernel_ir(roots, ir, ok, message)
    if (.not. ok) error stop message
    spec%name = str('krook_scale_leaf')
    spec%args = [str('nu'), str('kappa'), str('temperature'), str('mass')]
    spec%outputs = [str('scale')]
    spec%generator = str('gen/gen_krook.f90')
    spec%generator_revision = str(trim(revision))
    spec%regenerate_command = &
        str('gen_krook COMMON/generated/krook_leaf.inc FORTSYM_REVISION')
    source = emit_fortran_kernel_ir(ir, spec, ok, message)
    if (.not. ok) error stop message
    open (newunit=unit, file=trim(output), status='replace')
    write (unit, '(a)') chars(source)
    speed = sym(arena, 'speed')
    alpha = sym(arena, 'alpha')
    beta = sym(arena, 'beta')
    prj = sym(arena, 'prj')
    expansion = sym(arena, 'expansion')
    weights(2) = speed**alpha*exp(-(1.0_real64 + beta)*speed**2)*prj*expansion
    weights(1) = speed**3*weights(2)
    call lower_kernel_ir(weights, ir, ok, message)
    if (.not. ok) error stop message
    spec%name = str('krook_weight_leaf')
    spec%args = [str('speed'), str('alpha'), str('beta'), str('prj'), &
        str('expansion')]
    spec%outputs = [str('constant_weight'), str('inverse_cube_weight')]
    source = emit_fortran_kernel_ir(ir, spec, ok, message)
    if (.not. ok) error stop message
    write (unit, '(a)') chars(source)
    h_phi = sym(arena, 'h_phi')
    bhat = sym(arena, 'bhat')
    roots(1) = -kappa/(h_phi*bhat)
    call lower_kernel_ir(roots, ir, ok, message)
    if (.not. ok) error stop message
    spec%name = str('krook_pitch_leaf')
    spec%args = [str('kappa'), str('h_phi'), str('bhat')]
    spec%outputs = [str('coefficient')]
    source = emit_fortran_kernel_ir(ir, spec, ok, message)
    if (.not. ok) error stop message
    write (unit, '(a)', advance='no') chars(source)
    close (unit)
end program gen_krook
