program gen_covariant_curl
    use iso_fortran_env, only: output_unit
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t, VERDICT_TRUE
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_kernel_ir, only: kernel_ir_t, lower_kernel_ir
    use fortsym_kernel_emit, only: kernel_emit_spec_t, emit_fortran_kernel_ir
    use fortsym_string, only: str_t, str, chars
    implicit none
    character(*), parameter :: revision = &
        '7e8ef7abdd90c49bb5f75e2c119e7d181656bcc6'
    character(8), parameter :: names(12) = [character(8) :: &
        'b', 'c', 'g', 'i', 'bs', 'bt', 'bp', 'ct', 'cp', 'is', 'gs', 'jac']
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(expr_t) :: v(12), tangent(12), weight(3), h(3), partial(3,3)
    type(expr_t) :: curl(3), expanded(3), jacobian(3,12), jvp(3), vjp(12), pairing
    type(engine_result_t) :: checked
    type(str_t) :: arguments(12), directions(12), weights(3), outputs(12)
    character(2048) :: path
    integer :: j, k, unit

    call get_command_argument(1, path)
    if (len_trim(path) == 0) error stop 'usage: gen_covariant_curl OUTPUT'
    call arena%init()
    engine = make_native_engine(arena)
    do j = 1, 12
        v(j) = sym(arena, trim(names(j)))
        arguments(j) = str(trim(names(j)))
        directions(j) = str('d'//trim(names(j)))
        tangent(j) = sym(arena, chars(directions(j)))
    end do
    do j = 1, 3
        weights(j) = str('w'//integer_text(j))
        weight(j) = sym(arena, chars(weights(j)))
    end do

    ! One mathematical source: returned CGS covariant unit field (s,phi,theta).
    ! Domain B>0, final signed CGS jac /= 0. G and I are flux functions.
    h = [100*v(2)/v(1), 100*v(3)/v(1), 100*v(4)/v(1)]
    do j = 1, 3
        partial(j,1) = diff(h(j),v(1))*v(5) &
            + diff(h(j),v(3))*v(11) + diff(h(j),v(4))*v(10)
        partial(j,2) = diff(h(j),v(1))*v(7) + diff(h(j),v(2))*v(9)
        partial(j,3) = diff(h(j),v(1))*v(6) + diff(h(j),v(2))*v(8)
    end do
    ! C_s is deliberately absent: the Levi-Civita contraction never uses it.
    curl = [(partial(3,2)-partial(2,3))/v(12), &
            (partial(1,3)-partial(3,1))/v(12), &
            (partial(2,1)-partial(1,2))/v(12)]
    expanded = [100*(v(3)*v(6)-v(4)*v(7))/(v(1)**2*v(12)), &
        100*(v(1)*(v(8)-v(10))-v(2)*v(6)+v(4)*v(5))/(v(1)**2*v(12)), &
        100*(v(1)*(v(11)-v(9))-v(3)*v(5)+v(2)*v(7))/(v(1)**2*v(12))]
    do j = 1, 3
        checked = engine%zero_test(curl(j)-expanded(j))
        if (checked%verdict /= VERDICT_TRUE) error stop 'curl identity unproved'
        curl(j) = simplified(curl(j))
        do k = 1, 12
            jacobian(j,k) = simplified(diff(curl(j), v(k)))
        end do
        jvp(j) = jacobian(j,1)*tangent(1)
        do k = 2, 12
            jvp(j) = jvp(j) + jacobian(j,k)*tangent(k)
        end do
        jvp(j) = simplified(jvp(j))
    end do
    do k = 1, 12
        vjp(k) = jacobian(1,k)*weight(1)
        do j = 2, 3
            vjp(k) = vjp(k) + jacobian(j,k)*weight(j)
        end do
        vjp(k) = simplified(vjp(k))
    end do

    pairing = weight(1)*jvp(1) + weight(2)*jvp(2) + weight(3)*jvp(3)
    do k = 1, 12
        pairing = pairing - tangent(k)*vjp(k)
    end do
    checked = engine%zero_test(pairing)
    if (checked%verdict /= VERDICT_TRUE) error stop 'JVP/VJP pairing unproved'

    open(newunit=unit, file=trim(path), status='replace')
    write(unit,'(a)') '! Generated; mathematical source: gen/gen_covariant_curl.f90'
    write(unit,'(a)') '! Domain: finite inputs, B>0, signed final CGS jacobian /=0.'
    write(unit,'(a)') 'module neo_covariant_curl_generated'
    write(unit,'(a)') '    implicit none'
    write(unit,'(a)') '    private'
    write(unit,'(a)') '    public :: covariant_curl, covariant_curl_jvp, covariant_curl_vjp'
    write(unit,'(a)') 'contains'
    do j = 1, 3
        outputs(j) = str('curl'//integer_text(j))
    end do
    call emit('covariant_curl', arguments, outputs(:3), curl)
    call emit('covariant_curl_jvp', [arguments,directions], outputs(:3), jvp)
    do j = 1, 12
        outputs(j) = str('bar_'//trim(names(j)))
    end do
    call emit('covariant_curl_vjp', [arguments,weights], outputs, vjp)
    write(unit,'(a)') 'end module neo_covariant_curl_generated'
    close(unit)
    write(output_unit,'(a)') 'PASS symbolic curl identities and JVP/VJP pairing (B>0, jac/=0)'
contains
    function simplified(expression) result(value)
        type(expr_t), intent(in) :: expression
        type(expr_t) :: value
        type(engine_result_t) :: result
        result = engine%simplify(expression)
        if (.not. result%ok) error stop 'native simplification failed'
        value = result%value
    end function simplified

    function integer_text(number) result(text)
        integer, intent(in) :: number
        character(:), allocatable :: text
        character(10) :: buffer
        write(buffer,'(i0)') number
        text = trim(buffer)
    end function integer_text

    subroutine emit(name, args, results, roots)
        character(*), intent(in) :: name
        type(str_t), intent(in) :: args(:), results(:)
        type(expr_t), intent(in) :: roots(:)
        type(kernel_ir_t) :: ir
        type(kernel_emit_spec_t) :: spec
        type(str_t) :: source
        character(:), allocatable :: message
        logical :: ok
        call lower_kernel_ir(roots, ir, ok, message)
        if (.not. ok) error stop 'curl IR lowering failed'
        spec%name = str(name)
        spec%args = args
        spec%outputs = results
        spec%temp_prefix = str('t')
        spec%pure_procedure = .true.
        spec%generator = str('gen_covariant_curl')
        spec%generator_revision = str(revision)
        spec%regenerate_command = str('see gen/README.md')
        source = emit_fortran_kernel_ir(ir, spec, ok, message)
        if (.not. ok) error stop 'curl Fortran emission failed'
        write(unit,'(a)') chars(source)
    end subroutine emit
end program gen_covariant_curl
