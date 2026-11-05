program gen_phi_ma
    ! Pointwise leaves of the prescribed Phi^MA drive (DOC/phi_ma_response.md).
    use, intrinsic :: iso_fortran_env, only: real64
    use fortsym_expr, only: expr_t, sym, num, rat, pi_expr, i_expr, exp, sqrt, &
        operator(*), operator(/), operator(+), operator(-), operator(**)
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortsym_string, only: chars, str, str_t
    implicit none
    type(expr_t) :: temperature, mass, charge, light_speed, bref, mtovr, iota, sg
    type(expr_t) :: g, dn_ds, density, dt_ds, er, a1, a2, vt
    type(expr_t) :: n, m, bcov_theta, bcov_phi, phase, phi_ma, us
    type(expr_t) :: dlambda, rho, bhat, h_phi, jpar, potential
    type(expr_t) :: phi, theta, theta_beg, phi_beg, x, alpha, beta, prj, ii
    character(len=2048) :: qlfile, commonfile, flowfile, revision
    type(expr_t) :: eparb, volume, step_weight
    integer :: unit

    call get_command_argument(1, qlfile)
    call get_command_argument(2, commonfile)
    call get_command_argument(3, flowfile)
    call get_command_argument(4, revision)
    if (len_trim(qlfile) == 0 .or. len_trim(commonfile) == 0 &
        .or. len_trim(flowfile) == 0 .or. len_trim(revision) == 0) &
        error stop 'usage: gen_phi_ma QL_INC COMMON_INC FLOW_INC FORTSYM_REVISION'
    ii = i_expr()
    temperature = sym('temperature')
    mass = sym('mass')
    charge = sym('charge')
    light_speed = sym('light_speed')
    bref = sym('bref')
    mtovr = sym('mtovr')
    iota = sym('iota')
    sg = sym('sg')
    g = sym('g')
    dn_ds = sym('dn_ds')
    density = sym('density')
    dt_ds = sym('dt_ds')
    er = sym('er')
    a1 = sym('a1')
    n = sym('n')
    m = sym('m')
    bcov_theta = sym('bcov_theta')
    bcov_phi = sym('bcov_phi')
    phase = sym('phase')
    phi_ma = sym('phi_ma')
    us = sym('us')
    dlambda = sym('dlambda')
    rho = sym('rho')
    bhat = sym('bhat')
    h_phi = sym('h_phi')
    jpar = sym('jpar')
    potential = sym('potential')
    phi = sym('phi')
    theta = sym('theta')
    theta_beg = sym('theta_beg')
    phi_beg = sym('phi_beg')
    x = sym('x')
    alpha = sym('alpha')
    beta = sym('beta')
    prj = sym('prj')
    eparb = sym('eparb')
    volume = sym('volume')
    step_weight = sym('step_weight')

    open (newunit=unit, file=trim(qlfile), status='replace')
    vt = sqrt(num(2)*temperature/mass)
    ! Thermal speed, gyroradius at Bref, and E_r from ntv_mod's
    ! Om_tE = c*E_r/(iota*sg) with Om_tE = MtOvR*vt.
    call emit_real(unit, 'phi_ma_scales_leaf', &
        [vt, vt*mass*light_speed/(charge*bref), mtovr*vt*iota*sg/light_speed], &
        [str('temperature'), str('mass'), str('charge'), str('light_speed'), &
        str('bref'), str('mtovr'), str('iota'), str('sg')], &
        [str('vt'), str('rho'), str('er')])
    ! NEO-2's thermodynamic forces in the effective radius, Phi_0' = -E_r.
    a2 = g*dt_ds/temperature
    call emit_real(unit, 'thermodynamic_forces_leaf', &
        [g*dn_ds/density - rat(3, 2)*a2 - charge*er/temperature, a2], &
        [str('g'), str('dn_ds'), str('density'), str('dt_ds'), &
        str('temperature'), str('charge'), str('er')], [str('a1'), str('a2')])
    call emit_real(unit, 'field_line_theta_leaf', &
        [theta_beg + iota*(phi - phi_beg)], &
        [str('theta_beg'), str('iota'), str('phi'), str('phi_beg')], [str('theta')])
    call emit_complex(unit, 'phi_ma_phase_leaf', [exp(ii*(m*theta + n*phi))], &
        [str('m'), str('theta'), str('n'), str('phi')], [str('phase')], &
        [str('real(dp)'), str('real(dp)'), str('real(dp)'), str('real(dp)')])
    ! Radial ExB velocity of the misaligned harmonic in r (note Eq. (6)),
    ! times g = ds/dr because the stored solution is per unit s.
    call emit_complex(unit, 'phi_ma_drift_leaf', &
        [g*light_speed*ii*(n*bcov_theta - m*bcov_phi)*phase*phi_ma &
        /(sg*(iota*bcov_theta + bcov_phi))], &
        [str('g'), str('light_speed'), str('n'), str('bcov_theta'), str('m'), &
        str('bcov_phi'), str('phase'), str('phi_ma'), str('sg'), str('iota')], &
        [str('us')], [str('real(dp)'), str('real(dp)'), str('real(dp)'), &
        str('real(dp)'), str('real(dp)'), str('real(dp)'), str('complex(dp)'), &
        str('complex(dp)'), str('real(dp)'), str('real(dp)')])
    ! Pitch-integrated drive of one turning band (D7).
    call emit_complex(unit, 'phi_ma_band_leaf', &
        [num(-2)*us*dlambda/(sym('vt')*rho*bhat*h_phi)], &
        [str('us'), str('dlambda'), str('vt'), str('rho'), str('bhat'), &
        str('h_phi')], [str('q')], [str('complex(dp)'), str('real(dp)'), &
        str('real(dp)'), str('real(dp)'), str('real(dp)'), str('real(dp)')])
    ! J_hat = j_par/(charge*n*v_T,contract*A_1*Phi/E_r), v_T,contract = vt/sqrt(2).
    call emit_complex(unit, 'phi_ma_jhat_leaf', &
        [jpar*er/(charge*density*sqrt(temperature/mass)*a1*potential)], &
        [str('jpar'), str('er'), str('charge'), str('density'), &
        str('temperature'), str('mass'), str('a1'), str('potential')], [str('jhat')], &
        [str('complex(dp)'), str('real(dp)'), str('real(dp)'), str('real(dp)'), &
        str('real(dp)'), str('real(dp)'), str('real(dp)'), str('complex(dp)')])
    close (unit)

    open (newunit=unit, file=trim(commonfile), status='replace')
    ! Integrand of the x^-1 speed moment of the ExB source, with the basis
    ! measure combined before evaluation at x = 0.
    call emit_real(unit, 'exb_weight_leaf', &
        [pi_expr()**rat(-3, 2)*x**(num(3) + alpha) &
        *exp(-(num(1) + beta)*x**2)*prj], &
        [str('x'), str('alpha'), str('beta'), str('prj')], [str('value')])
    close (unit)

    open (newunit=unit, file=trim(flowfile), status='replace')
    ! A3 of compute_ParFlow from <E_par B>/<B^2>.
    call emit_real(unit, 'parallel_field_force_leaf', &
        [eparb*charge/temperature], &
        [str('eparb'), str('charge'), str('temperature')], [str('a3')])
    ! Surface-average weight of one RK step, local D3j per unit surface
    ! volume, and the conversion of the parallel flow v_par*B to a current.
    call emit_real(unit, 'local_parflow_scales_leaf', &
        [step_weight/(bhat*h_phi), volume*bhat*h_phi, charge/(bref*bhat)], &
        [str('step_weight'), str('volume'), str('bhat'), str('h_phi'), &
        str('charge'), str('bref')], [str('weight'), str('d3_scale'), str('j_scale')])
    close (unit)

contains

    subroutine emit_real(out, name, roots, args, outputs)
        integer, intent(in) :: out
        character(*), intent(in) :: name
        type(expr_t), intent(in) :: roots(:)
        type(str_t), intent(in) :: args(:), outputs(:)
        type(kernel_spec_t) :: spec

        call fill_spec(spec, name, args, outputs)
        call write_kernel(out, roots, spec)
    end subroutine emit_real

    subroutine emit_complex(out, name, roots, args, outputs, types)
        integer, intent(in) :: out
        character(*), intent(in) :: name
        type(expr_t), intent(in) :: roots(:)
        type(str_t), intent(in) :: args(:), outputs(:), types(:)
        type(kernel_spec_t) :: spec

        call fill_spec(spec, name, args, outputs)
        spec%scalar_type = str('complex(dp)')
        spec%arg_types = types
        call write_kernel(out, roots, spec)
    end subroutine emit_complex

    subroutine fill_spec(spec, name, args, outputs)
        type(kernel_spec_t), intent(inout) :: spec
        character(*), intent(in) :: name
        type(str_t), intent(in) :: args(:), outputs(:)

        spec%name = str(name)
        spec%args = args
        spec%outputs = outputs
        spec%elemental_procedure = .true.
        spec%generator = str('gen/gen_phi_ma.f90')
        spec%generator_revision = str(trim(revision))
        spec%regenerate_command = str('gen_phi_ma NEO-2-QL/generated/phi_ma_leaves.inc '// &
            'COMMON/generated/exb_weight_leaf.inc '// &
            'NEO-2-QL/generated/local_parflow_leaves.inc FORTSYM_REVISION')
    end subroutine fill_spec

    subroutine write_kernel(out, roots, spec)
        integer, intent(in) :: out
        type(expr_t), intent(in) :: roots(:)
        type(kernel_spec_t), intent(in) :: spec
        type(str_t) :: source
        character(:), allocatable :: message
        logical :: ok

        source = emit_kernel(roots, spec, ok, message=message)
        if (.not. ok) then
            print *, chars(spec%name), ': ', message
            error stop 'gen_phi_ma: emission failed'
        end if
        write (out, '(a)') chars(source)
    end subroutine write_kernel

end program gen_phi_ma
