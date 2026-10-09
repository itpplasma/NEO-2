! A regular analytic Fourier field on two surfaces: independent interpolation
! oracle, including a sine mode. Before the fix the second allocation aborts.
program test_magfie_multisurface
  use nrtype, only: dp
  use neo_control, only: in_file,inp_swi,lab_swi,write_progress,fluxs_interp, &
    max_m_mode,max_n_mode,INP_SWI_TOK,theta_n,phi_n
  use neo_sub_mod, only: neo_read,neo_prep,neo_init_spline
  use neo_spline_data, only: lsw_linear_boozer
  use neo_magfie, only: magfie_sarray,magfie_spline,magfie_result,neo_magfie_a
  use ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: u,i,j
  real(dp) :: s,t,x(3),b,jac,bd(3),hc(3),ht(3),curl(3),truth,pi
  in_file='analytic_multisurface.bc'
  open(newunit=u,file=trim(in_file),status='replace')
  write(u,'(a)') 'CC analytic', 'CC two surfaces', 'CC SI', 'CC axisymmetric', 'm n ns nfp flux a R'
  write(u,*) 1,0,40,1,2._dp,1._dp,4._dp
  do i=1,40
    s=(i-.5_dp)/40
    write(u,'(a)') 's iota Jpol Itor pprime sqrtg', 'SI'
    write(u,'(6es25.16)') s,.7_dp+.1_dp*s,-2e7_dp,-1e6_dp,0._dp,1._dp
    write(u,'(a)') 'm n Rcos Rsin Zcos Zsin vcos vsin Bcos Bsin'
    write(u,'(2i5,8es25.16)') 0,0,4._dp,0._dp,0._dp,0._dp,0._dp,0._dp,4+.25_dp*s,0._dp
    write(u,'(2i5,8es25.16)') 1,0,sqrt(s),0._dp,0._dp,sqrt(s),0._dp,0._dp,.2_dp*sqrt(s),-.05_dp*sqrt(s)
  end do
  close(u)
  inp_swi=INP_SWI_TOK; lab_swi=10; write_progress=0; fluxs_interp=1
  max_m_mode=99; max_n_mode=99; theta_n=129; phi_n=9
  lsw_linear_boozer=.false.
  call neo_read()
  call neo_prep()
  call neo_init_spline()
  allocate(magfie_sarray(2))
  magfie_sarray=[.2_dp,.8_dp]
  magfie_spline=1; magfie_result=0
  pi=4*atan(1._dp)
  do i=1,2
    do j=0,31
      s=magfie_sarray(i); t=2*pi*(j+.37_dp)/32
      x=[s,.123_dp,t]
      call neo_magfie_a(x,b,jac,bd,hc,ht,curl)
      truth=4+.25_dp*s+sqrt(s)*(.2_dp*cos(t)-.05_dp*sin(t))
      if(abs(b-truth)>5e-9_dp) error stop 'wrong analytic field'
      if (.not.all(ieee_is_finite([b,jac,bd,hc,ht,curl]))) error stop 'nonfinite geometry'
      if (abs(ht(3)/ht(2)-(.7_dp+.1_dp*s))>1e-12_dp) error stop 'wrong field pitch'
    end do
  end do
  open(newunit=u,file=trim(in_file),status='old')
  close(u,status='delete')
  print *, 'Two-surface analytic magfie test passed'
end program
