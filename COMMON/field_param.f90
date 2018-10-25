module field_param
  use nrtype, only : dp
  ! No input file for a simple field without perturbation:
  logical :: prop=.false.
  integer       :: mmx, nmx, m_isl,n_isl
  real(kind=dp) :: br0,delta_m,delta_n,r0,aiota0,aiota_pr
  real(kind=dp) :: eps_isl
  real(kind=dp), dimension(:),   allocatable :: sinmt,cosmt,sinnp,cosnp
  real(kind=dp), dimension(:,:), allocatable :: expmn,sipsmn,copsmn
end module field_param
