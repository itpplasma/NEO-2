!> Units and Formats
module neo_units

  integer, parameter ::   r_u1   = 3
  integer, parameter ::   r_u2   = 4
  integer, parameter ::   r_us   = 5
  integer, parameter ::   r_u23  = 23
  integer, parameter ::   r_ua   = 21
  integer, parameter ::   w_us   = 6
  integer, parameter ::   w_u1   = 7
  integer, parameter ::   w_u2   = 8
  integer, parameter ::   w_u3   = 9
  integer, parameter ::   w_u4   = 10
  integer, parameter ::   w_u5   = 11
  integer, parameter ::   w_u6   = 12
  integer, parameter ::   w_u7   = 13
  integer, parameter ::   w_u8   = 14
  integer, parameter ::   w_u9   = 15
  integer, parameter ::   w_u10  = 16
  integer, parameter ::   w_u11  = 17
  integer, parameter ::   w_u12  = 18
  integer, parameter ::   w_u13  = 19
  integer, parameter ::   w_u14  = 20
  integer, parameter ::   w_u15  = 21
  integer, parameter ::   w_u16  = 22
  integer, parameter ::   w_u17  = 23

  integer            ::   w_u6_open
  character(20), parameter :: format220="(500d18.5)"

  character(30)                      :: base_file
  character(30)                      :: out_file
  character(30)                      :: chk_file
  character(30)                      :: epslog_file
  character(30)                      :: epscon_file
  character(30)                      :: epsdia_file
  character(30)                      :: curcon_file
  character(30)                      :: curint_file
  character(30)                      :: curdis_file
  character(30)                      :: epsadd_file
  character(30)                      :: cur_file
  character(30)                      :: pla_file
  character(30)                      :: sbc_file

end module neo_units
