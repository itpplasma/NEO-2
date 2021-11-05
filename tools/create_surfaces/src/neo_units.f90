MODULE neo_units
! Units and Formats

  INTEGER, PARAMETER ::   r_u1   = 3
  INTEGER, PARAMETER ::   r_u2   = 4
  INTEGER, PARAMETER ::   r_us   = 5
  INTEGER, PARAMETER ::   r_u23  = 23
  INTEGER, PARAMETER ::   r_ua   = 21
  INTEGER, PARAMETER ::   w_us   = 6
  INTEGER, PARAMETER ::   w_u1   = 7
  INTEGER, PARAMETER ::   w_u2   = 8
  INTEGER, PARAMETER ::   w_u3   = 9
  INTEGER, PARAMETER ::   w_u4   = 10
  INTEGER, PARAMETER ::   w_u5   = 11
  INTEGER, PARAMETER ::   w_u6   = 12
  INTEGER, PARAMETER ::   w_u7   = 13
  INTEGER, PARAMETER ::   w_u8   = 14
  INTEGER, PARAMETER ::   w_u9   = 15
  INTEGER, PARAMETER ::   w_u10  = 16
  INTEGER, PARAMETER ::   w_u11  = 17
  INTEGER, PARAMETER ::   w_u12  = 18
  INTEGER, PARAMETER ::   w_u13  = 19
  INTEGER, PARAMETER ::   w_u14  = 20
  INTEGER, PARAMETER ::   w_u15  = 21
  INTEGER, PARAMETER ::   w_u16  = 22
  INTEGER, PARAMETER ::   w_u17  = 23

  INTEGER            ::   w_u6_open
  CHARACTER(20),PARAMETER :: format220="(500d18.5)"

  CHARACTER(30)                      :: base_file
  CHARACTER(30)                      :: out_file
  CHARACTER(30)                      :: chk_file
  CHARACTER(30)                      :: epslog_file
  CHARACTER(30)                      :: epscon_file
  CHARACTER(30)                      :: epsdia_file
  CHARACTER(30)                      :: curcon_file
  CHARACTER(30)                      :: curint_file
  CHARACTER(30)                      :: curdis_file
  CHARACTER(30)                      :: epsadd_file
  CHARACTER(30)                      :: cur_file
  CHARACTER(30)                      :: pla_file
  CHARACTER(30)                      :: sbc_file

END MODULE neo_units
