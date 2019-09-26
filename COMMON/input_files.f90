module input_files
  character(len=1024) :: eqfile, cfile, gfile,pfile,convexfile,fluxdatapath
  integer :: iunit=1738
  integer :: ieqfile=1

  data eqfile  /'ASDEX/d3d-087506.03687.equ'/
  data cfile   /'DATA/ccoil.dat'/
!  data gfile   /'gfiles/shot115452/g115452.03525'/
!  data pfile   /'Conly/probe_g129_bfield.out'/
end module input_files
