character(len=*), parameter :: CMake_Compiler = 'GNU'
character(len=*), parameter :: CMake_Compiler_Version = '4.7'
character(len=*), parameter :: CMake_Build_Type = 'Debug'
character(len=*), parameter :: CMake_Flags = '-cpp -pthread'
character(len=*), parameter :: CMake_Flags_Release = '-O3'
character(len=*), parameter :: CMake_Flags_Debug = '-g'
character(len=*), parameter :: CMake_System = 'Linux-3.14.25'
character(len=*), parameter :: CMake_SuiteSparse_Dir = '/proj/plasma/Libs/SuiteSparse_4.2.1/'
character(len=*), parameter :: CMake_Blas_Lib = '-llapack -lblas'
character(len=*), parameter :: Neo2_Version = '219060fe8b7907067e3ed44162059439d431aee2'
character(len=*), parameter :: Neo2_Version_Additional = "WARNING, &
&THERE ARE UNCOMMITTED CHANGES. Run may not be reproduceable: &
&CMakeLists.txt &
&Scripts/tag_version.sh &
&hdf5_tools_module.f90 &
&neo2.f90 &
&"
