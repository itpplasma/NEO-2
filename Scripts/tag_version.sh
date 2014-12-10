#!/bin/bash
GITHASH=`git rev-parse HEAD`
GITSHORTHASH=`git rev-parse --short HEAD`
GITCHANGEDFILES=`git diff-index --name-only HEAD`

export GITSHORTHASH=$GITSHORTHASH

echo "character(len=*), parameter :: CMake_Compiler = '@CMAKE_Fortran_COMPILER_ID@'" > ./version.f90.in
echo "character(len=*), parameter :: CMake_Compiler_Version = '@GCC_MAJOR@.@GCC_MINOR@'" >> ./version.f90.in
echo "character(len=*), parameter :: CMake_Build_Type = '@CMAKE_BUILD_TYPE@'" >> ./version.f90.in
echo "character(len=*), parameter :: CMake_Flags = '@CMAKE_Fortran_FLAGS@'" >> ./version.f90.in
echo "character(len=*), parameter :: CMake_Flags_Release = '@CMAKE_Fortran_FLAGS_RELEASE@'" >> ./version.f90.in
echo "character(len=*), parameter :: CMake_Flags_Debug = '@CMAKE_Fortran_FLAGS_DEBUG@'" >> ./version.f90.in
echo "character(len=*), parameter :: CMake_System = '@CMAKE_SYSTEM@'" >> ./version.f90.in
echo "character(len=*), parameter :: CMake_SuiteSparse_Dir = '@SUITESPARSE_DIR@'" >> ./version.f90.in
echo "character(len=*), parameter :: CMake_Blas_Lib = '@BLAS_lib@'" >> ./version.f90.in

echo "character(len=*), parameter :: Neo2_Version = '${GITHASH}'" >> ./version.f90.in
echo ""
echo "Git Hash is ${GITHASH}"

if [ -n "$GITCHANGEDFILES" ]; then
    echo 'character(len=*), parameter :: Neo2_Version_Additional = "WARNING, &' >> ./version.f90.in
    echo "&THERE ARE UNCOMMITTED CHANGES. Run may not be reproduceable: &" >> ./version.f90.in
    echo "THERE ARE UNCOMMITTED CHANGES!"
    echo ""
    while read -r line; do
        echo "&${line} &" >> ./version.f90.in
    done <<< "$GITCHANGEDFILES"
    echo '&"' >> ./version.f90.in

else

    echo 'character(len=*), parameter :: Neo2_Version_Additional = ""' >> ./version.f90.in

fi
