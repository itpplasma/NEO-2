==== How to build the neo2 code ====

This documents describes the basic steps to get the neo2 executable(s).
It does not go into details on specific machines/clusters.

The build system used by neo-2 is cmake.

=== Dependencies ===
Here is a list of libraries and tools required for building/running neo2.
  * git
    Used for version control.
  * cmake
    Used for configuration, has so far only been used in connection with
    GNUMake.
  * SuiteSparse
  * Either OpenBLAS, MKL, or ATLAS for LAPACK/BLAS
  * fortnum
    Numerical backend (integration, special functions, B-splines). Fetched
    together with libneo; replaces the former GSL/FGSL dependency.
  * MPI
    Used for parallelization.
  * libneo
    Contains common to different codes of the group. E.g. interface to
    rusage. Its Fortio-backed compatibility modules provide the HDF5 and
    NetCDF subset used by NEO-2, so system HDF5/NetCDF libraries are not
    required.
  * mc_efit
    Part of libneo repository. Neo-2 uses magfie code of the library
  * hdf5_tools
    Part of libneo. Stable Fortran compatibility interface backed by Fortio.
  * MyMPILib
    Part of the NEO-2 repository. Fortran interface to MPI used by the
    parallel code paths.

Note, that as compiler, so far mainly gfortran has been used. The intel
compiler should work, but this is so far not certain.

=== Configuration ===
Configuration, i.e. setting of library paths, is done with the file
ProjectConfig.cmake.in (should be made configurable) in the COMMON/
folder.
Setting of the paths for the required libraries, should be done in this
file.

=== Building ===
This has to be done for both versions of the code seperately, and also
for the documentation.

== Documentation ==
Lets start with the latter.
Change into the DOC folder, and create a directory called 'Build', if it
not already exists. Then go into the folder. The first time the
directory is created, you have to execute 'cmake ..'. This will run
cmake, and tells it to look for its configuration file ('CMakeLists.txt')
in the top directory. Cmake will then to the configuration, i.e. it will
check for a compiler, linker, libraries and maybe more and store the
results. It will also create a 'Makefile'. Once you have this, there is
normally no need to run 'cmake ..' again, as it will be called
automatically, if required.
If you have the Makefile, it can be run with 'make'. This will now
actually create the documenation. You may see the different tasks
printed to the screen, along with and indication of how much is already
done.
Once the process is done, you can find the user manual as 'main.pdf' in
the Build folder. For developers of interest is the doxygen
documentation which can be found in the 'html' and 'latex' subfolders.
For the latter it is 'refman.pdf', while the starting point for the
html version is 'index.html'.

== Coding Style ==
The coding style is an extra document in the DOC folder, which is not
made by default. If a cmake build directory has been set up, you can
create it by 'make coding_style_safepdf'.

== Neo-2 ==
To build either of the code variants, you first have to configure the
code (see previous section).
If this is done, change from the main directory of neo2, into the
appropriate subdirectory (i.e. 'NEO-2-PAR' or 'NEO-2-QL') and create a
directory Build-Debug (other ones, like Build-Release, might be created
later if needed).
Change into this directory and run 'cmake ..', to start the
configuration. As for the documentation, this has to be done only once
(normally), as the make command will also execute this if necessary.
When the configuration is done, there should be a 'Makefile' which can
be used by executing 'make'.
Compiling and linking of the code may take some time, but once it is
finished, the result will be an executable named 'neo_2.x'.

There are five possibilities to use the executable for your simulations.
  * Do simulation in the Build-folder (not recommended).
  * Copy into the folder where you want to run the simulation.
  * Make a link to the executable, where you want to run the simulation.
  * Set your PATH variable to include the build folder (not recommended).
  * Copy the executable into a folder to which the PATH variable points
  (not really recommended, as both code variants have the same
  executable name).
Note that some tools might have troubles with some of the methods, which
can put further restrictions on your posibilities.

With this the build process is done, and you can now do simulations.

== Octave/Python/Shell scripts ==
To be able to use these, the coresponding directory probably needs to be
added to the coresponding path variable, so that the scripts can be
found.

== CINECA Leonardo ==
Built and tested on the Leonardo CPU login nodes (GCC 12.2.0, CMake 3.27.9,
Ninja 1.11.1, OpenBLAS 0.3.26, OpenMPI 4.1.6, GSL 2.7.1) at NEO-2 117b7ea,
2026-10-19; the NEO-2 test suite (below) passed 6 of 6, also as a batch job on
the budget-free serial partition (lrd_all_serial, one task, four cores).
Pitagora access was not available, so nothing is verified there. GitHub is
reachable from the login nodes, so CMake fetches libneo, Fortnum and the
other dependencies as usual. No golden-record or example run was done.

Leonardo has no SuiteSparse module. The automatic fallback download of
SuiteSparse does not work: it needs MPFR >= 4.0.2 (the system has 3.1.6),
and COMMON/CMakeLists.txt links the target SuiteSparse::umfpack, which the
in-tree build does not define. Build SuiteSparse once into a prefix and
point NEO-2's FindSuiteSparse at it. KLU, BTF and SPQR are built because
cmake/FindSuiteSparse.cmake requires them. Do not put long -I or -L paths
into CMAKE_Fortran_FLAGS: the flags are written into a generated Fortran
line that is truncated beyond 132 columns.

----
module purge
module load profile/base gcc/12.2.0 cmake/3.27.9 ninja/1.11.1 \
    openblas/0.3.26--gcc--12.2.0 openmpi/4.1.6--gcc--12.2.0-cuda-12.2 \
    gsl/2.7.1--gcc--12.2.0-spack0.22
export B=$CINECA_SCRATCH/builds TMPDIR=$CINECA_SCRATCH/builds/tmp
mkdir -p $TMPDIR; cd $B
export CMAKE_BUILD_PARALLEL_LEVEL=4

# SuiteSparse 7.6.1 (UMFPACK and the components NEO-2's FindSuiteSparse needs)
git clone -q --depth 1 -b v7.6.1 https://github.com/DrTimothyAldenDavis/SuiteSparse.git
nice cmake -S SuiteSparse -B ss-build -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$B/suitesparse -DSUITESPARSE_USE_CUDA=OFF \
    -DSUITESPARSE_USE_FORTRAN=OFF -DBUILD_STATIC_LIBS=OFF \
    "-DSUITESPARSE_ENABLE_PROJECTS=suitesparse_config;amd;btf;camd;ccolamd;colamd;cholmod;klu;spqr;umfpack"
nice cmake --build ss-build -j 4 && cmake --install ss-build

git clone -q https://github.com/itpplasma/NEO-2.git && cd NEO-2
SS=$B/suitesparse
nice cmake -S . -B build-lrd -G Ninja -DCMAKE_BUILD_TYPE=Release \
    "-DCMAKE_EXE_LINKER_FLAGS=-Wl,-rpath,$SS/lib64" \
    "-DSUITESPARSE_INCLUDE_DIR_HINTS=$SS/include/suitesparse" \
    "-DSUITESPARSE_LIBRARY_DIR_HINTS=$SS/lib64"
nice cmake --build build-lrd -j 4
ctest --test-dir build-lrd/TEST --output-on-failure
----

The SuiteSparse step takes about 1 minute and NEO-2 (with the fetched
dependencies) about 2 minutes with -j 4 under nice. The build produces
build-lrd/NEO-2-QL/neo_2_ql.x and build-lrd/NEO-2-PAR/neo_2_par.x. ctest runs
boozmn_read_test, nrutil_test, er_rotation_test, spline_sparse_assembly_test,
arnoldi_breakdown_test and test_mympilib (all passed; run longer cases in
a batch job, not on the login node).
