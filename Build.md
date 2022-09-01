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
  * SuiteSparse, SuperLU
  * ATLAS, OpenBLAS
  * GSL, FGSL
  * MPI
    Used for parallelization.
  * HDF5
    Used for output.
  * libneo
    Contains common to different codes of the group. E.g. interface to
    rusage.
  * mc_efit
    Part of libneo repository. Neo-2 uses magfie code of the library
  * hdf5_tools, MyMPILib
    Part of libneo repository.Fortran interfaces to the respective
    library.

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
