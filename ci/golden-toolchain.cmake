include("${CMAKE_CURRENT_LIST_DIR}/cluster-toolchain.cmake")
# Preserve the golden workflow's original project compiler flags. The unit
# tier's initial flags would change the checked legacy CMake_Flags metadata.
set(CMAKE_Fortran_FLAGS_INIT "")
# The stable release and older main revisions build libneo revisions that pin
# fortnum 1a52e13, which the lazy-fortran history rewrite of 14 November 2026
# removed; 79d94dd is the same commit after the rewrite. Those libneo revisions
# keep a FORTNUM_REF that is already set; libneo 94f59c1 and later ignore it.
set(FORTNUM_REF 79d94dd2fe54048e177f3a1e58a9d76f06dad0b7 CACHE STRING
    "fortnum revision for libneo revisions before 94f59c1")
