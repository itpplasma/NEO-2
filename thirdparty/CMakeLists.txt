include(ExternalProject)

include(libneo.cmake)
include(fgsl.cmake)

add_dependencies(LIBNEO FGSL)
