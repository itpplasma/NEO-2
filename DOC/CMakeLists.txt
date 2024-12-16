cmake_minimum_required (VERSION 3.20)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CmakeScripts/")

set(LATEX_DEFAULT_BUILD "SafePdf")

project(neo-2-documentation)

find_package(LATEX)
find_package(Doxygen
             REQUIRED dot
             OPTIONAL_COMPONENTS mscgen dia)


add_custom_command(
  OUTPUT ${CMAKE_BINARY_DIR}/Doxygen_Reference_Manual.pdf
  COMMAND cmake -E chdir ${PROJECT_BINARY_DIR}/latex/ make && cp ${PROJECT_BINARY_DIR}/latex/refman.pdf ${CMAKE_BINARY_DIR}/Doxygen_Reference_Manual.pdf
  DEPENDS latex
  )

file(GLOB COMMON_SOURCES ${CMAKE_SOURCE_DIR}/COMMON/*.f90)
file(GLOB NEO_2_PAR_SOURCES ${CMAKE_SOURCE_DIR}/NEO-2-PAR/*.f90)
file(GLOB NEO_2_QL_SOURCES ${CMAKE_SOURCE_DIR}/NEO-2-QL/*.f90)

add_custom_command(
  OUTPUT latex
  COMMAND cmake -E env ${DOXYGEN_EXECUTABLE} ${PROJECT_SOURCE_DIR}/Doxyfile
  DEPENDS ${COMMON_SOURCES} ${NEO_2_PAR_SOURCES} ${NEO_2_QL_SOURCES}
  )

add_custom_target(doc ALL DEPENDS ${CMAKE_BINARY_DIR}/Doxygen_Reference_Manual.pdf)

# Include script that allows easy compilation of LaTeX files with cmake.
# Not from us, but added to the repository for easy acess.
include(UseLATEX)

# Command for the main documentation.
add_latex_document(UserDoc/main.tex
  INPUTS
  UserDoc/Bugs.tex
  UserDoc/Building.tex
  UserDoc/Input.tex
  UserDoc/Introduction.tex
  UserDoc/Output.tex
  UserDoc/Programflow.tex
  UserDoc/Running.tex
  UserDoc/Theory.tex
  UserDoc/Todo.tex
)
