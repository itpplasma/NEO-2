#cmake -DLATEX_DEFAULT_BUILD=SafePdf ..

cmake_minimum_required (VERSION 3.0)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../CmakeScripts/")

set(LATEX_DEFAULT_BUILD "SafePdf")

project(neo-2-documentation)

enable_language (Fortran)

find_package(LATEX)
find_package(Doxygen
             REQUIRED dot
             OPTIONAL_COMPONENTS mscgen dia)

add_custom_command(
  OUTPUT ../Doxygen_Referenz_Manual.pdf
  COMMAND cmake -E chdir ${CMAKE_SOURCE_DIR}/Build/latex/ make && cp ${CMAKE_SOURCE_DIR}/Build/latex/refman.pdf ${CMAKE_SOURCE_DIR}/Doxygen_Referenz_Manual.pdf
  DEPENDS latex
  )

add_custom_command(
  OUTPUT latex html
  COMMAND cmake -E env ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/Doxyfile
  DEPENDS ${CMAKE_SOURCE_DIR}/../COMMON/*.f90 ${CMAKE_SOURCE_DIR}/../NEO-2-PAR/*.f90 ${CMAKE_SOURCE_DIR}/../NEO-2-QL/*.f90
  )

add_custom_target(documentation ALL DEPENDS ${CMAKE_SOURCE_DIR}/Doxygen_Referenz_Manual.pdf)

# Include script that allows easy compilation of LaTeX files with cmake.
# Not from us, but added to the repository for easy acess.
include(UseLATEX)

# Command for the main documentation.
add_latex_document(UserDoc/main.tex
  INPUTS
  UserDoc/Commands.sty
  UserDoc/Bugs.tex
  UserDoc/Building.tex
  UserDoc/Input.tex
  UserDoc/Introduction.tex
  UserDoc/Output.tex
  UserDoc/Programflow.tex
  UserDoc/Running.tex
  UserDoc/Todo.tex
)

# Command for the coding style document. This is not done by default.
add_latex_document(ExtraDocumentsTex/coding_style.tex
  EXCLUDE_FROM_ALL
)
    #BIBFILES mybib.bib
    #IMAGE_DIRS images
