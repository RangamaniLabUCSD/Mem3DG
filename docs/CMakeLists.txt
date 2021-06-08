# Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)
# 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
# 
# Copyright (c) 2020:
#     Laboratory for Computational Cellular Mechanobiology
#     Cuncheng Zhu (cuzhu@eng.ucsd.edu)
#     Christopher T. Lee (ctlee@ucsd.edu)
#     Ravi Ramamoorthi (ravir@cs.ucsd.edu)
#     Padmini Rangamani (prangamani@eng.ucsd.edu)
# 


#####################################################################
# Target to generate Doxygen documentation
#####################################################################
find_package(Doxygen)
set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/cpp")
set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "README.md")
set(DOXYGEN_PROJECT_NUMBER ${VERSION})
doxygen_add_docs(mem3dg_doxy_docs
    ${CMAKE_SOURCE_DIR}/README.md
    ${DDG_HEADERS}
)

# if(DOXYGEN_FOUND)
#     # Parse Doxyfile.in and replace CMake @macros@.
#     configure_file(${CMAKE_SOURCE_DIR}/docs/Doxyfile.in ${CMAKE_BINARY_DIR}/Doxyfile @ONLY)
#     add_custom_target(mem3dg_doxy_docs
#         ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile
#         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#         COMMENT "Generating API documentation with Doxygen" VERBATIM
#     )
# else()
#     message(NOTICE "CMake could not find Doxygen. Please install Doxygen or help me find your Doxygen binary to build the documentation!")
# endif(DOXYGEN_FOUND)

find_package(Sphinx)
SET(SPHINXOPTS CACHE STRING "Extra options to pass to Sphinx")

if(SPHINX_FOUND)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
        "${CMAKE_CURRENT_BINARY_DIR}/conf.py"
        @ONLY)

    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/conf.py"
        @ONLY)

    add_custom_target(mem3dg_sphinx_docs
        ${SPHINX_EXECUTABLE}
            -M html
            "${CMAKE_CURRENT_SOURCE_DIR}/src"
            "${CMAKE_CURRENT_BINARY_DIR}/py"
            -c "${CMAKE_CURRENT_BINARY_DIR}"
            ${SPHINXOPTS}
        # DEPENDS pymem3dg
        COMMENT "Building HTML documentation with Sphinx")
    message(STATUS "Build the Sphinx documentation using `make mem3dg_sphinx_docs`")
else()
    message(STATUS "CMake could not find Sphinx. Please install Sphinx or help me find your Sphinx binary to build the documentation!")
endif(SPHINX_FOUND)

set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/src/_pythonapi"
        "${CMAKE_CURRENT_BINARY_DIR}/doctrees"
        "${CMAKE_CURRENT_BINARY_DIR}/html"
)
