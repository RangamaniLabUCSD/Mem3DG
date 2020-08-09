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
if(DOXYGEN_FOUND)
    # Parse Doxyfile.in and replace CMake @macros@.
    configure_file(${CMAKE_SOURCE_DIR}/docs/Doxyfile.in ${CMAKE_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(mem3dg_docs
        ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
else()
    message(NOTICE "CMake could not find Doxygen. Please install Doxygen or help me find your Doxygen binary to build the documentation!")
endif(DOXYGEN_FOUND)
