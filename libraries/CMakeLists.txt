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
#     Padmini Rangmani (prangamani@eng.ucsd.edu)
# 

include(FetchContent)

set(CMAKE_POLICY_DEFAULT_CMP0063 NEW)

add_subdirectory(pcg)

# ##############################################################################
# GET GEOMETRY CENTRAL
# ##############################################################################
FetchContent_Declare(
  geometrycentrallib
  GIT_REPOSITORY https://github.com/ctlee/geometry-central.git
  GIT_TAG origin/raw_access
  GIT_SHALLOW TRUE
  SOURCE_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/geometrycentral-src"
  BINARY_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/geometrycentral-build"
)
FetchContent_GetProperties(geometrycentrallib)
if(NOT geometrycentrallib_POPULATED)
  FetchContent_Populate(geometrycentrallib)
  add_subdirectory(
    ${geometrycentrallib_SOURCE_DIR} ${geometrycentrallib_BINARY_DIR}
  )
  set_target_properties(geometry-central PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden)
endif()

# ##############################################################################
# GET polyscope
# ##############################################################################
FetchContent_Declare(
  polyscopelib
  GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
  GIT_SHALLOW TRUE
  SOURCE_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/polyscope-src"
  BINARY_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/polyscope-build"
)
FetchContent_GetProperties(polyscopelib)
if(NOT polyscopelib_POPULATED)
  FetchContent_Populate(polyscopelib)
  add_subdirectory(${polyscopelib_SOURCE_DIR} ${polyscopelib_BINARY_DIR})
  set_target_properties(polyscope PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden)
endif()
