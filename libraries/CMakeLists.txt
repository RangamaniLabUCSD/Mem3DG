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

include(FetchContent)

add_subdirectory(pcg)

# ##############################################################################
# GET EIGEN
# ##############################################################################
FetchContent_Declare(
  eigen
  GIT_REPOSITORY  https://gitlab.com/libeigen/eigen.git
  GIT_TAG         3.3.7
  GIT_SHALLOW     TRUE
  SOURCE_DIR      "${CMAKE_CURRENT_BINARY_DIR}/eigen-src"
  BINARY_DIR      "${CMAKE_CURRENT_BINARY_DIR}/eigen-build"
)
FetchContent_GetProperties(eigen)
if(NOT eigen_POPULATED)
    FetchContent_Populate(eigen)
    include(EigenChecker)
    eigen3checker(${eigen_SOURCE_DIR} 3.3)
endif()
add_library (eigen INTERFACE)
add_library (Eigen3::Eigen ALIAS eigen)
target_include_directories(eigen INTERFACE
    $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>
    )

# ##############################################################################
# GET GEOMETRY CENTRAL
# ##############################################################################
FetchContent_Declare(
  geometrycentrallib
  GIT_REPOSITORY https://github.com/nmwsharp/geometry-central.git
  GIT_TAG origin/master
  GIT_SHALLOW TRUE
  SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/geometrycentral-src" 
  BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/geometrycentral-build"
)
FetchContent_GetProperties(geometrycentrallib)
if(NOT geometrycentrallib_POPULATED)
  FetchContent_Populate(geometrycentrallib)
  add_subdirectory(
    ${geometrycentrallib_SOURCE_DIR} ${geometrycentrallib_BINARY_DIR}
  )
  set_target_properties(
    geometry-central PROPERTIES POSITION_INDEPENDENT_CODE ON
  )
endif()


if(BUILD_PYMEM3DG)
  # ############################################################################
  # GET polyscope
  # ############################################################################
  FetchContent_Declare(
    polyscopelib
    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
    GIT_SHALLOW TRUE
    SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/polyscope-src" 
    BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/polyscope-build"
  )
  FetchContent_GetProperties(polyscopelib)
  if(NOT polyscopelib_POPULATED)
    FetchContent_Populate(polyscopelib)
    add_subdirectory(${polyscopelib_SOURCE_DIR} ${polyscopelib_BINARY_DIR})
    set_target_properties(polyscope PROPERTIES POSITION_INDEPENDENT_CODE ON)
  endif()
endif()


# ##############################################################################
# GET netcdf-cxx4
# ##############################################################################
# FetchContent_Declare(
#   netcdf-cxx4
#   GIT_REPOSITORY https://github.com/ctlee/netcdf-cxx4.git
#   GIT_SHALLOW TRUE
#   SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/netcdf-cxx4-src" 
#   BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/netcdf-cxx4-build"
# )
# FetchContent_GetProperties(netcdf-cxx4)
# if(NOT netcdf-cxx4_POPULATED)
#   FetchContent_Populate(netcdf-cxx4)
#   add_subdirectory(${netcdf-cxx4_SOURCE_DIR} ${netcdf-cxx4_BINARY_DIR})
#   set_target_properties(netcdf-cxx4 PROPERTIES POSITION_INDEPENDENT_CODE ON)
#   if(NOT TARGET NetCDF::NetCDF-cxx4)
#     add_library(NetCDF::NetCDF-cxx4 UNKNOWN IMPORTED)
#     set_target_properties(
#       NetCDF::NetCDF-cxx4 PROPERTIES INTERFACE_LINK_LIBRARIES "netcdf-cxx4"
#     )
#   endif()
# endif()
