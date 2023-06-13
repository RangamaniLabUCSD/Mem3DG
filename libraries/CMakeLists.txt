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

if(M3DG_GET_OWN_EIGEN)
  # ############################################################################
  # GET EIGEN
  # ############################################################################
  FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE
  )
  FetchContent_GetProperties(eigen)
  if(NOT eigen_POPULATED)
    FetchContent_Populate(eigen)
    include(EigenChecker)
    eigen3checker(${eigen_SOURCE_DIR} 3.3)
  endif()
  add_library(eigen INTERFACE)
  add_library(Eigen3::Eigen ALIAS eigen)
  target_include_directories(
    eigen SYSTEM INTERFACE $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIR}>
                           $<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>
  )
else()
  find_package(Eigen3 3.3 REQUIRED NO_MODULE)
endif()

# ##############################################################################
# GET GEOMETRY CENTRAL
# ##############################################################################
# FetchContent_Declare( geometrycentrallib GIT_REPOSITORY
# https://github.com/nmwsharp/geometry-central.git GIT_TAG origin/master
# GIT_SHALLOW TRUE SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/geometrycentral-src"
# BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/geometrycentral-build" )
# FetchContent_GetProperties(geometrycentrallib) if(NOT
# geometrycentrallib_POPULATED) FetchContent_Populate(geometrycentrallib)
# add_subdirectory( ${geometrycentrallib_SOURCE_DIR}
# ${geometrycentrallib_BINARY_DIR} ) set_target_properties( geometry-central
# PROPERTIES POSITION_INDEPENDENT_CODE ON ) endif()
add_subdirectory(geometry-central SYSTEM EXCLUDE_FROM_ALL)
set_target_properties(geometry-central PROPERTIES POSITION_INDEPENDENT_CODE ON)

# ##############################################################################
# Configure libigl
# ##############################################################################
add_subdirectory(libigl EXCLUDE_FROM_ALL SYSTEM)
