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

cmake_minimum_required(VERSION 3.11)
cmake_policy(SET CMP0063 NEW)

# Set default build type
set(CMAKE_BUILD_TYPE_INIT Release)
# Add path to custom modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

include(DefaultBuildTypeAndFlags)
include(Macros)

get_version_from_git()
project(
  Mem3DG
  VERSION ${VERSION_SHORT}
  LANGUAGES CXX
)
message(STATUS "Mem3DG version: ${VERSION}")

# Configure version.cpp to give access to version in code
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.cpp.in
  ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
)
set(VERSION_FILE "${CMAKE_CURRENT_BINARY_DIR}/version.cpp")

# Require c++14 and standard libraries
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Initialize output directory locations
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# Check if DDG is being used directly or via add_subdirectory
set(DDG_MASTER_PROJECT OFF)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(DDG_MASTER_PROJECT ON)
endif()

print_debug_messages()

include(FetchContent)

option(BUILD_PYMEM3DG "Build the python extensions?" ON)
option(WITH_NETCDF "Build with NetCDF (binary trajectory output)?" OFF)

# ##############################################################################
# BUNDLED LIBRARIES & EXTERNAL LIBS
# ##############################################################################
add_subdirectory(libraries)

set(LINKED_LIBS geometry-central pcg::pcg)

if(WITH_NETCDF)
  find_package(NetCDF_cxx4 REQUIRED)
  list(APPEND LINKED_LIBS NetCDF::NetCDF_cxx4)
endif()

# ##############################################################################
# DDG SOLVER LIBRARY
# ##############################################################################
add_subdirectory(include) # Get list of headers
add_subdirectory(src) # Get list of sources
list(APPEND DDG_SOURCES ${VERSION_FILE})

# ddgsolver object library
add_library(ddgsolver_objlib OBJECT ${DDG_SOURCES} ${DDG_HEADERS})
target_include_directories(
  ddgsolver_objlib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                          $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(ddgsolver_objlib PUBLIC ${LINKED_LIBS})
set_target_properties(
  ddgsolver_objlib PROPERTIES POSIION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET
                                                          hidden
) # VISIBILITY_INLINES_HIDDEN ON)
if(WITH_NETCDF)
  target_compile_definitions(ddgsolver_objlib PUBLIC -DMEM3DG_WITH_NETCDF)
endif()

# ddgsolver library
add_library(ddgsolver SHARED $<TARGET_OBJECTS:ddgsolver_objlib>)
target_link_libraries(ddgsolver PUBLIC ddgsolver_objlib)

if(BUILD_PYMEM3DG)
  add_subdirectory(pymem3dg)
endif()

# ##############################################################################
# Main DDG target
# ##############################################################################
# add_executable(ddg tests/scratch.cpp) target_link_libraries(ddg PRIVATE
# ddgsolver polyscope) set_target_properties(ddg PROPERTIES
# CXX_VISIBILITY_PRESET hidden)


enable_testing()
add_subdirectory(tests)

include(GNUInstallDirs)
# Install targets and headers
install(
  TARGETS ddgsolver
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(
  DIRECTORY include/ddgsolver
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING
  PATTERN "*.h"
)
